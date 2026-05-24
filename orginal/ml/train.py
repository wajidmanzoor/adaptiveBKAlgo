#!/usr/bin/env python3
"""
train.py — train XGBoost and Random Forest classifiers for pruning-rule selection.

Usage:
    python train.py <dataset_csv> [test_size]
    test_size defaults to 0.2

Expected CSV columns:
    graph                    (dropped — identifier only)
    <feature columns>        (all structural graph features)
    prune1 prune2 sp1 sp2 sp3 sp4 sp5 sp6   (binary labels)

Outputs:
    ml/weights/best_model.pkl    best model (sklearn pipeline)
    ml/weights/best_model_info.txt  human-readable summary
    ml/plots/*.png               accuracy, ROC, feature importance, confusion
"""

import sys
import pickle
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.ensemble import (
    RandomForestClassifier, ExtraTreesClassifier,
    GradientBoostingClassifier,
)
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.multioutput import MultiOutputClassifier
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import (
    accuracy_score, hamming_loss, f1_score,
    roc_auc_score, roc_curve, confusion_matrix,
)
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

try:
    from xgboost import XGBClassifier
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
    print("Warning: xgboost not installed — skipping XGBoost.")

warnings.filterwarnings("ignore")

# ── Constants ─────────────────────────────────────────────────────────────────

LABEL_COLS = ["prune1", "prune2", "sp1", "sp2", "sp3", "sp4", "sp5", "sp6"]
SKIP_COLS  = ["graph"]

# ── Paths ─────────────────────────────────────────────────────────────────────

SCRIPT_DIR  = Path(__file__).parent
WEIGHTS_DIR = SCRIPT_DIR / "weights"
PLOTS_DIR   = SCRIPT_DIR / "plots"
WEIGHTS_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# ── Data loading ──────────────────────────────────────────────────────────────

def load_data(csv_path):
    df = pd.read_csv(csv_path)

    missing = [c for c in LABEL_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"CSV is missing label columns: {missing}")

    feat_cols = [c for c in df.columns if c not in LABEL_COLS + SKIP_COLS]
    if not feat_cols:
        raise ValueError("No feature columns found in CSV.")

    X = df[feat_cols].astype(float).values
    y = df[LABEL_COLS].astype(int).values

    print(f"Dataset   : {len(df)} graphs,  {len(feat_cols)} features,  {len(LABEL_COLS)} labels")
    print(f"Features  : {feat_cols}")
    print(f"Label prevalence (fraction ON per rule):")
    for i, col in enumerate(LABEL_COLS):
        print(f"  {col:<10} {y[:, i].mean():.2%}")
    print()

    return X, y, feat_cols

# ── Scoring helpers ───────────────────────────────────────────────────────────

def hamming_acc(y_true, y_pred):
    """Per-label accuracy averaged across all labels."""
    return 1.0 - hamming_loss(y_true, y_pred)

def exact_match(y_true, y_pred):
    """All 8 labels correct for a sample."""
    return accuracy_score(y_true, y_pred)

def per_label_f1(y_true, y_pred):
    return f1_score(y_true, y_pred, average=None, zero_division=0)

def per_label_auc(y_true, y_prob_list):
    """y_prob_list: list of (n_samples, 2) arrays from MultiOutputClassifier."""
    aucs = []
    for i, probs in enumerate(y_prob_list):
        col = y_true[:, i]
        if len(np.unique(col)) < 2:
            aucs.append(float("nan"))
        else:
            aucs.append(roc_auc_score(col, probs[:, 1]))
    return np.array(aucs)

# ── Grid search configs ───────────────────────────────────────────────────────

RF_PARAM_GRID = {
    "estimator__n_estimators":      [100, 300],
    "estimator__max_depth":         [None, 10, 20],
    "estimator__min_samples_split": [2, 5],
    "estimator__min_samples_leaf":  [1, 2],
}

ET_PARAM_GRID = {
    "estimator__n_estimators":      [100, 300],
    "estimator__max_depth":         [None, 10, 20],
    "estimator__min_samples_split": [2, 5],
    "estimator__min_samples_leaf":  [1, 2],
}

GB_PARAM_GRID = {
    "estimator__n_estimators":  [100, 200],
    "estimator__max_depth":     [3, 5],
    "estimator__learning_rate": [0.05, 0.1, 0.2],
    "estimator__subsample":     [0.8, 1.0],
}

XGB_PARAM_GRID = {
    "estimator__n_estimators":  [100, 300],
    "estimator__max_depth":     [3, 6],
    "estimator__learning_rate": [0.05, 0.1, 0.2],
    "estimator__subsample":     [0.8, 1.0],
}

DT_PARAM_GRID = {
    "estimator__max_depth":         [None, 5, 10, 15],
    "estimator__min_samples_split": [2, 5, 10],
    "estimator__min_samples_leaf":  [1, 2, 4],
    "estimator__criterion":         ["gini", "entropy"],
}

LR_PARAM_GRID = {
    "estimator__C":        [0.01, 0.1, 1, 10],
    "estimator__solver":   ["lbfgs", "liblinear"],
    "estimator__max_iter": [2000],
}

SVM_PARAM_GRID = {
    "estimator__C":      [0.1, 1, 10],
    "estimator__kernel": ["rbf", "linear"],
    "estimator__gamma":  ["scale", "auto"],
}

KNN_PARAM_GRID = {
    "estimator__n_neighbors": [3, 5, 7, 9],
    "estimator__weights":     ["uniform", "distance"],
    "estimator__metric":      ["euclidean", "manhattan"],
}

NB_PARAM_GRID = {
    "estimator__var_smoothing": [1e-9, 1e-7, 1e-5],
}

# ── Model training ────────────────────────────────────────────────────────────

def build_pipeline(base_clf):
    return Pipeline([
        ("scaler", StandardScaler()),
        ("clf",    MultiOutputClassifier(base_clf, n_jobs=-1)),
    ])

def prefix_grid(grid, pipeline_step="clf"):
    """Prefix param grid keys for use inside a Pipeline."""
    return {f"{pipeline_step}__{k}": v for k, v in grid.items()}

def run_grid_search(pipeline, param_grid, X_train, y_train, cv=5):
    from sklearn.metrics import make_scorer
    scorer = make_scorer(
        lambda yt, yp: 1.0 - hamming_loss(yt, yp),
        greater_is_better=True,
    )
    gs = GridSearchCV(
        pipeline,
        param_grid,
        scoring=scorer,
        cv=cv,
        n_jobs=-1,
        verbose=1,
        refit=True,
    )
    gs.fit(X_train, y_train)
    return gs

def evaluate(model, X_test, y_test, name):
    y_pred       = model.predict(X_test)
    y_prob_list  = model.predict_proba(X_test)   # list of (n, 2) per label

    ham_acc  = hamming_acc(y_test, y_pred)
    ex_acc   = exact_match(y_test, y_pred)
    f1_per   = per_label_f1(y_test, y_pred)
    auc_per  = per_label_auc(y_test, y_prob_list)

    print(f"── {name} ──────────────────────────────")
    print(f"  Hamming accuracy  : {ham_acc:.4f}")
    print(f"  Exact-match acc   : {ex_acc:.4f}")
    print(f"  Per-label F1:")
    for col, f1, auc in zip(LABEL_COLS, f1_per, auc_per):
        auc_str = f"{auc:.3f}" if not np.isnan(auc) else "  N/A"
        print(f"    {col:<10}  F1={f1:.3f}  AUC={auc_str}")
    print()

    return {
        "hamming_acc": ham_acc,
        "exact_acc":   ex_acc,
        "f1_per":      f1_per,
        "auc_per":     auc_per,
        "y_pred":      y_pred,
        "y_prob_list": y_prob_list,
    }

# ── Feature importance ────────────────────────────────────────────────────────

def get_importances(model, feat_cols):
    """Average feature importance across all label-specific estimators."""
    clf = model.named_steps["clf"]
    imps = np.array([e.feature_importances_ for e in clf.estimators_])
    return imps.mean(axis=0)

# ── Plots ─────────────────────────────────────────────────────────────────────

COLORS = {
    "Random Forest":    "#2196F3",
    "Extra Trees":      "#4CAF50",
    "Gradient Boost":   "#9C27B0",
    "XGBoost":          "#FF5722",
    "Decision Tree":    "#FF9800",
    "Logistic Reg":     "#607D8B",
    "SVM":              "#E91E63",
    "KNN":              "#00BCD4",
    "Naive Bayes":      "#795548",
}

def plot_model_comparison(results, test_size):
    """Bar chart: Hamming acc and Exact-match acc for each model."""
    models  = list(results.keys())
    ham_acc = [results[m]["hamming_acc"] for m in models]
    ex_acc  = [results[m]["exact_acc"]   for m in models]

    x      = np.arange(len(models))
    width  = 0.35
    fig, ax = plt.subplots(figsize=(7, 5))

    bars1 = ax.bar(x - width/2, ham_acc, width, label="Hamming Accuracy",
                   color=[COLORS.get(m, "#888") for m in models], alpha=0.9)
    bars2 = ax.bar(x + width/2, ex_acc,  width, label="Exact-Match Accuracy",
                   color=[COLORS.get(m, "#888") for m in models], alpha=0.5,
                   edgecolor="black")

    ax.set_ylim(0, 1.1)
    ax.set_xticks(x)
    ax.set_xticklabels(models, fontsize=12)
    ax.set_ylabel("Accuracy")
    ax.set_title("Model Comparison (test set)")
    ax.legend()
    for bar in list(bars1) + list(bars2):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                f"{bar.get_height():.3f}", ha="center", va="bottom", fontsize=9)

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "model_comparison.png", dpi=150)
    plt.close()
    print("Saved: model_comparison.png")


def plot_per_rule_metrics(results):
    """Grouped bar: F1 per rule for each model."""
    models = list(results.keys())
    n_labels = len(LABEL_COLS)
    x = np.arange(n_labels)
    width = 0.35 / max(len(models) - 1, 1) * 2

    fig, ax = plt.subplots(figsize=(12, 5))
    offsets = np.linspace(-width * (len(models)-1)/2,
                           width * (len(models)-1)/2, len(models))

    for offset, m in zip(offsets, models):
        ax.bar(x + offset, results[m]["f1_per"], width,
               label=m, color=COLORS.get(m, "#888"), alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(LABEL_COLS, fontsize=11)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("F1 Score")
    ax.set_title("Per-Rule F1 Score (test set)")
    ax.legend()
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "per_rule_f1.png", dpi=150)
    plt.close()
    print("Saved: per_rule_f1.png")


def plot_roc_curves(results, y_test):
    """ROC curve per label for each model."""
    n_labels = len(LABEL_COLS)
    cols     = 4
    rows     = (n_labels + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 3.5))
    axes = axes.flatten()

    for i, col in enumerate(LABEL_COLS):
        ax = axes[i]
        y_col = y_test[:, i]

        if len(np.unique(y_col)) < 2:
            ax.text(0.5, 0.5, "Only one class", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(col)
            continue

        for m_name, res in results.items():
            probs = res["y_prob_list"][i][:, 1]
            fpr, tpr, _ = roc_curve(y_col, probs)
            auc = res["auc_per"][i]
            ax.plot(fpr, tpr, label=f"{m_name} (AUC={auc:.2f})",
                    color=COLORS.get(m_name, "#888"), linewidth=1.8)

        ax.plot([0, 1], [0, 1], "k--", linewidth=0.8)
        ax.set_title(col, fontsize=11)
        ax.set_xlabel("FPR", fontsize=9)
        ax.set_ylabel("TPR", fontsize=9)
        ax.legend(fontsize=8)

    for j in range(n_labels, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle("ROC Curves per Rule", fontsize=13, y=1.01)
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "roc_curves.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved: roc_curves.png")


def plot_feature_importance(best_model_name, best_model, feat_cols):
    """Horizontal bar chart of averaged feature importance."""
    try:
        imps   = get_importances(best_model, feat_cols)
        order  = np.argsort(imps)
        labels = [feat_cols[i] for i in order]
        vals   = imps[order]

        fig, ax = plt.subplots(figsize=(8, max(4, len(feat_cols) * 0.35)))
        bars = ax.barh(labels, vals, color=COLORS.get(best_model_name, "#888"),
                       alpha=0.85, edgecolor="white")
        ax.set_xlabel("Mean Importance (avg across rules)")
        ax.set_title(f"Feature Importance — {best_model_name}")
        for bar, val in zip(bars, vals):
            ax.text(val + 0.001, bar.get_y() + bar.get_height()/2,
                    f"{val:.3f}", va="center", fontsize=8)
        plt.tight_layout()
        plt.savefig(PLOTS_DIR / "feature_importance.png", dpi=150)
        plt.close()
        print("Saved: feature_importance.png")
    except AttributeError:
        print("Feature importance not available for this model type.")


def plot_confusion_matrices(best_model_name, y_test, y_pred):
    """2×4 grid of confusion matrices, one per rule."""
    n_labels = len(LABEL_COLS)
    cols = 4
    rows = 2
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 3, rows * 3))
    axes = axes.flatten()

    for i, col in enumerate(LABEL_COLS):
        ax = axes[i]
        cm = confusion_matrix(y_test[:, i], y_pred[:, i], labels=[0, 1])
        im = ax.imshow(cm, interpolation="nearest", cmap=plt.cm.Blues)
        ax.set_title(col, fontsize=11)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["OFF", "ON"])
        ax.set_yticks([0, 1]); ax.set_yticklabels(["OFF", "ON"])
        ax.set_xlabel("Predicted"); ax.set_ylabel("Actual")
        for r in range(2):
            for c in range(2):
                val = cm[r, c] if r < cm.shape[0] and c < cm.shape[1] else 0
                ax.text(c, r, str(val), ha="center", va="center",
                        color="white" if val > cm.max() / 2 else "black",
                        fontsize=12)

    plt.suptitle(f"Confusion Matrices — {best_model_name}", fontsize=13)
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "confusion_matrices.png", dpi=150)
    plt.close()
    print("Saved: confusion_matrices.png")


def plot_auc_heatmap(results):
    """Heatmap: AUC per model × rule."""
    models = list(results.keys())
    auc_matrix = np.array([results[m]["auc_per"] for m in models])

    fig, ax = plt.subplots(figsize=(10, 2.5 + 0.5 * len(models)))
    im = ax.imshow(auc_matrix, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
    plt.colorbar(im, ax=ax, label="AUC")
    ax.set_xticks(range(len(LABEL_COLS))); ax.set_xticklabels(LABEL_COLS, fontsize=11)
    ax.set_yticks(range(len(models)));   ax.set_yticklabels(models, fontsize=11)
    ax.set_title("ROC-AUC per Model × Rule")
    for i in range(len(models)):
        for j in range(len(LABEL_COLS)):
            val = auc_matrix[i, j]
            text = f"{val:.2f}" if not np.isnan(val) else "N/A"
            ax.text(j, i, text, ha="center", va="center", fontsize=9,
                    color="black" if 0.35 < val < 0.75 else "white")
    plt.tight_layout()
    plt.savefig(PLOTS_DIR / "auc_heatmap.png", dpi=150)
    plt.close()
    print("Saved: auc_heatmap.png")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print("Usage: python train.py <dataset_csv> [test_size]")
        sys.exit(1)

    csv_path  = Path(sys.argv[1])
    test_size = float(sys.argv[2]) if len(sys.argv) >= 3 else 0.2

    # ── Load ─────────────────────────────────────────────────────────────────
    X, y, feat_cols = load_data(csv_path)

    # Stratify on sum of labels (0-8) as a proxy for class balance
    strat_key = y.sum(axis=1)
    # Merge rare strata so stratify doesn't fail
    counts = pd.Series(strat_key).value_counts()
    rare   = counts[counts < 2].index
    strat_key_safe = np.where(np.isin(strat_key, rare), -1, strat_key)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=test_size,
        random_state=42,
        stratify=strat_key_safe,
    )
    print(f"Train: {len(X_train)}  Test: {len(X_test)}\n")

    # ── Build candidates ──────────────────────────────────────────────────────
    # Each entry: (display_name, base_clf, param_grid)
    model_specs = [
        ("Random Forest",
         RandomForestClassifier(random_state=42),
         RF_PARAM_GRID),
        ("Extra Trees",
         ExtraTreesClassifier(random_state=42),
         ET_PARAM_GRID),
        ("Gradient Boost",
         GradientBoostingClassifier(random_state=42),
         GB_PARAM_GRID),
        ("Decision Tree",
         DecisionTreeClassifier(random_state=42),
         DT_PARAM_GRID),
        ("Logistic Reg",
         LogisticRegression(random_state=42),
         LR_PARAM_GRID),
        ("SVM",
         SVC(probability=True, random_state=42),
         SVM_PARAM_GRID),
        ("KNN",
         KNeighborsClassifier(),
         KNN_PARAM_GRID),
        ("Naive Bayes",
         GaussianNB(),
         NB_PARAM_GRID),
    ]
    if HAS_XGB:
        model_specs.append((
            "XGBoost",
            XGBClassifier(random_state=42, eval_metric="logloss",
                          verbosity=0, use_label_encoder=False),
            XGB_PARAM_GRID,
        ))

    candidates = {}
    for name, clf, grid in model_specs:
        pipe = build_pipeline(clf)
        prefixed = prefix_grid(grid, pipeline_step="clf")
        print(f"=== Grid search: {name} ===")
        try:
            gs = run_grid_search(pipe, prefixed, X_train, y_train)
            candidates[name] = gs.best_estimator_
            print(f"Best params: {gs.best_params_}\n")
        except Exception as exc:
            print(f"FAILED ({exc}) — skipping {name}\n")

    # ── Evaluate ──────────────────────────────────────────────────────────────
    print("=== Test-set evaluation ===\n")
    results = {}
    for name, model in candidates.items():
        results[name] = evaluate(model, X_test, y_test, name)

    # Pick best by Hamming accuracy
    best_name  = max(results, key=lambda m: results[m]["hamming_acc"])
    best_model = candidates[best_name]
    best_res   = results[best_name]

    print(f"Best model: {best_name}  "
          f"(Hamming acc={best_res['hamming_acc']:.4f}, "
          f"Exact acc={best_res['exact_acc']:.4f})\n")

    # ── Save model ────────────────────────────────────────────────────────────
    model_path = WEIGHTS_DIR / "best_model.pkl"
    with open(model_path, "wb") as f:
        pickle.dump({"name": best_name, "model": best_model,
                     "feat_cols": feat_cols, "label_cols": LABEL_COLS}, f)
    print(f"Saved model → {model_path}")

    info_path = WEIGHTS_DIR / "best_model_info.txt"
    with open(info_path, "w") as f:
        f.write(f"Best model    : {best_name}\n")
        f.write(f"Hamming acc   : {best_res['hamming_acc']:.4f}\n")
        f.write(f"Exact-match   : {best_res['exact_acc']:.4f}\n")
        f.write(f"Features ({len(feat_cols)}): {feat_cols}\n")
        f.write(f"Labels ({len(LABEL_COLS)}):  {LABEL_COLS}\n\n")
        for col, f1, auc in zip(LABEL_COLS, best_res["f1_per"], best_res["auc_per"]):
            f.write(f"  {col:<10}  F1={f1:.3f}  AUC={auc:.3f}\n")
    print(f"Saved info  → {info_path}\n")

    # ── Plots ─────────────────────────────────────────────────────────────────
    print("=== Generating plots ===")
    plot_model_comparison(results, test_size)
    plot_per_rule_metrics(results)
    plot_roc_curves(results, y_test)
    plot_feature_importance(best_name, best_model, feat_cols)
    plot_confusion_matrices(best_name, y_test, best_res["y_pred"])
    plot_auc_heatmap(results)

    print(f"\nAll plots saved → {PLOTS_DIR}")
    print("Done.")


if __name__ == "__main__":
    main()
