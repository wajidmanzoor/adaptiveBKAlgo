#include <bits/stdc++.h>
using namespace std;

struct Graph {
    int n;
    vector<vector<int>> adj;
    bool connected(int u, int v) const {
        return binary_search(adj[u].begin(), adj[u].end(), v);
    }
};

vector<bool> skip_mask;   // marks vertices already covered
vector<int> order;        // global vertex order

// Check if R forms a clique
bool is_clique(const Graph& G, const vector<int>& R) {
    for (size_t i = 0; i < R.size(); i++)
        for (size_t j = i + 1; j < R.size(); j++)
            if (!G.connected(R[i], R[j])) return false;
    return true;
}

// Report and reorder after finding a maximal clique
void handle_clique(const vector<int>& R) {
    cout << "Found maximal clique: { ";
    for (int v : R) cout << v << " ";
    cout << "}\n";

    // Mark covered vertices
    for (int v : R) skip_mask[v] = true;

    // Reorder: uncovered first, covered last
    vector<int> uncovered, covered;
    for (int v : order) {
        if (skip_mask[v]) covered.push_back(v);
        else uncovered.push_back(v);
    }
    order.clear();
    order.insert(order.end(), uncovered.begin(), uncovered.end());
    order.insert(order.end(), covered.begin(), covered.end());
}

// Simple recursive exploration
void enumerate(const Graph& G, vector<int>& R, int start) {
    bool expanded = false;

    for (int idx = start; idx < (int)order.size(); idx++) {
        int v = order[idx];
        if (skip_mask[v]) continue; // skip redundant vertex

        vector<int> Rnext = R;
        Rnext.push_back(v);

        if (is_clique(G, Rnext)) {
            expanded = true;
            enumerate(G, Rnext, idx + 1);
        }
    }

    // Leaf node → maximal clique
    if (!expanded && !R.empty()) handle_clique(R);
}