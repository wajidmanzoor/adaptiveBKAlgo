import os
import sys


def convert_file(input_path, output_path):
    edges = set()
    with open(input_path, 'r') as f:
        lines = f.read().splitlines()

    # skip empty lines
    lines = [l.strip() for l in lines if l.strip()]

    # first line: num_vertices num_edges (informational, not used)
    idx = 1

    for line in lines[idx:]:
        parts = list(map(int, line.split()))
        if len(parts) < 1:
            continue
        u = parts[0]
        for v in parts[1:]:
            if u != v:
                edge = (min(u, v), max(u, v))
                edges.add(edge)

    with open(output_path, 'w') as f:
        for u, v in sorted(edges):
            f.write(f"{u} {v}\n")

    return len(edges)


def main():
    if len(sys.argv) < 2:
        print("Usage: python convert_to_edgelist.py <input_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    if not os.path.isdir(input_dir):
        print(f"Error: '{input_dir}' is not a directory.")
        sys.exit(1)

    output_dir = os.path.join(os.getcwd(), "edgeListCleaned")
    os.makedirs(output_dir, exist_ok=True)

    files = sorted(f for f in os.listdir(input_dir)
                   if os.path.isfile(os.path.join(input_dir, f)))

    if not files:
        print(f"No files found in '{input_dir}'.")
        sys.exit(1)

    print(f"Output dir : {output_dir}")
    print(f"Files found: {len(files)}")
    print("-" * 40)

    for fname in files:
        input_path = os.path.join(input_dir, fname)
        output_path = os.path.join(output_dir, fname + ".clean")
        try:
            num_edges = convert_file(input_path, output_path)
            print(f"  {fname} -> {fname}.clean  ({num_edges} edges)")
        except Exception as e:
            print(f"  ERROR on {fname}: {e}")

    print("-" * 40)
    print("Done.")


if __name__ == "__main__":
    main()
