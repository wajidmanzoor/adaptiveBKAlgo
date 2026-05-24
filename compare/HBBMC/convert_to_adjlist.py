import os
import sys
from collections import defaultdict


def convert_file(input_path, output_path):
    adj = defaultdict(set)

    with open(input_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            u, v = map(int, line.split())
            if u != v:
                adj[u].add(v)
                adj[v].add(u)

    if not adj:
        return 0, 0

    vertices = sorted(adj.keys())
    num_vertices = max(vertices) + 1
    num_edges = sum(len(nbrs) for nbrs in adj.values()) // 2

    with open(output_path, 'w') as f:
        f.write(f"{num_vertices} {num_edges}\n")
        for u in range(num_vertices):
            neighbors = sorted(adj.get(u, set()))
            f.write(" ".join(map(str, [u] + neighbors)) + "\n")

    return num_vertices, num_edges


def main():
    if len(sys.argv) < 2:
        print("Usage: python convert_to_adjlist.py <input_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    if not os.path.isdir(input_dir):
        print(f"Error: '{input_dir}' is not a directory.")
        sys.exit(1)

    output_dir = os.path.join(os.getcwd(), "adjListConverted")
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
        # strip .clean suffix if present, otherwise keep name as-is
        out_name = fname[:-6] if fname.endswith(".clean") else fname
        output_path = os.path.join(output_dir, out_name)
        try:
            v, e = convert_file(input_path, output_path)
            print(f"  {fname} -> {out_name}  (|V|={v}, |E|={e})")
        except Exception as e:
            print(f"  ERROR on {fname}: {e}")

    print("-" * 40)
    print("Done.")


if __name__ == "__main__":
    main()
