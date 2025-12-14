
import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

from Cassiopeia.ILP.tree_utils import edge_cost, state_to_list, list_to_state


TOKEN_RE = re.compile(r"\(|\)|,|:|;|[^(),:;]+")

# newick parser
def parse_newick(s: str):
    tokens = TOKEN_RE.findall(s.strip())
    stack = []
    next_id = 0

    parent = {}
    children = {}
    label = {}

    def new_node():
        nonlocal next_id
        nid = next_id
        next_id += 1
        children[nid] = []
        return nid

    cur = new_node()
    root = cur
    stack.append(cur)

    i = 0
    while i < len(tokens):
        t = tokens[i]
        if t == "(":
            # start new internal node:
            internal = new_node()
            p = stack[-1]
            children[p].append(internal)
            parent[internal] = p
            stack.append(internal)
            i += 1
        elif t == ",":
            i += 1
        elif t == ")":
            stack.pop()
            i += 1
        elif t == ":":
            i += 2
        elif t == ";":
            i += 1
        else:
            leaf = new_node()
            p = stack[-1]
            children[p].append(leaf)
            parent[leaf] = p
            label[leaf] = t
            i += 1

    return parent, children, root, label


def leaf_label_to_state(label: str) -> Tuple[int, ...]:
    first = label.split("|")[0]
    state_str = first.split("_")[-1]
    return tuple(state_to_list(state_str))


# build state for each node:
def compute_parsimony(newick_path: Path, missing: int = -1) -> int:
    s = newick_path.read_text(encoding="utf-8").strip()
    parent, children, root, labels = parse_newick(s)


    node_state: Dict[int, Tuple[int, ...]] = {}

    order = []
    stack = [root]
    while stack:
        u = stack.pop()
        order.append(u)
        for v in children.get(u, []):
            stack.append(v)
    order.reverse()

    leaf_ids = [nid for nid in labels.keys()]
    if not leaf_ids:
        raise ValueError(f"No leaves found in {newick_path}")
    L = len(leaf_label_to_state(labels[leaf_ids[0]]))

    def consensus(child_states: List[Tuple[int, ...]]) -> Tuple[int, ...]:
        out = []
        for j in range(L):
            counts = {}
            for st in child_states:
                a = st[j]
                if a == missing:
                    continue
                counts[a] = counts.get(a, 0) + 1
            if counts:
                out.append(max(counts, key=counts.get))
            else:
                out.append(missing)
        return tuple(out)

    for u in order:
        if u in labels:
            node_state[u] = leaf_label_to_state(labels[u])
        else:
            child_states = [node_state[v] for v in children.get(u, []) if v in node_state]
            if child_states:
                node_state[u] = consensus(child_states)
            else:
                node_state[u] = tuple([0] * L)

    total = 0
    for v, p in parent.items():
        if p is None:
            continue
        total += edge_cost(node_state[p], node_state[v], missing=missing)
    return total

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("newick", type=str, nargs="+", help="Path(s) to Newick .nwk files")
    ap.add_argument("--missing", type=int, default=-1)
    args = ap.parse_args()

    for p in args.newick:
        path = Path(p)
        cost = compute_parsimony(path, missing=args.missing)
        print(f"{path.name}\tparsimony_cost={cost}")

if __name__ == "__main__":
    main()
