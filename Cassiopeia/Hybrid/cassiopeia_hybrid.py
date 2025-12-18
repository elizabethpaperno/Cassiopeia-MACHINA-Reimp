import argparse
from collections import defaultdict

from Cassiopeia.ILP.tree_utils import (
    parse_character_matrix,
    state_to_str,
    hamming_distance,
    compute_consensus,
    edge_cost,
    to_newick,
    edges_to_rooted_tree
)
from Cassiopeia.ILP.potential_graph import build_potential_graph_from_cells
from Cassiopeia.ILP.steiner_ilp import solve_steiner_tree_ilp

# adapted methods from greedy implementation

# compute consensus for set of samples
def compute_parent_consensus(samples, states, missing=-1):
    if not samples:
        return tuple()
    
    state_seqs = [states[s] for s in samples]
    return compute_consensus(state_seqs, missing=missing)

# split into non-empty groups
def best_split_hybrid(samples, states, missing=-1):
    num_sites = len(next(iter(states.values())))
    best_site, best_val = None, None
    best_score = -1

    for i in range(num_sites):
        counts = defaultdict(list)

        for s in samples:
            v = states[s][i]
            if v != missing:
                counts[v].append(s)

        # take at least 2 distinct values
        if len(counts) < 2:
            continue

        for val, left in counts.items():
            right_size = len(samples) - len(left)
            if right_size == 0:
                continue

            score = min(len(left), right_size)
            if score > best_score:
                best_score = score
                best_site = i
                best_val = val

    return best_site, best_val

# split samples into left/right based on mutation at site
def split_samples_hybrid(samples, states, site, state, missing=-1):
    left, right, missing_cells = [], [], []
    for s in samples:
        val = states[s][site]
        if val == missing:
            missing_cells.append(s)
        elif val == state:
            left.append(s)
        else:
            right.append(s)
    if missing_cells:
        if len(left) >= len(right):
            left.extend(missing_cells)
        else:
            right.extend(missing_cells)
    
    return left, right


# recursive hybrid newick construction
def greedy_newick_hybrid(samples, states, ilp_threshold=10, missing=-1, unmutated=0, max_imputations=256, include_all_pair_edges=False, parent_consensus=None, depth=0, debug=False,):
    if len(samples) == 0:
        return "", tuple()

    consensus = compute_parent_consensus(samples, states, missing)
    if parent_consensus is None:
        parent_consensus = consensus

    # base case
    if len(samples) == 1:
        s = samples[0]
        if debug:
            print(f"Leaf: {s}")
        return f"{s}_{state_to_str(states[s])}", states[s]

    # ILP
    if ilp_threshold > 0 and len(samples) <= ilp_threshold:
        if debug:
            print(f"ILP: {len(samples)} cells")

        subset_states = {s: states[s] for s in samples}

        try:
            # build potential graph
            G, _ = build_potential_graph_from_cells(
                states_by_cell=subset_states,
                missing=missing,
                unmutated=unmutated,
                max_imputations_per_terminal=max_imputations,
                include_all_pair_edges=include_all_pair_edges,
            )

            res = solve_steiner_tree_ilp(G, msg=False)

            if res.status == "Infeasible":
                raise RuntimeError("ILP infeasible")

            rooted = edges_to_rooted_tree(G.root, res.selected_edges)

            # create new leaf nodes:
            CELL = -777777
            labels = {}
            for cell, st in subset_states.items():
                leaf = (CELL, int(cell))
                labels[leaf] = f"{cell}_{state_to_str(st)}"
                rooted.children.setdefault(st, []).append(leaf)
                rooted.children[leaf] = []
                rooted.parent[leaf] = st

            ilp_newick = to_newick(
                rooted,
                leaf_label_fn=lambda n: labels.get(n, ""),
                internal_label_fn=lambda _: "",
                branch_length_fn=lambda u, v: 0 if isinstance(v, tuple) else edge_cost(u, v, missing),
            ).rstrip(";")

            if debug:
                print(f"ILP solved")

            return ilp_newick, G.root

        except Exception as e:
            leaves = [
                f"{s}_{state_to_str(states[s])}:{hamming_distance(states[s], parent_consensus, missing)}"
                for s in samples
            ]
            return f"({','.join(leaves)})", consensus

    # greedy
    site, state_val = best_split_hybrid(samples, states, missing)
    if site is None:
        # cannot split further, return all leaves
        leaves = [
            f"{s}_{state_to_str(states[s])}:{hamming_distance(states[s], parent_consensus, missing)}"
            for s in samples
        ]
        return f"({','.join(leaves)})", consensus

    left, right = split_samples_hybrid(samples, states, site, state_val, missing)

    if not left or not right:
        leaves = [
            f"{s}_{state_to_str(states[s])}:{hamming_distance(states[s], parent_consensus, missing)}"
            for s in samples
        ]
        return f"({','.join(leaves)})", consensus

    if debug:
        print(f"split site {site}={state_val}: {len(left)}|{len(right)}")

    # compute child consensus and recurse with child as new parent
    left_newick, left_state = greedy_newick_hybrid(
        left, states, ilp_threshold, missing, unmutated,
        max_imputations, include_all_pair_edges,
        compute_parent_consensus(left, states, missing),
        depth + 1, debug
    )
    right_newick, right_state = greedy_newick_hybrid(
        right, states, ilp_threshold, missing, unmutated,
        max_imputations, include_all_pair_edges,
        compute_parent_consensus(right, states, missing),
        depth + 1, debug
    )

    left_branch_len = hamming_distance(left_state, consensus, missing)
    right_branch_len = hamming_distance(right_state, consensus, missing)

    return f"({left_newick}:{left_branch_len},{right_newick}:{right_branch_len})", consensus

def run_hybrid(txt_file, output_newick, ilp_threshold=0, missing=-1, unmutated=0, max_imputations=256, include_all_pair_edges=False, debug=False):
    
    states_by_cell = parse_character_matrix(txt_file, missing=missing)
    print(f"Loaded {len(states_by_cell)} cells")
    
    # build hybrid tree
    samples = list(states_by_cell.keys())
    newick, _ = greedy_newick_hybrid(
        samples=samples,
        states=states_by_cell,
        ilp_threshold=ilp_threshold,
        missing=missing,
        unmutated=unmutated,
        max_imputations=max_imputations,
        include_all_pair_edges=include_all_pair_edges,
        parent_consensus=None,
        depth=0,
        debug=debug,
    )
    
    newick += ";"

    with open(output_newick, "w") as f:
        f.write(newick)
    
    print(f"\noutput: {output_newick}")

# adapted from cassiopeia_ilp.py
def main():
    parser = argparse.ArgumentParser(
        description="Cassiopeia Hybrid: Greedy splitting + ILP optimization"
    )
    parser.add_argument(
        "txt_file",
        help="Input tab-delimited file with index=cell and column 'state'"
    )
    parser.add_argument(
        "-o", "--output",
        default="tree_hybrid.newick",
        help="Output Newick filename"
    )
    parser.add_argument(
        "-t", "--threshold",
        type=int,
        default=0,
        help="ILP threshold: group size to switch to ILP (pure greedy=0)"
    )
    parser.add_argument(
        "-m", "--missing",
        type=int,
        default=-1,
        help="Missing value sentinel (default: -1)"
    )
    parser.add_argument(
        "--unmutated",
        type=int,
        default=0,
        help="Founder/unmutated allele (default: 0)"
    )
    parser.add_argument(
        "--max-imputations",
        type=int,
        default=256,
        help="Max imputations per terminal/cell state with missing entries"
    )
    parser.add_argument(
        "--dense",
        action="store_true",
        help="Use dense potential graph (slower but may reduce infeasibility)"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print debug summary"
    )
    
    args = parser.parse_args()
    
    run_hybrid(
        txt_file=args.txt_file,
        output_newick=args.output,
        ilp_threshold=args.threshold,
        missing=args.missing,
        unmutated=args.unmutated,
        max_imputations=args.max_imputations,
        include_all_pair_edges=args.dense,
        debug=args.debug,
    )


if __name__ == "__main__":
    main()