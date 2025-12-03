import pandas as pd
from collections import defaultdict
import argparse

# --- Utilities ---

def state_to_list(state_str):
    return [int(c) for c in state_str]

def compute_parent_consensus(samples, states, missing=-1):
    """Compute consensus of a set of samples."""
    if not samples:
        return []
    num_sites = len(next(iter(states.values())))
    consensus = []
    for i in range(num_sites):
        counts = defaultdict(int)
        for s in samples:
            val = states[s][i]
            if val != missing:
                counts[val] += 1
        if counts:
            consensus.append(max(counts, key=counts.get))
        else:
            consensus.append(missing)
    return consensus

def best_split(samples, states, missing=-1):
    """Choose the site/state that best splits the samples into two non-empty groups."""
    num_sites = len(next(iter(states.values())))
    best_site, best_state = None, None
    best_score = -1

    for i in range(num_sites):
        counts = defaultdict(list)
        for s in samples:
            val = states[s][i]
            counts[val].append(s)
        for val, group in counts.items():
            if val == missing:
                continue
            left = group
            right = [s for s in samples if s not in left]
            if len(left) > 0 and len(right) > 0:
                score = min(len(left), len(right))
                if score > best_score:
                    best_score = score
                    best_site = i
                    best_state = val
    return best_site, best_state

def split_samples(samples, states, site, state, missing=-1):
    """Split samples into left/right based on mutation at site."""
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

def hamming_distance(state1, state2, missing=-1):
    """Compute Hamming distance between two states, ignoring missing values."""
    dist = 0
    for s1, s2 in zip(state1, state2):
        if s1 != missing and s2 != missing and s1 != s2:
            dist += 1
    return dist

# --- Recursive Newick Construction ---

def greedy_newick(samples, states, missing=-1, parent_consensus=None):
    if len(samples) == 0:
        return "", []

    if parent_consensus is None:
        parent_consensus = compute_parent_consensus(samples, states, missing)

    # Leaf node
    if len(samples) == 1:
        s = samples[0]
        state_str = "".join(str(c) for c in states[s])
        branch_len = hamming_distance(states[s], parent_consensus, missing)
        return f"{s}_{state_str}", states[s]

    # Internal node
    site, state = best_split(samples, states, missing)
    if site is None:
        # Cannot split further, return all leaves
        leaf_strings = []
        for s in samples:
            state_str = "".join(str(c) for c in states[s])
            branch_len = hamming_distance(states[s], parent_consensus, missing)
            leaf_strings.append(f"{s}_{state_str}:{branch_len}")
        return "(" + ",".join(leaf_strings) + ")", parent_consensus

    left, right = split_samples(samples, states, site, state, missing)
    if len(left) == 0 or len(right) == 0:
        leaf_strings = []
        for s in samples:
            state_str = "".join(str(c) for c in states[s])
            branch_len = hamming_distance(states[s], parent_consensus, missing)
            leaf_strings.append(f"{s}_{state_str}:{branch_len}")
        return "(" + ",".join(leaf_strings) + ")", parent_consensus

    # Compute child consensus states
    left_consensus = compute_parent_consensus(left, states, missing)
    right_consensus = compute_parent_consensus(right, states, missing)

    # Recurse with child consensus as new parent
    left_newick, left_state = greedy_newick(left, states, missing, left_consensus)
    right_newick, right_state = greedy_newick(right, states, missing, right_consensus)

    # Branch lengths from parent to children
    left_branch_len = hamming_distance(left_consensus, parent_consensus, missing)
    right_branch_len = hamming_distance(right_consensus, parent_consensus, missing)

    newick = f"({left_newick}:{left_branch_len},{right_newick}:{right_branch_len})"
    return newick, parent_consensus

# --- Main ---

def main():
    parser = argparse.ArgumentParser(description="Greedy perfect phylogeny reconstruction")
    parser.add_argument("txt_file", help="Tab-delimited TXT with columns: cell, state")
    parser.add_argument("-o", "--output", default="tree.newick", help="Output Newick file")
    parser.add_argument("-m", "--missing", type=int, default=-1, help="Missing state indicator")
    args = parser.parse_args()

    df = pd.read_csv(args.txt_file, sep="\t", index_col=0, dtype=str)
    df.fillna(str(args.missing), inplace=True)

    states = {s: state_to_list(df.loc[s, "state"]) for s in df.index}

    newick, _ = greedy_newick(list(df.index), states, args.missing)
    newick += ";"

    print(newick)
    with open(args.output, "w") as f:
        f.write(newick)
    print(f"Saved to {args.output}")

if __name__ == "__main__":
    main()