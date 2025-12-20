import pandas as pd
from collections import defaultdict
import argparse

def state_to_list(state_str):
    return [int(c) for c in state_str]

def compute_parent_consensus(samples, states, missing=-1):
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
    finalsite, finalstate = None, None
    best_support = -1

    num_sites = len(next(iter(states.values())))

    for i in range(num_sites):
        counts = defaultdict(int)
        for s in samples:
            val = states[s][i]
            if val != missing:
                counts[val] += 1
        for val, count in counts.items():
            if val == missing:
                continue
            if 0 < count < len(samples):
                if count > best_support:
                    best_support = count
                    finalsite = i
                    finalstate = val

    return finalsite, finalstate


def split_samples(samples, states, site, state, missing=-1):
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
    dist = 0
    for s1, s2 in zip(state1, state2):
        if s1 != missing and s2 != missing and s1 != s2:
            dist += 1
    return dist


def greedy_newick(samples, states, missing=-1, parent_consensus=None):
    if len(samples) == 0:
        return "", []

    if parent_consensus is None:
        parent_consensus = compute_parent_consensus(samples, states, missing)

    if len(samples) == 1:
        s = samples[0]
        state_str = "".join(str(c) for c in states[s])
        branch_len = hamming_distance(states[s], parent_consensus, missing)
        return f"{s}_{state_str}", states[s]

    site, state = best_split(samples, states, missing)
    if site is None:
        leaf_strs = []
        for s in samples:
            state_str = "".join(str(c) for c in states[s])
            branch_len = hamming_distance(states[s], parent_consensus, missing)
            leaf_strs.append(f"{s}_{state_str}:{branch_len}")
        return "(" + ",".join(leaf_strs) + ")", parent_consensus

    left, right = split_samples(samples, states, site, state, missing)
    if len(left) == 0 or len(right) == 0:
        leaf_strs = []
        for s in samples:
            state_str = "".join(str(c) for c in states[s])
            branch_len = hamming_distance(states[s], parent_consensus, missing)
            leaf_strs.append(f"{s}_{state_str}:{branch_len}")
        return "(" + ",".join(leaf_strs) + ")", parent_consensus

    lconsensus = compute_parent_consensus(left, states, missing)
    rconsensus = compute_parent_consensus(right, states, missing)

    ltree, lstate = greedy_newick(left, states, missing, lconsensus)
    rtree, rstate = greedy_newick(right, states, missing, rconsensus)

    branchl = hamming_distance(lconsensus, parent_consensus, missing)
    branchr = hamming_distance(rconsensus, parent_consensus, missing)

    newick = f"({ltree}:{branchl},{rtree}:{branchr})"
    return newick, parent_consensus

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