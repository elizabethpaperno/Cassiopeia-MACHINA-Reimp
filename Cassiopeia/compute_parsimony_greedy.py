#!/usr/bin/env python3
"""
compute_parsimony_greedy.py

Compute parsimony cost for *Greedy* trees by summing the branch lengths
encoded in the Newick string.

Your Greedy Newick uses a pattern like:
  leaf:5:0
or
  ):1
or
  ):1:0

We interpret the FIRST number after ':' for each edge as the edge cost,
and ignore any additional ':<num>' that appears before the next delimiter
(',', ')', or ';').

Usage:
  python compute_parsimony_greedy.py path/to/greedy_tree.nwk [...]
"""

import argparse
from pathlib import Path


DELIMS = {",", ")", ";"}


def sum_first_branch_lengths(newick: str) -> float:
    """
    Scan the Newick string and sum the first branch length per edge.

    Implementation detail:
    After a label or ')', Newick may have ':<num>' (branch length).
    Your files sometimes have ':<num>:<num>' — we count only the first.
    We reset the "already counted a length for this edge" flag at delimiters.
    """
    s = newick.strip()
    total = 0.0

    i = 0
    counted_for_edge = False  # have we already counted the first :len for the current edge?

    while i < len(s):
        ch = s[i]

        if ch == ":":
            # parse number after ':'
            j = i + 1
            # number can be int or float
            while j < len(s) and (s[j].isdigit() or s[j] == "."):
                j += 1

            num_str = s[i + 1 : j]
            if num_str and not counted_for_edge:
                try:
                    total += float(num_str)
                    counted_for_edge = True
                except ValueError:
                    pass

            i = j
            continue

        if ch in DELIMS:
            # next edge starts after delimiter
            counted_for_edge = False
            i += 1
            continue

        # normal character
        i += 1

    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("newick", nargs="+", help="Greedy .nwk files")
    args = ap.parse_args()

    for p in args.newick:
        path = Path(p)
        txt = path.read_text(encoding="utf-8").strip()
        cost = sum_first_branch_lengths(txt)
        # cost should be integer-valued in your data, but keep float-safe printing
        if abs(cost - round(cost)) < 1e-9:
            cost_str = str(int(round(cost)))
        else:
            cost_str = f"{cost:.6f}"
        print(f"{path.name}\tparsimony_cost={cost_str}")


if __name__ == "__main__":
    main()

