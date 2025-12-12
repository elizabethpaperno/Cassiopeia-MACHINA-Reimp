# main runner for cassiopeia =-ILP
# parses character matrix --> builds potential graph from cells --> solves steiner tree --> convertes selected edges --> Newick tree

from __future__ import annotations

import argparse
import os
from typing import Optional

from .tree_utils import (
    parse_character_matrix,
    edges_to_rooted_tree,
    to_newick,
    build_leaf_label_fn_from_cells,
    edge_cost,
)
from .potential_graph import build_potential_graph_from_cells
from .steiner_ilp import solve_steiner_tree_ilp


def _write_text(path: str, s: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(s)


def run_ilp(
    txt_file: str,
    output_newick: str,
    missing: int = -1,
    unmutated: int = 0,
    max_imputations_per_terminal: int = 256,
    include_all_pair_edges: bool = False,
    time_limit_seconds: Optional[int] = None,
    mip_gap: Optional[float] = None,
    solver_msg: bool = False,
    write_debug: bool = False,
) -> None:
    # parse input
    states_by_cell = parse_character_matrix(txt_file, missing=missing)

    # build potential graph:
    G, cell_to_terminal_states = build_potential_graph_from_cells(
        states_by_cell=states_by_cell,
        missing=missing,
        unmutated=unmutated,
        max_imputations_per_terminal=max_imputations_per_terminal,
        include_all_pair_edges=include_all_pair_edges,
    )

    # solve Steiner ILP
    res = solve_steiner_tree_ilp(
        graph=G,
        terminals=None,
        time_limit_seconds=time_limit_seconds,
        msg=solver_msg,
        mip_gap=mip_gap,
    )

    if res.status == "Infeasible":
        raise RuntimeError("ILP returned Infeasible. Try include_all_pair_edges=True or increase max_imputations.")

    # build rooted tree + Newick
    rooted = edges_to_rooted_tree(G.root, res.selected_edges)

    # leaf labels:
    leaf_label_fn = build_leaf_label_fn_from_cells(states_by_cell)

    # branch lengths:
    def branch_len(u, v):
        return edge_cost(u, v, missing=missing)

    newick = to_newick(
        rooted,
        leaf_label_fn=leaf_label_fn,
        internal_label_fn=lambda _st: "",
        branch_length_fn=branch_len,
    )

    _write_text(output_newick, newick)

    # write all the current states out for debugging;
    if write_debug:
        dbg_path = output_newick + ".debug.txt"
        lines = []
        lines.append(f"input_file: {txt_file}")
        lines.append(f"missing: {missing}, unmutated: {unmutated}")
        lines.append(f"max_imputations_per_terminal: {max_imputations_per_terminal}")
        lines.append(f"include_all_pair_edges: {include_all_pair_edges}")
        lines.append("")
        lines.append(f"num_cells: {len(states_by_cell)}")
        lines.append(f"num_terminals (after imputation): {len(G.terminals)}")
        lines.append(f"num_nodes (potential graph): {len(G.nodes)}")
        lines.append(f"num_edges (potential graph): {len(G.costs)}")
        lines.append("")
        lines.append(f"ILP status: {res.status}")
        lines.append(f"ILP objective: {res.objective_value}")
        lines.append(f"selected_edges: {len(res.selected_edges)}")
        lines.append("")
        max_terms = max(len(v) for v in cell_to_terminal_states.values()) if cell_to_terminal_states else 0
        lines.append(f"max imputed terminals for any cell: {max_terms}")
        _write_text(dbg_path, "\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Cassiopeia ILP (Steiner tree) reconstruction")
    parser.add_argument("txt_file", type=str, help="Input tab-delimited file with index=cell and column 'state'")
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="tree_ilp.newick",
        help="Output Newick filename"
    )
    parser.add_argument(
        "-m", "--missing",
        type=int,
        default=-1,
        help="Missing value sentinel (must match data preprocessing; Greedy uses -1 by default)"
    )
    parser.add_argument(
        "--unmutated",
        type=int,
        default=0,
        help="Founder/unmutated allele (default 0)"
    )
    parser.add_argument(
        "--max-imputations",
        type=int,
        default=256,
        help="Max number of imputations per terminal/cell state with missing entries"
    )
    parser.add_argument(
        "--dense",
        action="store_true",
        help="Use dense potential graph edges (any irreversible transition). Slower but can reduce infeasibility."
    )
    parser.add_argument(
        "--time-limit",
        type=int,
        default=None,
        help="ILP solver time limit in seconds"
    )
    parser.add_argument(
        "--mip-gap",
        type=float,
        default=None,
        help="Relative MIP gap"
    )
    parser.add_argument(
        "--solver-msg",
        action="store_true",
        help="Print ILP solver logs"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Write debug summary next to output"
    )

    args = parser.parse_args()

    run_ilp(
        txt_file=args.txt_file,
        output_newick=args.output,
        missing=args.missing,
        unmutated=args.unmutated,
        max_imputations_per_terminal=args.max_imputations,
        include_all_pair_edges=args.dense,
        time_limit_seconds=args.time_limit,
        mip_gap=args.mip_gap,
        solver_msg=args.solver_msg,
        write_debug=args.debug,
    )

    print(f"Wrote ILP tree to: {args.output}")


if __name__ == "__main__":
    main()
