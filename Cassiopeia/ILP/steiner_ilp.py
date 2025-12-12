from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .potential_graph import PotentialGraph
from .tree_utils import State


@dataclass
class SteinerILPResult:
    selected_edges: Set[Tuple[State, State]]
    objective_value: float
    status: str

# check that each terminal is reachable
def _check_terminals_reachable(graph: PotentialGraph) -> None:
    root = graph.root
    terminals = graph.terminals

    stack = [root]
    seen = {root}
    while stack:
        u = stack.pop()
        for v in graph.out_edges.get(u, []):
            if v not in seen:
                seen.add(v)
                stack.append(v)

    unreachable = [t for t in terminals if t not in seen]
    if unreachable:
        raise ValueError(
            f"Infeasible: {len(unreachable)}/{len(terminals)} terminals unreachable from root in potential graph. "
            f"Example unreachable terminal: {unreachable[0]}"
        )

# main function to solve the steiner tree
#   takes a given PotentialGraph and returns a SteinerILPTree with selected_edges and objective
def solve_steiner_tree_ilp(
    graph: PotentialGraph,
    terminals: Optional[Sequence[State]] = None,
    time_limit_seconds: Optional[int] = None,
    msg: bool = False,
    mip_gap: Optional[float] = None,
) -> SteinerILPResult:
    try:
        import pulp
    except Exception as e:
        raise ImportError(
            "PuLP is required for steiner_ilp.py. Install with: pip install pulp"
        ) from e

    if terminals is None:
        terminals_list = list(graph.terminals)
    else:
        terminals_list = list(terminals)

    if not terminals_list:
        return SteinerILPResult(selected_edges=set(), objective_value=0.0, status="NO_TERMINALS")

    _check_terminals_reachable(graph)

    edges: List[Tuple[State, State]] = list(graph.costs.keys())
    costs: Dict[Tuple[State, State], int] = graph.costs

    root = graph.root
    nodes: List[State] = list(graph.nodes)
    node_set = set(nodes)

    # in/out adjacency
    in_edges: Dict[State, List[Tuple[State, State]]] = {v: [] for v in nodes}
    out_edges: Dict[State, List[Tuple[State, State]]] = {u: [] for u in nodes}
    for (u, v) in edges:
        if u not in node_set or v not in node_set:
            continue
        out_edges[u].append((u, v))
        in_edges[v].append((u, v))

    prob = pulp.LpProblem("Cassiopeia_Steiner_Tree", pulp.LpMinimize)

    # edge selection variables:
    y: Dict[Tuple[State, State], pulp.LpVariable] = {
        (u, v): pulp.LpVariable(f"y_{hash(u)}_{hash(v)}", lowBound=0, upBound=1, cat=pulp.LpBinary)
        for (u, v) in edges
    }
    # flow variables:
    f: Dict[Tuple[State, State, int], pulp.LpVariable] = {}
    for ti, t in enumerate(terminals_list):
        for (u, v) in edges:
            f[(u, v, ti)] = pulp.LpVariable(
                f"f_{hash(u)}_{hash(v)}_{ti}",
                lowBound=0,
                upBound=1,
                cat=pulp.LpContinuous,
            )

    # objective::
    prob += pulp.lpSum(costs[(u, v)] * y[(u, v)] for (u, v) in edges)

    for ti, _t in enumerate(terminals_list):
        for (u, v) in edges:
            prob += f[(u, v, ti)] <= y[(u, v)], f"cap_{hash(u)}_{hash(v)}_{ti}"

 
    # flow constraints for each terminal
    #  for each terminal, enforce these constraints:
        #   root: sum_out f - sum_in f = 1
        #   terminal: sum_in f - sum_out f = 1
        #   everywhere else: sum_out f - sum_in f = 0
    for ti, t in enumerate(terminals_list):
        if t not in node_set:
            raise ValueError("Provided terminal not present in graph.nodes.")

        for v in nodes:
            out_sum = pulp.lpSum(f[(u, w, ti)] for (u, w) in out_edges.get(v, []))
            in_sum = pulp.lpSum(f[(u, w, ti)] for (u, w) in in_edges.get(v, []))

            if v == root:
                prob += out_sum - in_sum == 1, f"flow_root_{ti}"
            elif v == t:
                prob += in_sum - out_sum == 1, f"flow_term_{ti}"
            else:
                prob += out_sum - in_sum == 0, f"flow_cons_{hash(v)}_{ti}"


    for t in terminals_list:
        if t == root:
            continue
        prob += pulp.lpSum(y[(u, v)] for (u, v) in in_edges.get(t, [])) >= 1, f"inc_{hash(t)}"

    solver = pulp.PULP_CBC_CMD(msg=msg)

    if time_limit_seconds is not None:
        try:
            solver.timeLimit = int(time_limit_seconds)
        except Exception:
            pass
    if mip_gap is not None:
        try:
            solver.gapRel = float(mip_gap)
        except Exception:
            pass

    status_code = prob.solve(solver)
    status_str = pulp.LpStatus.get(status_code, str(status_code))

    if status_str not in ("Optimal", "Not Solved", "Infeasible", "Unbounded", "Undefined"):
        pass

    if status_str == "Infeasible":
        return SteinerILPResult(selected_edges=set(), objective_value=float("inf"), status="Infeasible")

    # extract the selected edges:
    selected: Set[Tuple[State, State]] = set()
    for (u, v) in edges:
        val = y[(u, v)].value()
        if val is None:
            continue
        if val >= 0.5:
            selected.add((u, v))

    obj_val = float(pulp.value(prob.objective)) if prob.objective is not None else float("nan")
    return SteinerILPResult(selected_edges=selected, objective_value=obj_val, status=status_str)


# converts the selected edges into children map
def selected_edges_to_arborescence(
    root: State,
    selected_edges: Set[Tuple[State, State]],
) -> Dict[State, List[State]]:
    from collections import defaultdict, deque

    out_adj = defaultdict(list)
    for u, v in selected_edges:
        out_adj[u].append(v)

    children = defaultdict(list)
    parent: Dict[State, State] = {}

    q = deque([root])
    seen = {root}
    while q:
        u = q.popleft()
        for v in out_adj.get(u, []):
            if v not in parent and v != root:
                parent[v] = u
                children[u].append(v)
            if v not in seen:
                seen.add(v)
                q.append(v)

    for n in list(seen):
        children.setdefault(n, []) # every node should have a key

    return dict(children)
