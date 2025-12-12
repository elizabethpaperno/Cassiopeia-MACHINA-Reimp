from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from collections import defaultdict, deque
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

from .tree_utils import (
    State,
    list_to_state,
    edge_cost,
    is_irreversible_transition,
)


@dataclass
class PotentialGraph:
    root: State
    terminals: Set[State]
    nodes: Set[State]
    out_edges: Dict[State, List[State]]
    costs: Dict[Tuple[State, State], int]

    def edges(self) -> Iterable[Tuple[State, State, int]]:
        for (u, v), c in self.costs.items():
            yield u, v, c


# collect alleles in data for every site
def _collect_site_alleles(
    observed_states: Iterable[State],
    unmutated: int = 0,
    missing: int = -1,
) -> Dict[int, Set[int]]:
    
    obs = list(observed_states)
    if not obs:
        return {}

    L = len(obs[0])
    alleles: Dict[int, Set[int]] = {i: set([unmutated]) for i in range(L)}
    for st in obs:
        for i, a in enumerate(st):
            if a == missing:
                continue
            alleles[i].add(int(a))
    return alleles


# turn state with possible 'missing' entries into list of fully-imputed states
def _impute_missing_states(
    state: State,
    site_alleles: Dict[int, Set[int]],
    missing: int = -1,
    max_imputations: int = 256,
) -> List[State]:
    
    missing_positions = [i for i, a in enumerate(state) if a == missing]
    if not missing_positions:
        return [state]

    # for every missing site, build a choice list
    choices: List[List[int]] = []
    for i in missing_positions:
        opts = sorted(site_alleles.get(i, {0}))
        choices.append(opts)

    total = 1
    for ch in choices:
        total *= max(1, len(ch))
        if total > max_imputations:
            break

    results: List[State] = []
    # count products until reaching max_imputations
    for combo in product(*choices):
        lst = list(state)
        for idx, site in enumerate(missing_positions):
            lst[site] = int(combo[idx])
        results.append(list_to_state(lst))
        if len(results) >= max_imputations:
            break
    return results


# generates set of plausible ancesstors
#     repeatedly reverts mutated sites back to unmuted until reaching root
def _ancestor_chain_states(
    state: State,
    root: State,
    unmutated: int = 0,
) -> Set[State]:
    ancestors: Set[State] = set()
    ancestors.add(state)

    if state == root:
        return ancestors

    #  queue over states:
    q = deque([state])
    seen = {state}

    L = len(state)
    while q:
        cur = q.popleft()
        if cur == root:
            continue # if cur = root - don't need to keep expanding 
        for i in range(L):
            if cur[i] != unmutated:
                nxt = list(cur)
                nxt[i] = unmutated
                nxt_t = list_to_state(nxt)
                if nxt_t not in seen:
                    seen.add(nxt_t)
                    ancestors.add(nxt_t)
                    q.append(nxt_t)

    ancestors.add(root) # make sure that the root is included
    return ancestors


# main function to build potential graph, returns a PotentialGraph with nodes, directed edges, and edge costs.
def build_potential_graph(
    observed_terminals: Iterable[State],
    missing: int = -1,
    unmutated: int = 0,
    max_imputations_per_terminal: int = 256,
    include_all_pair_edges: bool = False,
) -> PotentialGraph:

    observed_terminals = list(observed_terminals)
    if not observed_terminals:
        raise ValueError("No observed terminal states provided.")

    L = len(observed_terminals[0])
    root = list_to_state([unmutated] * L)

    site_alleles = _collect_site_alleles(observed_terminals, unmutated=unmutated, missing=missing)

    terminals: Set[State] = set()
    for st in observed_terminals:
        for imp in _impute_missing_states(st, site_alleles, missing=missing, max_imputations=max_imputations_per_terminal):
            terminals.add(imp)

    # build nodes:
    nodes: Set[State] = set()
    nodes.add(root)
    for t in terminals:
        nodes.update(_ancestor_chain_states(t, root=root, unmutated=unmutated))

    # create adjacency and costs:
    out_edges: Dict[State, List[State]] = defaultdict(list)
    costs: Dict[Tuple[State, State], int] = {}

    node_list = list(nodes)

    if include_all_pair_edges: # dense --> connects irreversible parent-child pairs
        for u in node_list:
            for v in node_list:
                if u == v:
                    continue
                if is_irreversible_transition(u, v, missing=missing, unmutated=unmutated):
                    c = edge_cost(u, v, missing=missing)
                    if c <= 0: # try to avoid zero-cost edges unless identical
                        continue
                    out_edges[u].append(v)
                    costs[(u, v)] = c
    else: # sparse --> single-site mutation edges only 
        index: Dict[Tuple[int, ...], List[State]] = defaultdict(list)

        SENT = -999999 # for each node create signature with 1 position replaced by sentinel
        for st in node_list:
            for i in range(L):
                sig = list(st)
                sig[i] = SENT
                index[tuple(sig)].append(st)

        # connect pairs that have one position difference:
        for bucket in index.values():
            if len(bucket) < 2:
                continue
            # for each pair, check irreversibility and add edge
            for u in bucket:
                for v in bucket:
                    if u == v:
                        continue
                    diff = 0
                    for a, b in zip(u, v):
                        if a != b:
                            diff += 1
                            if diff > 1:
                                break
                    if diff != 1:
                        continue

                    if not is_irreversible_transition(u, v, missing=missing, unmutated=unmutated):
                        continue
                    c = edge_cost(u, v, missing=missing)
                    if c <= 0:
                        continue
                    out_edges[u].append(v)
                    costs[(u, v)] = c

    for n in nodes:
        out_edges.setdefault(n, []) # check that every node should have a key in out_edges

    return PotentialGraph(
        root=root,
        terminals=terminals,
        nodes=nodes,
        out_edges=dict(out_edges),
        costs=costs,
    )


# wrapper for starting from parse_character_matrix output:
def build_potential_graph_from_cells(
    states_by_cell: Dict[str, State],
    missing: int = -1,
    unmutated: int = 0,
    max_imputations_per_terminal: int = 256,
    include_all_pair_edges: bool = False,
) -> Tuple[PotentialGraph, Dict[str, Set[State]]]:
  

    site_alleles = _collect_site_alleles(states_by_cell.values(), unmutated=unmutated, missing=missing)

    cell_to_terminals: Dict[str, Set[State]] = {}
    all_terminals: List[State] = []
    for cell, st in states_by_cell.items():
        imps = set(_impute_missing_states(st, site_alleles, missing=missing, max_imputations=max_imputations_per_terminal))
        cell_to_terminals[cell] = imps
        all_terminals.extend(list(imps))

    G = build_potential_graph(
        observed_terminals=all_terminals,
        missing=missing,
        unmutated=unmutated,
        max_imputations_per_terminal=max_imputations_per_terminal,
        include_all_pair_edges=include_all_pair_edges,
    )
    return G, cell_to_terminals
