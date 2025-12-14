from __future__ import annotations

from dataclasses import dataclass
from collections import defaultdict, deque
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Set, Callable, Union

import pandas as pd


State = Tuple[int, ...]



# converts state to list (same as greedy)
def state_to_list(state_str: str) -> List[int]:
    return [int(c) for c in state_str]


def list_to_state(x: Sequence[int]) -> State:
    return tuple(int(v) for v in x)


def state_to_str(state: Sequence[int]) -> str:
    return "".join(str(int(c)) for c in state)


# parses char matrix - same parsing pattern as greedy, returns states by cell
def parse_character_matrix(
    txt_file: str,
    missing: int = -1,
) -> Dict[str, State]:
    
    df = pd.read_csv(txt_file, sep="\t", index_col=0, dtype=str)
    if "state" not in df.columns:
        raise ValueError(f"Expected a column named 'state' in {txt_file}, got columns={list(df.columns)}")

    df.fillna(str(missing), inplace=True)

    states_by_cell: Dict[str, State] = {}
    lengths: Set[int] = set()

    for cell in df.index:
        s = df.loc[cell, "state"]
        if not isinstance(s, str):
            s = str(s)
        lst = state_to_list(s)
        st = list_to_state(lst)
        states_by_cell[str(cell)] = st
        lengths.add(len(st))

    if len(lengths) != 1:
        raise ValueError(f"Inconsistent state lengths found: {sorted(lengths)} in {txt_file}")

    return states_by_cell


# comptue hamming distance (same as for greedy)
def hamming_distance(state1: Sequence[int], state2: Sequence[int], missing: int = -1) -> int:
    dist = 0
    for s1, s2 in zip(state1, state2):
        if s1 != missing and s2 != missing and s1 != s2:
            dist += 1
    return dist

# computes consensus - for each site counts non-missing alleles and picks most frequent allele
def compute_consensus(states: Iterable[Sequence[int]], missing: int = -1) -> State:
    states = list(states)
    if not states:
        return tuple()

    num_sites = len(states[0])
    consensus: List[int] = []
    for i in range(num_sites):
        counts = defaultdict(int)
        for st in states:
            val = int(st[i])
            if val != missing:
                counts[val] += 1
        if counts:
            consensus.append(max(counts, key=counts.get))
        else:
            consensus.append(missing)

    return list_to_state(consensus)


# returns edge cost - hamming distance between child and parent consensus
def edge_cost(parent: Sequence[int], child: Sequence[int], missing: int = -1) -> int:
    return hamming_distance(parent, child, missing=missing)


# helper func for potential graph:
def is_irreversible_transition(
    parent: Sequence[int],
    child: Sequence[int],
    missing: int = -1,
    unmutated: int = 0,
) -> bool:
    for p, c in zip(parent, child):
        if p == missing or c == missing:
            continue
        if p != unmutated and c == unmutated:
            return False
        if p != unmutated and c != unmutated and p != c:
            return False
    return True


# build rooted tree from directed edges:
@dataclass(frozen=True)
class RootedTree:
    root: State
    children: Dict[State, List[State]]  # adjacency (rooted)
    parent: Dict[State, Optional[State]]  # parent pointers (root is None)


# converts directed edges --> rooted tree
def edges_to_rooted_tree(
    root: State,
    edges: Iterable[Tuple[State, State]],
) -> RootedTree:
    
    children: Dict[State, List[State]] = defaultdict(list)
    parent: Dict[State, Optional[State]] = {root: None}

    # build adjacency from the edges:
    out_adj: Dict[State, List[State]] = defaultdict(list)
    in_deg: Dict[State, int] = defaultdict(int)
    nodes: Set[State] = {root}

    for u, v in edges:
        out_adj[u].append(v)
        in_deg[v] += 1
        nodes.add(u)
        nodes.add(v)

    # bfs to assign parent:
    q = deque([root])
    visited: Set[State] = {root}

    while q:
        u = q.popleft()
        for v in out_adj.get(u, []):
            if v not in parent:
                parent[v] = u
                children[u].append(v)
            if v not in visited:
                visited.add(v)
                q.append(v)

    for n in nodes:
        children.setdefault(n, [])

    return RootedTree(root=root, children=dict(children), parent=parent)


# Newick serialization:
def to_newick(
    tree: RootedTree,
    leaf_label_fn: Optional[Callable[[State], str]] = None,
    internal_label_fn: Optional[Callable[[State], str]] = None,
    branch_length_fn: Optional[Callable[[State, State], Union[int, float]]] = None,
) -> str:
    if leaf_label_fn is None:
        leaf_label_fn = lambda st: state_to_str(st)
    if internal_label_fn is None:
        internal_label_fn = lambda st: ""
    if branch_length_fn is None:
        branch_length_fn = lambda u, v: None

    children = tree.children

    def rec(u: State) -> str:
        kids = children.get(u, [])
        if not kids:
            return leaf_label_fn(u)
        # internal node:
        sub = ",".join(
            _with_len(rec(v), branch_length_fn(u, v)) for v in kids
        )
        label = internal_label_fn(u)
        return f"({sub}){label}"

    def _with_len(s: str, bl) -> str:
        if bl is None:
            return s
        return f"{s}:{bl}"

    newick = rec(tree.root) + ";"
    return newick


# Greedy labels leaves as '{cell}_{state_str}' 
def greedy_style_leaf_labels(
    states_by_cell: Dict[str, State]
) -> Dict[State, List[str]]:
    m: Dict[State, List[str]] = defaultdict(list)
    for cell, st in states_by_cell.items():
        m[st].append(f"{cell}_{state_to_str(st)}")
    return dict(m)


# joins multiple cell labels
def build_leaf_label_fn_from_cells(
    states_by_cell: Dict[str, State]
) -> Callable[[State], str]:
    labels_by_state = greedy_style_leaf_labels(states_by_cell)

    def label_fn(st: State) -> str:
        labs = labels_by_state.get(st, [])
        if not labs:
            return state_to_str(st)
        if len(labs) == 1:
            return labs[0]
        return "|".join(labs)

    return label_fn
