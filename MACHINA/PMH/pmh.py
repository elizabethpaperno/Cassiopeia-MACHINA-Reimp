"""
pmh.py

Restricted and unrestricted Parsimonious Migration History (PMH), following the
MACHINA supplement:

- Section B.3.1  Unconstrained PMH via Sankoff (equation (11), Algorithms 1-3)
- Section B.3.2  Constrained PMH via integer linear programming (equations (12-32))

Pattern sets P:
    {PS}, {PS,S}, {PS,S,M}, {PS,S,M,R}

Input format:

Clone tree file: edge list, one edge per line: "parent child"
    A A1
    A A2
    A B
    B B1
    ...

Leaf labeling file: leaf to anatomical site labels, one per line:
    "leaf_vertex anatomical_site"
    A1 Om
    A2 SBwl
    B1 LOv
    ...

- All labeled vertices are leaves in the clone tree.
- The primary anatomical site (for example "LOv") is given as a label and
  must appear in the leaf labeling.
- The primary site is treated as site index zero, matching the convention
  in the MACHINA supplement where site one is primary.

Output layout (all under results/):

    results/PRIMARY_PATTERNSET/
        summary.txt
        labelings_all_mp.txt       # all maximum parsimony labelings with statistics
        labelings_opt.txt          # only (migration, comigration) optimal labelings
        site_graphs_all_mp/
            site_graph_0.txt
            site_graph_1.txt
            ...
        site_graphs_opt/
            site_graph_0.txt
            ...
"""

import argparse
from pathlib import Path
from collections import Counter
import pulp

# Base data structures
class CloneTree:
    """
    Rooted clone tree used for parsimonious migration history.

    Attributes
        vertices: List of vertex labels in the clone tree
        index_of: Map from vertex label to integer index
        children: Adjacency list of child indices for each vertex
        parent: Parent index of each vertex or None for the root
        root: Index of the root vertex
        leaves: Set of indices of leaf vertices
    """

    def __init__(self, edges):
        vertex_labels = set()
        for u, v in edges:
            vertex_labels.add(u)
            vertex_labels.add(v)

        self.vertices = sorted(vertex_labels)
        self.index_of = {lab: i for i, lab in enumerate(self.vertices)}
        n = len(self.vertices)

        self.children = [[] for _ in range(n)]
        self.parent = [None] * n

        # parent and children arrays
        for u_lab, v_lab in edges:
            u = self.index_of[u_lab]
            v = self.index_of[v_lab]
            self.children[u].append(v)
            if self.parent[v] is not None:
                raise ValueError(
                    "Vertex {} has more than one parent in the clone tree".format(v_lab)
                )
            self.parent[v] = u

        # root
        roots = [i for i, p in enumerate(self.parent) if p is None]
        if len(roots) != 1:
            raise ValueError(
                "Expected exactly one root, but found {}".format(len(roots))
            )
        self.root = roots[0]

        # leaves
        self.leaves = {i for i in range(n) if not self.children[i]}

class SiteIndex:
    """
    Indexing of anatomical sites with primary site at index zero

    Attributes:
        labels: List of site labels ordered with primary first
        index_of: Map from site label to integer index
        primary_index: Index of the primary site (zero)
        m: Number of anatomical sites
    """

    def __init__(self, site_labels, primary_site_label):
        unique_sites = sorted(set(site_labels))
        if primary_site_label not in unique_sites:
            raise ValueError(
                "Primary site {} is not among sites {}".format(
                    primary_site_label, unique_sites
                )
            )
        ordered = [primary_site_label] + [
            s for s in unique_sites if s != primary_site_label
        ]
        self.labels = ordered
        self.index_of = {lab: i for i, lab in enumerate(self.labels)}
        self.primary_index = 0
        self.m = len(self.labels)

# Input files parsing
def read_clone_tree(path):
    """
    Read a clone tree edge list into a CloneTree

    Arguments:
        path: Path to a file of parent child edges

    Returns
        tree: CloneTree built from the edge list
    """
    edges = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise ValueError("Invalid edge line: {}".format(line))
            edges.append((parts[0], parts[1]))
    return CloneTree(edges)


def read_leaf_labeling(path):
    """
    Read leaf to anatomical site labels from file

    Arguments:
        path: Path to a file of leaf and site label pairs

    Returns
        labeling: Map from leaf label to anatomical site label
    """
    labeling = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise ValueError("Invalid labeling line: {}".format(line))
            leaf, site = parts
            labeling[leaf] = site
    return labeling

# Modified Sankoff (unconstrained PMH)

# Algorithm 1 in Supplementary Info
def sankoff(tree, sites, leaf_labeling):
    """
    Finds vertex labelings of tree with the minimum migration number

    Arguments:
        tree: CloneTree to label
        sites: SiteIndex describing anatomical sites
        leaf_labeling: Map from leaf label to site label

    Returns
        all_labelings: List of all maximum parsimony labelings
        mu_star: Minimum migration number for the tree
        M: Dynamic program cost table
        Delta: Dynamic program choice table
    """
    M, Delta, mu_star = solve(tree, sites, leaf_labeling)
    all_labelings = backtrace(tree, sites, leaf_labeling, M, Delta)
    return all_labelings, mu_star, M, Delta


# Algorithm 2 in Supplementary Info
def solve(tree, sites, leaf_labeling):
    """
    Run the Sankoff dynamic program for unconstrained PMH

    Arguments:
        tree: CloneTree
        sites: SiteIndex describing anatomical sites
        leaf_labeling: Map from leaf label to site label

    Returns
        M: Minimal migration cost table for each vertex and site
        Delta: Allowed child site choices for each vertex and site
        mu_star: Optimal migration number
    """
    m = sites.m

    leaf_site = {}
    for leaf_label, site_label in leaf_labeling.items():
        if leaf_label not in tree.index_of:
            raise ValueError(
                "Labeled vertex {} is not present in the clone tree".format(leaf_label)
            )
        u = tree.index_of[leaf_label]
        if u not in tree.leaves:
            raise ValueError(
                "Labeled vertex {} is not a leaf in the clone tree".format(leaf_label)
            )
        leaf_site[u] = sites.index_of[site_label]

    big_number = 10 ** 9
    M = {}
    Delta = {}

    # post order traversal
    order = []

    def dfs(u):
        for v in tree.children[u]:
            dfs(v)
        order.append(u)

    dfs(tree.root)

    for u in order:
        if u in tree.leaves:
            # leaf
            for s in range(m):
                if leaf_site.get(u, None) == s:
                    M[(u, s)] = 0
                else:
                    M[(u, s)] = big_number
                Delta[(u, s)] = {}
        else:
            # internal node
            for s in range(m):
                total_cost = 0
                child_choices = {}
                for v in tree.children[u]:
                    best_cost = big_number
                    best_sites = []
                    for t in range(m):
                        extra = 0 if s == t else 1
                        cost = extra + M[(v, t)]
                        if cost < best_cost:
                            best_cost = cost
                            best_sites = [t]
                        elif cost == best_cost:
                            best_sites.append(t)
                    total_cost += best_cost
                    child_choices[v] = best_sites
                M[(u, s)] = total_cost
                Delta[(u, s)] = child_choices

    mu_star = M[(tree.root, sites.primary_index)]
    return M, Delta, mu_star

# Algorithm 3 in Supplementary Info
def backtrace(tree, sites, leaf_labeling, M, Delta):
    """
    Enumerate all minimum migration labelings by backtrace

    Arguments:
        tree: CloneTree
        sites: SiteIndex describing anatomical sites
        M: Minimal migration cost table
        Delta: Child site choice table

    Returns
        labelings: List of vertex to site index labelings with mu equal to mu star
    """
    root = tree.root
    primary_index = sites.primary_index

    memo = {}

    def enum_subtree(u, s):
        key = (u, s)
        if key in memo:
            return memo[key]

        if not tree.children[u]:
            result = [{u: s}]
            memo[key] = result
            return result

        child_labelings_per_child = []
        for v in tree.children[u]:
            labelings_for_child = []
            for t in Delta[(u, s)][v]:
                for sublab in enum_subtree(v, t):
                    labelings_for_child.append(sublab)
            child_labelings_per_child.append(labelings_for_child)

        import itertools

        combined = []
        for combo in itertools.product(*child_labelings_per_child):
            merged = {u: s}
            for lab in combo:
                merged.update(lab)
            combined.append(merged)
        return combined
    
    return enum_subtree(root, primary_index)

# ILP for constrained PMH 
# Modeled after Section B.3.2 in MACHINA Supplementary Information
def solve_pmh_ilp_mode(tree, sites, leaf_labeling, mode):
    """
    Solve constrained PMH with a single ILP mode

    Arguments:
        tree: CloneTree to label
        sites: SiteIndex describing anatomical sites
        leaf_labeling: Map from leaf label to site label
        mode: One of "PS", "S", or "M" selecting the pattern

    Returns
        labeling: Optimal vertex to site index labeling
        objective_value: Value of the ILP objective at optimum
        mu: Migration number of the solution
        phi: Comigration number of the solution
        sigma: Seeding site number of the solution
    """
    if pulp is None:
        raise ImportError(
            "The pulp package is required for ILP modes"
            "Install it with: pip install pulp"
        )

    if mode not in {"PS", "S", "M"}:
        raise ValueError("ILP mode must be 'PS', 'S', or 'M', got {}".format(mode))

    n = len(tree.vertices)
    m = sites.m
    root = tree.root
    primary_index = sites.primary_index

    # Map labeled leaves to site indices
    leaf_site = {}
    for leaf_label, site_label in leaf_labeling.items():
        if leaf_label not in tree.index_of:
            raise ValueError(
                "Labeled vertex {} is not in the clone tree".format(leaf_label)
            )
        u = tree.index_of[leaf_label]
        if u not in tree.leaves:
            raise ValueError(
                "Labeled vertex {} is not a leaf".format(leaf_label)
            )
        leaf_site[u] = sites.index_of[site_label]

    # Edge list of the tree
    edges = []
    for i in range(n):
        for j in tree.children[i]:
            edges.append((i, j))

    model = pulp.LpProblem("PMH_ILP", pulp.LpMinimize)

    # Variables x[i, s] in {0, 1} indicate vertex i uses site s
    x = pulp.LpVariable.dicts(
        "x",
        ((i, s) for i in range(n) for s in range(m)),
        lowBound=0,
        upBound=1,
        cat=pulp.LpBinary,
    )

    # Variables z[i, j, s] in {0, 1} for edge (i, j) and site s
    z = pulp.LpVariable.dicts(
        "z",
        ((i, j, s) for (i, j) in edges for s in range(m)),
        lowBound=0,
        upBound=1,
        cat=pulp.LpBinary,
    )

    # Variables y[i, j] in {0, 1} indicate whether edge (i, j) is a migration
    y = pulp.LpVariable.dicts(
        "y",
        ((i, j) for (i, j) in edges),
        lowBound=0,
        upBound=1,
        cat=pulp.LpBinary,
    )

    # Variables c[s, t] in {0, 1} indicate whether there is at least one
    # migration edge from site s to site t
    c = pulp.LpVariable.dicts(
        "c",
        ((s, t) for s in range(m) for t in range(m) if s != t),
        lowBound=0,
        upBound=1,
        cat=pulp.LpBinary,
    )

    # Variables d[s] in {0, 1} indicate whether site s is a seeding site
    d = pulp.LpVariable.dicts(
        "d",
        (s for s in range(m)),
        lowBound=0,
        upBound=1,
        cat=pulp.LpBinary,
    )

    # Equation (12): each vertex gets exactly one site
    for i in range(n):
        model += pulp.lpSum(x[(i, s)] for s in range(m)) == 1, "vertex_one_site_{}".format(i)

    # Equation (13): root is labeled by primary site
    model += x[(root, primary_index)] == 1, "root_primary"

    # Equations (14) and (15): leaves have fixed site labels
    for i in tree.leaves:
        if i not in leaf_site:
            raise ValueError(
                "Leaf {} is missing a label in the leaf labeling file".format(
                    tree.vertices[i]
                )
            )
        fixed_site = leaf_site[i]
        for s in range(m):
            if s == fixed_site:
                model += x[(i, s)] == 1, "leaf_{}_site_{}_fixed_one".format(i, s)
            else:
                model += x[(i, s)] == 0, "leaf_{}_site_{}_fixed_zero".format(i, s)

    # Equations (16) and (17): z[i, j, s] is only one if both endpoints use site s
    for (i, j) in edges:
        for s in range(m):
            model += z[(i, j, s)] <= x[(i, s)], "z_leq_xi_{}_{}_{}".format(i, j, s)
            model += z[(i, j, s)] <= x[(j, s)], "z_leq_xj_{}_{}_{}".format(i, j, s)

    # Equation (18): sum over s of z[i, j, s] equals one minus y[i, j]
    for (i, j) in edges:
        model += (
            pulp.lpSum(z[(i, j, s)] for s in range(m)) == 1 - y[(i, j)]
        ), "z_sum_y_{}_{}".format(i, j)

    # Equation (19): c[s, t] is at least x[i, s] plus x[j, t] minus one
    for (i, j) in edges:
        for s in range(m):
            for t in range(m):
                if s == t:
                    continue
                model += (
                    c[(s, t)] >= x[(i, s)] + x[(j, t)] - 1
                ), "c_ge_xs_xt_{}_{}_{}_{}".format(i, j, s, t)

    # Equation (21): d[s] is at least c[s, t] for all t not equal to s
    for s in range(m):
        for t in range(m):
            if s == t:
                continue
            model += d[s] >= c[(s, t)], "d_ge_c_{}_{}".format(s, t)

    # Equation (22): primary is forced to be a seeding site
    model += d[primary_index] == 1, "primary_seeding"

    # Topological constraints by mode (equations (28) through (32))
    if mode == "PS":
        # Equation (28): no incoming edges to primary site
        model += pulp.lpSum(
            c[(s, primary_index)] for s in range(m) if s != primary_index
        ) == 0, "ps_no_incoming_primary"

        # Equation (29): no metastasis to metastasis edges
        for s in range(m):
            if s == primary_index:
                continue
            for t in range(m):
                if t == primary_index or t == s:
                    continue
                model += c[(s, t)] == 0, "ps_no_met_to_met_{}_{}".format(s, t)

        # Equation (30): primary site seeds every metastasis site
        for t in range(m):
            if t == primary_index:
                continue
            model += c[(primary_index, t)] == 1, "ps_primary_to_{}".format(t)

    elif mode == "S":
        # Single source seeding mode uses equation (28) and (31)

        # Equation (28): no incoming edges to primary site
        model += pulp.lpSum(
            c[(s, primary_index)] for s in range(m) if s != primary_index
        ) == 0, "s_no_incoming_primary"

        # Equation (31): for each non primary site t, there is exactly
        # one incoming edge from some site s
        for t in range(m):
            if t == primary_index:
                continue
            model += (
                pulp.lpSum(c[(s, t)] for s in range(m) if s != t) == 1
            ), "s_indegree_one_{}".format(t)

    elif mode == "M":
        # Multi source seeding mode with directed acyclic site graph
        # uses cycle inequalities (equation (32)).

        cycles = []

        def generate_cycles():
            for start in range(m):
                stack = [start]
                visited = {start}

                def dfs(curr):
                    for nxt in range(m):
                        if nxt == curr:
                            continue
                        if nxt == start:
                            # Found a cycle with length at least two
                            if len(stack) >= 2:
                                cycle = stack + [start]
                                if start == min(cycle):
                                    cycles.append(cycle)
                            continue
                        if nxt in visited:
                            continue
                        visited.add(nxt)
                        stack.append(nxt)
                        dfs(nxt)
                        stack.pop()
                        visited.remove(nxt)

                dfs(start)

        generate_cycles()

        for cycle in cycles:
            edges_cycle = []
            for i in range(len(cycle) - 1):
                s = cycle[i]
                t = cycle[i + 1]
                if s != t:
                    edges_cycle.append((s, t))
            if not edges_cycle:
                continue
            model += (
                pulp.lpSum(c[(s, t)] for (s, t) in edges_cycle)
                <= len(edges_cycle) - 1
            ), "dag_no_cycle_{}".format("_".join(str(x) for x in cycle))

    # Objective fxn
    model += (
        pulp.lpSum(y[(i, j)] for (i, j) in edges)
        + (1.0 / n)
        * pulp.lpSum(
            c[(s, t)] for s in range(m) for t in range(m) if s != t
        )
        + (1.0 / (m * n)) * pulp.lpSum(d[s] for s in range(m))
    )

    # Solve the ILP
    model.solve(pulp.PULP_CBC_CMD(msg=False))
    status = pulp.LpStatus[model.status]
    if status != "Optimal":
        raise RuntimeError("ILP status {}, not optimal".format(status))

    # Extract the vertex labeling from x[i, s]
    labeling = {}
    for i in range(n):
        assigned = None
        for s in range(m):
            val = pulp.value(x[(i, s)])
            if val is not None and val > 0.5:
                assigned = s
                break
        if assigned is None:
            raise RuntimeError(
                "No site chosen for vertex {}".format(tree.vertices[i])
            )
        labeling[i] = assigned

    mu, phi, sigma = compute_migration_stats(tree, sites, labeling)
    objective_value = pulp.value(model.objective)

    return labeling, objective_value, mu, phi, sigma


def compute_migration_stats(tree, sites, labeling):
    """
    Compute migration, comigration, and seeding number

    Arguments:
        tree: CloneTree with fixed vertex labeling
        sites: SiteIndex describing anatomical sites
        labeling: Map from vertex index to site index

    Returns
        mu: Number of migration edges in the tree
        phi: Comigration number over all ordered site pairs
        sigma: Number of seeding sites with outgoing migrations
    """
    n = len(tree.vertices)
    m = sites.m
    root = tree.root

    # Group migration edges by ordered site pair
    edges_by_pair = {}
    for s in range(m):
        for t in range(m):
            if s != t:
                edges_by_pair[(s, t)] = []

    migration_edges = []

    for u in range(n):
        for v in tree.children[u]:
            su = labeling[u]
            sv = labeling[v]
            if su != sv:
                migration_edges.append((u, v, su, sv))
                edges_by_pair[(su, sv)].append((u, v))

    mu = len(migration_edges)

    seeding_sites = set()
    for (u, v, su, sv) in migration_edges:
        seeding_sites.add(su)
    sigma = len(seeding_sites)

    # comigration number phi
    phi = 0
    children = tree.children

    for pair, edges in edges_by_pair.items():
        if not edges:
            continue

        edge_set = set(edges)
        max_on_path = 0

        # dfs from the root
        stack = [(root, 0)]
        while stack:
            u, count_so_far = stack.pop()
            if count_so_far > max_on_path:
                max_on_path = count_so_far
            for v in children[u]:
                extra = 1 if (u, v) in edge_set else 0
                stack.append((v, count_so_far + extra))

        phi += max_on_path

    return mu, phi, sigma

# Constructs migration graph from a labeling
def site_graph_from_labeling(tree, labeling):
    """
    Build a site migration graph

    Arguments:
        tree: CloneTree whose edges are examined
        labeling: Map from vertex index to site index

    Returns
        counts: Counts num of migration edges between site pairs
    """
    counts = Counter()
    n = len(tree.vertices)
    for u in range(n):
        for v in tree.children[u]:
            su = labeling[u]
            sv = labeling[v]
            if su != sv:
                counts[(su, sv)] += 1
    return counts


def write_site_graph(path, sites, edge_counts):
    """
    Outputs site graph as a tab separated text file

    Arguments:
        path: Path for the file
        sites: SiteIndex
        edge_counts: Site pair migration counts

    """
    with open(path, "w") as f:
        f.write("source_site\ttarget_site\tmultiplicity\n")
        for (s, t), mult in sorted(edge_counts.items(),
                                   key=lambda x: (x[0][0], x[0][1])):
            src = sites.labels[s]
            dst = sites.labels[t]
            f.write("{}\t{}\t{}\n".format(src, dst, mult))


# Helper method
def dedup_labelings(labelings):
    seen = set()
    unique = []
    for lab in labelings:
        key = tuple(sorted(lab.items()))
        if key not in seen:
            seen.add(key)
            unique.append(lab)
    return unique


# main function to solve pmh given a mattern set
def solve_pmh_with_pattern_set(tree, sites, leaf_labeling, pattern_set_str):
    """
    Solve PMH for a chosen pattern set using Sankoff or ILP

    Arguments:
        tree: CloneTree to label
        sites: SiteIndex describing anatomical sites
        leaf_labeling: Map from leaf label to site label
        pattern_set_str: String describing the allowed pattern set

    Returns
        labelings: List of vertex labelings consistent with the pattern set
        mu_star: Minimum migration number over the pattern set
        phi_star: Minimum comigration number at mu star
        sigma_star: Minimum seeding site number at mu and phi star
    """
    normalized = pattern_set_str.replace(" ", "")
    valid_sets = {"PS", "PS,S", "PS,S,M", "PS,S,M,R", "S,M,R"}
    if normalized not in valid_sets:
        raise ValueError(
            "Pattern set must be one of: 'PS', 'PS,S', 'PS,S,M', "
            "'PS,S,M,R' (or 'S,M,R' as an alias for unrestricted)"
        )

    # unrestricted PMH -> use modified Sankoff Algorithm
    if normalized in {"PS,S,M,R", "S,M,R"}:
        all_labelings, mu_star, M, Delta = sankoff(tree, sites, leaf_labeling)

        # find optimal among all maximum parsimony labelings
        best_phi = None
        best_sigma = None
        for lab in all_labelings:
            mu, phi, sigma = compute_migration_stats(tree, sites, lab)
            if mu != mu_star:
                raise RuntimeError(
                    "Sankoff enumeration error: labeling has mu = {} "
                    "but mu_star = {}".format(mu, mu_star)
                )
            if (
                best_phi is None
                or phi < best_phi
                or (phi == best_phi and sigma < best_sigma)
            ):
                best_phi = phi
                best_sigma = sigma

        all_labelings = dedup_labelings(all_labelings)
        return all_labelings, mu_star, best_phi, best_sigma

    # restricted pattern sets -> use ILP modes
    if normalized == "PS":
        lab, obj, mu, phi, sigma = solve_pmh_ilp_mode(
            tree, sites, leaf_labeling, mode="PS"
        )
        labelings = dedup_labelings([lab])
        return labelings, mu, phi, sigma

    if normalized == "PS,S":
        lab_PS, obj_PS, mu_PS, phi_PS, sigma_PS = solve_pmh_ilp_mode(
            tree, sites, leaf_labeling, mode="PS"
        )
        lab_S, obj_S, mu_S, phi_S, sigma_S = solve_pmh_ilp_mode(
            tree, sites, leaf_labeling, mode="S"
        )

        key_PS = (mu_PS, phi_PS, sigma_PS)
        key_S = (mu_S, phi_S, sigma_S)

        if key_PS < key_S:
            labelings = [lab_PS]
            mu_star, phi_star, sigma_star = key_PS
        elif key_S < key_PS:
            labelings = [lab_S]
            mu_star, phi_star, sigma_star = key_S
        else:
            labelings = [lab_PS, lab_S]
            mu_star, phi_star, sigma_star = key_PS

        labelings = dedup_labelings(labelings)
        return labelings, mu_star, phi_star, sigma_star

    if normalized == "PS,S,M":
        lab_PS, obj_PS, mu_PS, phi_PS, sigma_PS = solve_pmh_ilp_mode(
            tree, sites, leaf_labeling, mode="PS"
        )
        lab_S, obj_S, mu_S, phi_S, sigma_S = solve_pmh_ilp_mode(
            tree, sites, leaf_labeling, mode="S"
        )
        lab_M, obj_M, mu_M, phi_M, sigma_M = solve_pmh_ilp_mode(
            tree, sites, leaf_labeling, mode="M"
        )

        key_PS = (mu_PS, phi_PS, sigma_PS)
        key_S = (mu_S, phi_S, sigma_S)
        key_M = (mu_M, phi_M, sigma_M)

        best_key = min(key_PS, key_S, key_M)

        labelings = []
        if key_PS == best_key:
            labelings.append(lab_PS)
        if key_S == best_key:
            labelings.append(lab_S)
        if key_M == best_key:
            labelings.append(lab_M)

        labelings = dedup_labelings(labelings)
        mu_star, phi_star, sigma_star = best_key
        return labelings, mu_star, phi_star, sigma_star

def main():
    parser = argparse.ArgumentParser(
        description=(
            "PMH implementation with pattern sets "
            "{PS}, {PS,S}, {PS,S,M}, {PS,S,M,R} "
            "following MACHINA's ILP and Sankoff formulations."
        )
    )
    parser.add_argument("--tree", required=True, help="Clone tree edge list file")
    parser.add_argument("--labels", required=True, help="Leaf labeling file")
    parser.add_argument(
        "--primary", required=True, help="Primary anatomical site label"
    )
    parser.add_argument(
        "--pattern-set",
        required=True,
        help="One of: 'PS', 'PS,S', 'PS,S,M', 'PS,S,M,R'",
    )

    args = parser.parse_args()

    # Load inputs
    tree = read_clone_tree(args.tree)
    leaf_labeling = read_leaf_labeling(args.labels)
    all_sites = set(leaf_labeling.values())
    sites = SiteIndex(all_sites, primary_site_label=args.primary)

    # Solve PMH for the chosen pattern set
    labelings, mu_star, phi_star, sigma_star = solve_pmh_with_pattern_set(
        tree, sites, leaf_labeling, args.pattern_set
    )

    normalized = args.pattern_set.replace(" ", "")
    is_unrestricted = normalized in {"PS,S,M,R", "S,M,R"}

    # Compute statistics for each labeling
    labeling_stats = []
    for idx, lab in enumerate(labelings):
        mu, phi, sigma = compute_migration_stats(tree, sites, lab)
        labeling_stats.append((idx, lab, mu, phi, sigma))

    num_mp = len(labelings)

    # Find optimal labelings from max parsimony labels
    if is_unrestricted:
        phi_opt = None
        sigma_opt = None
        for _, _, mu, phi, sigma in labeling_stats:
            if mu != mu_star:
                raise RuntimeError(
                    "Sankoff enumeration error: labeling has mu = {} "
                    "but mu_star = {}".format(mu, mu_star)
                )
            if (
                phi_opt is None
                or phi < phi_opt
                or (phi == phi_opt and sigma < sigma_opt)
            ):
                phi_opt = phi
                sigma_opt = sigma
        phi_star, sigma_star = phi_opt, sigma_opt
        opt_labelings = [
            (idx, lab, mu, phi, sigma)
            for (idx, lab, mu, phi, sigma) in labeling_stats
            if phi == phi_star and sigma == sigma_star
        ]
        num_opt = len(opt_labelings)
    else:
        opt_labelings = labeling_stats
        num_opt = num_mp

    # Save output files to results dir
    base_dir = Path("results")
    pattern_dir = base_dir / "{}_{}".format(args.primary, normalized)
    pattern_dir.mkdir(parents=True, exist_ok=True)

    # Summary file
    summary_path = pattern_dir / "summary.txt"
    with open(summary_path, "w") as f:
        f.write("pattern_set\t{}\n".format(normalized))
        f.write("primary_site\t{}\n".format(args.primary))
        f.write("num_sites\t{}\n".format(sites.m))
        f.write("num_vertices\t{}\n".format(len(tree.vertices)))
        f.write("mu_opt\t{}\n".format(mu_star))
        f.write("phi_opt\t{}\n".format(phi_star))
        f.write("sigma_opt\t{}\n".format(sigma_star))
        f.write("n_mp_labelings\t{}\n".format(num_mp))
        f.write("n_opt_labelings\t{}\n".format(num_opt))

    # Maximum parsimony labelings
    all_mp_path = pattern_dir / "labelings_all_mp.txt"
    with open(all_mp_path, "w") as f:
        for idx, lab, mu, phi, sigma in labeling_stats:
            f.write(
                "# labeling {}\tmu={}\tphi={}\tsigma={}\n".format(
                    idx, mu, phi, sigma
                )
            )
            for i in sorted(lab.keys(), key=lambda x: tree.vertices[x]):
                v = tree.vertices[i]
                site_index = lab[i]
                f.write("{}\t{}\n".format(v, sites.labels[site_index]))
            f.write("\n")

    # Optimal labeling output file
    opt_path = pattern_dir / "labelings_opt.txt"
    with open(opt_path, "w") as f:
        for idx, lab, mu, phi, sigma in opt_labelings:
            f.write(
                "# labeling {}\tmu={}\tphi={}\tsigma={}\n".format(
                    idx, mu, phi, sigma
                )
            )
            for i in sorted(lab.keys(), key=lambda x: tree.vertices[x]):
                v = tree.vertices[i]
                site_index = lab[i]
                f.write("{}\t{}\n".format(v, sites.labels[site_index]))
            f.write("\n")

    # Maximum parsimony graph output file
    sg_all_dir = pattern_dir / "site_graphs_all_mp"
    sg_all_dir.mkdir(exist_ok=True)
    for idx, lab, mu, phi, sigma in labeling_stats:
        edges = site_graph_from_labeling(tree, lab)
        sg_path = sg_all_dir / "site_graph_{}.txt".format(idx)
        write_site_graph(sg_path, sites, edges)

    # Opt site graph output file
    sg_opt_dir = pattern_dir / "site_graphs_opt"
    sg_opt_dir.mkdir(exist_ok=True)
    for idx, lab, mu, phi, sigma in opt_labelings:
        edges = site_graph_from_labeling(tree, lab)
        sg_path = sg_opt_dir / "site_graph_{}.txt".format(idx)
        write_site_graph(sg_path, sites, edges)

    # Printed to terminal
    print("P = {}, primary = {}".format(normalized, args.primary))
    print(
        "mu_star = {}, phi_opt = {}, sigma_opt = {}".format(
            mu_star, phi_star, sigma_star
        )
    )
    print(
        "Number of maximum parsimony labelings = {}, "
        "number of optimal labelings = {}".format(num_mp, num_opt)
    )
    print("Results written to: {}".format(pattern_dir))


if __name__ == "__main__":
    main()
