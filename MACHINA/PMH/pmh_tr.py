"""
pmh_tr.py

Parsimonious Migration History with Tree Resolution (PMH-TR), following MACHINA:

Key Innovation: Polytomy resolution is integrated directly into the ILP formulation
rather than enumerating all binary resolutions. Internal nodes are labeled with a 
single anatomical site, which implicitly optimizes the tree structure.

Pattern sets P:
    {PS}, {PS,S}, {PS,S,M}, {PS,S,M,R}
"""

import argparse
from pathlib import Path
from collections import Counter
import pulp


class CloneTree:
    """Rooted clone tree supporting polytomies (nodes with >2 children)."""

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

        for u_lab, v_lab in edges:
            u = self.index_of[u_lab]
            v = self.index_of[v_lab]
            self.children[u].append(v)
            if self.parent[v] is not None:
                raise ValueError(f"Vertex {v_lab} has multiple parents")
            self.parent[v] = u

        roots = [i for i, p in enumerate(self.parent) if p is None]
        if len(roots) != 1:
            raise ValueError(f"Expected 1 root, found {len(roots)}")
        self.root = roots[0]

        self.leaves = {i for i in range(n) if not self.children[i]}

    def has_polytomies(self):
        return any(len(ch) > 2 for ch in self.children)

    def get_polytomy_nodes(self):
        return [i for i, ch in enumerate(self.children) if len(ch) > 2]


class SiteIndex:
    """Indexing of anatomical sites with primary site at index zero."""

    def __init__(self, site_labels, primary_site_label):
        unique_sites = sorted(set(site_labels))
        if primary_site_label not in unique_sites:
            raise ValueError(f"Primary site {primary_site_label} not in {unique_sites}")
        ordered = [primary_site_label] + [s for s in unique_sites if s != primary_site_label]
        self.labels = ordered
        self.index_of = {lab: i for i, lab in enumerate(self.labels)}
        self.primary_index = 0
        self.m = len(self.labels)


def read_clone_tree(path):
    """Read clone tree edge list."""
    edges = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f"Invalid line: {line}")
            edges.append((parts[0], parts[1]))
    return CloneTree(edges)


def read_leaf_labeling(path):
    """Read leaf to site labels."""
    labeling = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f"Invalid line: {line}")
            labeling[parts[0]] = parts[1]
    return labeling


def compute_migration_stats(tree, sites, labeling):
    """Compute mu (migrations), phi (comigrations), sigma (seeding sites)."""
    n = len(tree.vertices)
    m = sites.m
    root = tree.root

    edges_by_pair = {(s, t): [] for s in range(m) for t in range(m) if s != t}
    migration_edges = []

    for u in range(n):
        for v in tree.children[u]:
            su = labeling[u]
            sv = labeling[v]
            if su != sv:
                migration_edges.append((u, v, su, sv))
                edges_by_pair[(su, sv)].append((u, v))

    mu = len(migration_edges)
    sigma = len(set(su for u, v, su, sv in migration_edges))

    # Comigration number
    phi = 0
    for pair, edges in edges_by_pair.items():
        if not edges:
            continue
        edge_set = set(edges)
        max_on_path = 0
        stack = [(root, 0)]
        while stack:
            u, count = stack.pop()
            max_on_path = max(max_on_path, count)
            for v in tree.children[u]:
                extra = 1 if (u, v) in edge_set else 0
                stack.append((v, count + extra))
        phi += max_on_path

    return mu, phi, sigma


def solve_pmh_tr_ilp_mode(tree, sites, leaf_labeling, mode):
    """
    Solve PMH-TR with tree resolution integrated into ILP.
    
    The key insight is that we label all vertices (including internal ones)
    with a single anatomical site, and the ILP automatically finds the labeling
    that minimizes migrations. No need to explicitly enumerate tree refinements.
    """
    if pulp is None:
        raise ImportError("pulp required: pip install pulp")

    if mode not in {"PS", "S", "M"}:
        raise ValueError(f"Invalid mode: {mode}")

    n = len(tree.vertices)
    m = sites.m
    root = tree.root
    primary_index = sites.primary_index

    # Map leaf labels to site indices
    leaf_site = {}
    for leaf_label, site_label in leaf_labeling.items():
        if leaf_label not in tree.index_of:
            continue
        u = tree.index_of[leaf_label]
        if u not in tree.leaves:
            continue
        leaf_site[u] = sites.index_of[site_label]

    edges = [(i, j) for i in range(n) for j in tree.children[i]]

    model = pulp.LpProblem("PMH_TR_ILP", pulp.LpMinimize)

    # x[i, s]: vertex i labeled with site s
    x = pulp.LpVariable.dicts(
        "x",
        ((i, s) for i in range(n) for s in range(m)),
        lowBound=0,
        upBound=1,
        cat=pulp.LpBinary,
    )

    # z[i, j, s]: edge (i,j) has no migration at site s
    z = pulp.LpVariable.dicts(
        "z",
        ((i, j, s) for (i, j) in edges for s in range(m)),
        lowBound=0,
        upBound=1,
        cat=pulp.LpBinary,
    )

    # y[i, j]: edge (i,j) is a migration
    y = pulp.LpVariable.dicts(
        "y",
        ((i, j) for (i, j) in edges),
        lowBound=0,
        upBound=1,
        cat=pulp.LpBinary,
    )

    # c[s, t]: migration exists from site s to site t
    c = pulp.LpVariable.dicts(
        "c",
        ((s, t) for s in range(m) for t in range(m) if s != t),
        lowBound=0,
        upBound=1,
        cat=pulp.LpBinary,
    )

    # d[s]: site s is seeding
    d = pulp.LpVariable.dicts(
        "d",
        (s for s in range(m)),
        lowBound=0,
        upBound=1,
        cat=pulp.LpBinary,
    )

    # Each vertex gets exactly one site
    for i in range(n):
        model += pulp.lpSum(x[(i, s)] for s in range(m)) == 1, f"site_{i}"

    # Root labeled by primary site
    model += x[(root, primary_index)] == 1, "root_primary"

    # Leaves with known labels are fixed
    for leaf_label, site_label in leaf_labeling.items():
        if leaf_label not in tree.index_of:
            continue
        i = tree.index_of[leaf_label]
        if i not in tree.leaves:
            continue
        fixed_site = sites.index_of[site_label]
        for s in range(m):
            model += x[(i, s)] == (1 if s == fixed_site else 0), f"leaf_{i}_{s}"

    # Edge constraints
    for (i, j) in edges:
        for s in range(m):
            model += z[(i, j, s)] <= x[(i, s)], f"z_i_{i}_{j}_{s}"
            model += z[(i, j, s)] <= x[(j, s)], f"z_j_{i}_{j}_{s}"

    # Migration edges: sum(z) + y[i,j] = 1
    for (i, j) in edges:
        model += pulp.lpSum(z[(i, j, s)] for s in range(m)) + y[(i, j)] == 1, f"migrate_{i}_{j}"

    # Migration graph edges
    for (i, j) in edges:
        for s in range(m):
            for t in range(m):
                if s != t:
                    model += c[(s, t)] >= x[(i, s)] + x[(j, t)] - 1, f"c_{i}_{j}_{s}_{t}"

    # Seeding constraints
    for s in range(m):
        for t in range(m):
            if s != t:
                model += d[s] >= c[(s, t)], f"seed_{s}_{t}"

    model += d[primary_index] == 1, "primary_is_seed"

    # Mode-specific topological constraints
    if mode == "PS":
        # No incoming to primary
        model += pulp.lpSum(
            c[(s, primary_index)] for s in range(m) if s != primary_index
        ) == 0, "ps_no_in"

        # No metastasis-to-metastasis
        for s in range(m):
            if s != primary_index:
                for t in range(m):
                    if t != primary_index and t != s:
                        model += c[(s, t)] == 0, f"ps_no_met_{s}_{t}"

        # Primary seeds all metastases
        for t in range(m):
            if t != primary_index:
                model += c[(primary_index, t)] == 1, f"ps_seed_{t}"

    elif mode == "S":
        # No incoming to primary
        model += pulp.lpSum(
            c[(s, primary_index)] for s in range(m) if s != primary_index
        ) == 0, "s_no_in"

        # Each non-primary has one incoming
        for t in range(m):
            if t != primary_index:
                model += pulp.lpSum(
                    c[(s, t)] for s in range(m) if s != t
                ) == 1, f"s_in_{t}"

    elif mode == "M":
        # No incoming to primary (DAG requirement)
        model += pulp.lpSum(
            c[(s, primary_index)] for s in range(m) if s != primary_index
        ) == 0, "m_no_in"

    # Objective: minimize migrations, then comigrations, then seeding
    model += (
        pulp.lpSum(y[(i, j)] for (i, j) in edges)
        + (1.0 / n) * pulp.lpSum(c[(s, t)] for s in range(m) for t in range(m) if s != t)
        + (1.0 / (m * n)) * pulp.lpSum(d[s] for s in range(m))
    )

    model.solve(pulp.PULP_CBC_CMD(msg=False))
    status = pulp.LpStatus[model.status]
    if status != "Optimal":
        raise RuntimeError(f"ILP status: {status}")

    labeling = {}
    for i in range(n):
        for s in range(m):
            if pulp.value(x[(i, s)]) and pulp.value(x[(i, s)]) > 0.5:
                labeling[i] = s
                break
        if i not in labeling:
            raise RuntimeError(f"No site for {tree.vertices[i]}")

    mu, phi, sigma = compute_migration_stats(tree, sites, labeling)
    return labeling, pulp.value(model.objective), mu, phi, sigma


def site_graph_from_labeling(tree, labeling):
    """Build migration graph."""
    counts = Counter()
    for u in range(len(tree.vertices)):
        for v in tree.children[u]:
            if labeling[u] != labeling[v]:
                counts[(labeling[u], labeling[v])] += 1
    return counts


def write_site_graph(path, sites, edge_counts):
    """Write site graph."""
    with open(path, "w") as f:
        for (s, t), mult in sorted(edge_counts.items()):
            f.write(f"{sites.labels[s]}\t{sites.labels[t]}\t{mult}\n")


def solve_pmh_tr_with_pattern_set(tree, sites, leaf_labeling, pattern_set_str):
    """Solve PMH-TR for a pattern set."""
    normalized = pattern_set_str.replace(" ", "")
    valid = {"PS", "PS,S", "PS,S,M", "PS,S,M,R"}
    if normalized not in valid:
        raise ValueError(f"Invalid pattern set: {normalized}")

    # Determine modes
    modes = {
        "PS": ["PS"],
        "PS,S": ["PS", "S"],
        "PS,S,M": ["PS", "S", "M"],
        "PS,S,M,R": ["PS", "S", "M"],
    }[normalized]

    best_key = None
    best_solutions = []

    for mode in modes:
        labeling, obj, mu, phi, sigma = solve_pmh_tr_ilp_mode(
            tree, sites, leaf_labeling, mode
        )
        current_key = (mu, phi, sigma)

        if best_key is None or current_key < best_key:
            best_key = current_key
            best_solutions = [(labeling, mu, phi, sigma)]
        elif current_key == best_key:
            best_solutions.append((labeling, mu, phi, sigma))

    mu_star, phi_star, sigma_star = best_key
    return best_solutions, mu_star, phi_star, sigma_star


def main():
    parser = argparse.ArgumentParser(description="PMH-TR: Parsimonious Migration History with Tree Resolution")
    parser.add_argument("--tree", required=True, help="Clone tree file")
    parser.add_argument("--labels", required=True, help="Leaf labeling file")
    parser.add_argument("--primary", required=True, help="Primary site")
    parser.add_argument("--pattern-set", default="PS,S,M,R", help="Pattern set")
    parser.add_argument("-o", "--output", default="patient1_tr/", help="Output directory")

    args = parser.parse_args()

    tree = read_clone_tree(args.tree)
    leaf_labeling = read_leaf_labeling(args.labels)
    sites = SiteIndex(set(leaf_labeling.values()), primary_site_label=args.primary)

    print(f"Loaded tree with {len(tree.vertices)} vertices")
    if tree.has_polytomies():
        print(f"Polytomies at: {[tree.vertices[i] for i in tree.get_polytomy_nodes()]}")

    pattern_sets = (
        ["PS", "PS,S", "PS,S,M", "PS,S,M,R"]
        if args.pattern_set.replace(" ", "") == "PS,S,M,R"
        else [args.pattern_set.replace(" ", "")]
    )

    results = []
    for ps in pattern_sets:
        best_solutions, mu_star, phi_star, sigma_star = solve_pmh_tr_with_pattern_set(
            tree, sites, leaf_labeling, ps
        )
        results.append((ps, best_solutions, mu_star, phi_star, sigma_star))

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # result.txt
    result_file = output_dir / "result.txt"
    with open(result_file, "w") as f:
        for ps, best_solutions, mu_star, phi_star, sigma_star in results:
            ps_norm = ps.replace(" ", "")
            pattern_type = "pR" if ps_norm == "PS,S,M,R" else "pPS"
            lb = ub = mu_star * 1000 + phi_star
            f.write(
                f"{args.primary}-\t({ps_norm})\t{mu_star}\t{phi_star}\t{sigma_star}\t"
                f"{pattern_type}\t{lb}\t{ub}\t0.0\n"
            )

    # Detailed results
    for ps, best_solutions, mu_star, phi_star, sigma_star in results:
        ps_norm = ps.replace(" ", "")
        pattern_dir = Path("results") / f"{args.primary}_{ps_norm}"
        pattern_dir.mkdir(parents=True, exist_ok=True)

        # Summary
        with open(pattern_dir / "summary.txt", "w") as f:
            f.write(f"pattern_set\t{ps_norm}\n")
            f.write(f"primary_site\t{args.primary}\n")
            f.write(f"num_sites\t{sites.m}\n")
            f.write(f"num_vertices\t{len(tree.vertices)}\n")
            f.write(f"mu_opt\t{mu_star}\n")
            f.write(f"phi_opt\t{phi_star}\n")
            f.write(f"sigma_opt\t{sigma_star}\n")

        # Labelings
        with open(pattern_dir / "labelings_opt.txt", "w") as f:
            for idx, (labeling, mu, phi, sigma) in enumerate(best_solutions):
                f.write(f"# solution {idx}\tmu={mu}\tphi={phi}\tsigma={sigma}\n")
                for i in sorted(labeling.keys(), key=lambda x: tree.vertices[x]):
                    f.write(f"{tree.vertices[i]}\t{sites.labels[labeling[i]]}\n")
                f.write("\n")

        # Site graphs
        sg_dir = pattern_dir / "site_graphs_opt"
        sg_dir.mkdir(exist_ok=True)
        for idx, (labeling, mu, phi, sigma) in enumerate(best_solutions):
            edges = site_graph_from_labeling(tree, labeling)
            write_site_graph(sg_dir / f"site_graph_{idx}.txt", sites, edges)

    # Terminal output
    print("\nResults:")
    for ps, best_solutions, mu_star, phi_star, sigma_star in results:
        ps_norm = ps.replace(" ", "")
        print(f"{args.primary}-\t({ps_norm})\t{mu_star}\t{phi_star}\t{sigma_star}")
    print(f"\nResults written to: {output_dir}")


if __name__ == "__main__":
    main()
