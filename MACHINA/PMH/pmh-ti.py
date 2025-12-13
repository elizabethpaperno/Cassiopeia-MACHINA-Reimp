"""
pmh_ti.py

PMH-TI (Parsimonious Migration History - Tree Inference) wrapper that builds on
the pmh.py implementation (Sankoff DP and ILP modes) to jointly infer:

  - frequency matrix F (descendant counts per site per clone node),
  - a resolved clone tree T (when the provided tree is ambiguous),
  - a vertex labeling l of T,

that minimizes migration number mu and then comigration number phi.

Strategy (practical, limited):
- Read unresolved tree + leaf labels using pmh.py helpers.
- Find polytomies (nodes with >2 children). For polytomies with degree <= MAX_POLYTOMY
  generate binary refinements (all valid full binary groupings of children).
- For each candidate resolved tree:
    - call pmh.solve_pmh_with_pattern_set(...) using the provided pattern-set.
    - record mu, phi and choose the candidate with minimal (mu, phi).
- Output: chosen tree, labeling(s), mu, phi, sigma, and frequency matrix F.

"""

import argparse
from pathlib import Path
from itertools import combinations
import itertools
import copy
import math

import pmh as pmh_og

MAX_POLYTOMY = 5  # maximum polytomy degree to attempt full refinement (avoid combinatorial explosion)

def find_polytomies(tree):
    """Return list of vertex indices that are polytomies (degree > 2)."""
    return [i for i in range(len(tree.vertices)) if len(tree.children[i]) > 2]

def make_edges_from_tree(tree):
    """Return list of (parent_label, child_label) pairs from a pmh.CloneTree."""
    return [(tree.vertices[u], tree.vertices[v]) for u in range(len(tree.vertices)) for v in tree.children[u]]

def build_clonetree_from_edges(edges):
    """Construct pmh.CloneTree from an edges list of labels."""
    return pmh_og.CloneTree(edges)

def generate_binary_refinements_for_node(parent_label, child_labels, internal_label_prefix="X"):
    """
    Given:
      - parent_label: label of the polytomy parent (string)
      - child_labels: list of child labels (strings) to be binary-refined
    Return:
      - list of edges_subsets (list of (parent_child) edges) that replace parent->child_labels
        with a small binary subtree connecting the children under parent_label.
    Note:
      - This function generates all ways to pair children into a binary tree by iteratively
        combining two elements. The number grows as (2k-3)!! for k children, so we limit k.
    """
    k = len(child_labels)
    if k <= 2:
        # nothing to refine
        return [[(parent_label, ch) for ch in child_labels]]

    if k > MAX_POLYTOMY:
        raise ValueError(f"Polytomy size {k} > MAX_POLYTOMY ({MAX_POLYTOMY}); refusing to refine automatically.")

    # We'll generate binary trees by repeatedly pairing two items into a new internal node.
    # Represent each active element as its label (strings). When pairing a and b,
    # create a new internal label like "X_0", etc., store edges (new -> a, new -> b),
    # and replace a,b in active list with the new label, continuing until one remains.
    refinements = []
    counter = {"c": 0}

    def recurse(active, edges_acc):
        # base
        if len(active) == 1:
            # Connect parent_label to the last active element
            final_edges = edges_acc + [(parent_label, active[0])]
            refinements.append(final_edges)
            return
        # choose unordered pairs (i<j)
        n = len(active)
        for i in range(n):
            for j in range(i+1, n):
                a = active[i]
                b = active[j]
                new_label = f"{internal_label_prefix}{counter['c']}"
                counter['c'] += 1
                # edges for new internal node
                new_edges = edges_acc + [(new_label, a), (new_label, b)]
                # new active list
                new_active = [x for idx, x in enumerate(active) if idx not in (i, j)]
                new_active.append(new_label)
                recurse(new_active, new_edges)

    recurse(list(child_labels), [])
    # refinements: each is edges inside the subtree; they implicitly create new internal node labels
    # but these labels won't be in the global vertex set yet — caller must integrate.
    # Return unique refinements (some may be isomorphic duplicates due to symmetric pairings)
    unique = []
    seen = set()
    for r in refinements:
        key = tuple(sorted(r))
        if key not in seen:
            seen.add(key)
            unique.append(r)
    return unique

def apply_refinement_to_edges(original_edges, parent_label, child_labels, subtree_edges):
    """
    Replace in original_edges the edges (parent_label, child) for child in child_labels
    with the new subtree_edges (which include connections among new internal nodes and child labels)
    and return a new edges list.
    """
    new_edges = [e for e in original_edges if not (e[0] == parent_label and e[1] in child_labels)]
    # ensure subtree internal labels are unique relative to existing labels
    # but subtree_edges already use unique internal labels like X_0...
    new_edges.extend(subtree_edges)
    return new_edges


# Frequency matrix builder:

def build_frequency_matrix(tree, labeling, sites):
    """
    Build a simple frequency matrix F: for each vertex v and each site s,
    count the number of descendant leaves of v that are labeled as site s
    (based on leaf_labeling). The result is a dict mapping vertex_label -> {site_label: count}.
    """
    n = len(tree.vertices)
    # find which leaves correspond to which site according to labeling (labeling uses node indices)
    # But labeling holds assigned site indices for all nodes; we want descendant leaf sites from original leaf labels.
    # So instead, we will compute for each vertex its descendant leaves' original leaf site labels.
    # First, find all leaf labels mapped to site from tree.leaf nodes (we assume pmh.read_leaf_labeling used)
    # We'll traverse tree and collect counts.
    descendant_counts = {v: {s: 0 for s in sites.labels} for v in tree.vertices}

    # Map leaf index -> site label (if leaf labeling is present in the original pmh context, we used that earlier).
    # However here we only have labeling (assignment of site indices to all vertices), not original leaf labels.
    # To compute descendant leaf counts we must use the original leaf labeling; the caller should pass it.
    # For usability, we'll instead compute F as counts of descendant nodes assigned to each site:
    # For node v, for each descendant u (including v), increment counts[site_of_u] by 1 if u is a leaf.
    # So this F approximates clone frequencies by counting leaf descendants per site.
    # We'll compute descendant leaf indices per vertex.
    children = tree.children
    # compute descendants via dfs
    def get_descendant_leaves(u):
        out = []
        stack = [u]
        while stack:
            curr = stack.pop()
            if not children[curr]:
                out.append(curr)
            else:
                for ch in children[curr]:
                    stack.append(ch)
        return out

    for vi, vlabel in enumerate(tree.vertices):
        desc_leaves = get_descendant_leaves(vi)
        # for each descendant leaf, find its assigned site via labeling[leaf_index]
        for leaf_idx in desc_leaves:
            s_idx = labeling[leaf_idx]
            site_label = sites.labels[s_idx]
            descendant_counts[vlabel][site_label] += 1

    return descendant_counts


# Main PMH-TI logic

def evaluate_candidate_tree(edges, leaf_labeling_map, sites, pattern_set):
    """
    Given edges (list of (parent_label, child_label)), build a CloneTree, convert
    leaf_labeling_map (label->site_label) into the required map, run solve_pmh_with_pattern_set,
    and return (tree_obj, labelings, mu_star, phi_star, sigma_star).
    """
    # Build tree object
    tree = build_clonetree_from_edges(edges)

    # run pmh.solve_pmh_with_pattern_set
    labelings, mu_star, phi_star, sigma_star = pmh_og.solve_pmh_with_pattern_set(tree, sites, leaf_labeling_map, pattern_set)
    # labelings returned by pmh_og are lists of dicts mapping node indices -> site_index (these indices correspond to tree.vertices ordering)
    return tree, labelings, mu_star, phi_star, sigma_star

def pmhti_main(tree_path, labels_path, primary_label, pattern_set, outdir):
    # Read inputs using pmh helpers
    tree_orig = pmh_og.read_clone_tree(tree_path)
    leaf_labeling_map = pmh_og.read_leaf_labeling(labels_path)
    all_sites = set(leaf_labeling_map.values())
    sites = pmh_og.SiteIndex(all_sites, primary_label)

    # Original edges (labels)
    original_edges = make_edges_from_tree(tree_orig)

    # Candidate trees list (edges lists)
    candidates = []
    # include original first
    candidates.append(("original", original_edges))

    # find polytomies
    polytomies = find_polytomies(tree_orig)
    if polytomies:
        # attempt to generate single-polytomy refinements (for complexity reasons)
        # We'll refine only one polytomy at a time (the largest), generating multiple candidate trees.
        # Pick the polytomy with largest degree
        degrees = [(i, len(tree_orig.children[i])) for i in polytomies]
        degrees.sort(key=lambda x: -x[1])
        node_to_refine = degrees[0][0]
        parent_label = tree_orig.vertices[node_to_refine]
        child_indices = tree_orig.children[node_to_refine]
        child_labels = [tree_orig.vertices[i] for i in child_indices]
        if len(child_labels) > MAX_POLYTOMY:
            print(f"Polytomy of size {len(child_labels)} > MAX_POLYTOMY ({MAX_POLYTOMY}). Skipping automatic refinements.")
        else:
            print(f"Generating binary refinements for polytomy at node {parent_label} with {len(child_labels)} children")
            subtree_refinements = generate_binary_refinements_for_node(parent_label, child_labels, internal_label_prefix="X_")
            # For each refinement, integrate into original edge list to produce a new candidate tree
            for idx, subedges in enumerate(subtree_refinements):
                new_edges = apply_refinement_to_edges(original_edges, parent_label, child_labels, subedges)
                # Additionally, ensure we don't accidentally duplicate labels: subedges may introduce internal labels "X_*"
                candidates.append((f"refinement_{parent_label}_{idx}", new_edges))
    else:
        print("No polytomies found; using original tree only.")

    # Evaluate each candidate
    best = None  # tuple: (mu, phi, candidate_name, tree_obj, labeling, sigma, obj)
    results = []
    for name, edges in candidates:
        try:
            tree_cand, labelings, mu_star, phi_star, sigma_star = evaluate_candidate_tree(edges, leaf_labeling_map, sites, pattern_set)
        except Exception as e:
            print(f"Candidate {name} raised an error during PMH solving: {e}")
            continue
        # labelings is a list of labelings (dictionaries)
        # Choose first labeling as canonical among returned optimums
        for lab in labelings:
            results.append((name, tree_cand, lab, mu_star, phi_star, sigma_star))
            # update best
            key = (mu_star, phi_star)
            if best is None or key < (best[0], best[1]):
                best = (mu_star, phi_star, name, tree_cand, lab, sigma_star)
    if best is None:
        raise RuntimeError("No feasible solution found for any candidate tree.")

    mu_opt, phi_opt, cand_name, best_tree, best_lab, sigma_opt = best
    print(f"Best candidate: {cand_name} with mu={mu_opt}, phi={phi_opt}, sigma={sigma_opt}")

    # build frequency matrix F using best_tree and best_lab
    F = build_frequency_matrix(best_tree, best_lab, sites)

    # Prepare output dir
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Write summary
    summary_path = outdir / "summary_pmh_ti.txt"
    with open(summary_path, "w") as f:
        f.write(f"chosen_candidate\t{cand_name}\n")
        f.write(f"pattern_set\t{pattern_set}\n")
        f.write(f"primary\t{primary_label}\n")
        f.write(f"num_sites\t{sites.m}\n")
        f.write(f"num_vertices\t{len(best_tree.vertices)}\n")
        f.write(f"mu_opt\t{mu_opt}\n")
        f.write(f"phi_opt\t{phi_opt}\n")
        f.write(f"sigma_opt\t{sigma_opt}\n")

    # Write chosen labeling
    lab_path = outdir / "labeling_pmh_ti.txt"
    with open(lab_path, "w") as f:
        f.write(f"# candidate {cand_name}, mu={mu_opt}, phi={phi_opt}, sigma={sigma_opt}\n")
        for i in sorted(best_lab.keys(), key=lambda idx: best_tree.vertices[idx]):
            v = best_tree.vertices[i]
            site_index = best_lab[i]
            f.write(f"{v}\t{sites.labels[site_index]}\n")

    # Write frequency matrix (tabular)
    F_path = outdir / "frequency_matrix_F.txt"
    with open(F_path, "w") as f:
        header = ["vertex"] + sites.labels
        f.write("\t".join(header) + "\n")
        for vlabel, site_counts in F.items():
            row = [vlabel] + [str(site_counts[s]) for s in sites.labels]
            f.write("\t".join(row) + "\n")

    # Write site graph
    sg = pmh_og.site_graph_from_labeling(best_tree, best_lab)
    sg_path = outdir / "site_graph_pmh_ti.txt"
    pmh_og.write_site_graph(sg_path, sites, sg)

    print(f"Wrote PMH-TI outputs to {outdir}")
    return outdir

# -------------------------
# CLI
# -------------------------

"""
def main():
    parser = argparse.ArgumentParser(description="PMH-TI wrapper using pmh.py helpers.")
    parser.add_argument("--tree", required=True, help="Clone tree edge list file (parent child per line)")
    parser.add_argument("--labels", required=True, help="Leaf labeling file (leaf site per line)")
    parser.add_argument("--primary", required=True, help="Primary anatomical site label")
    parser.add_argument("--pattern-set", required=True, help="Pattern set (e.g., 'PS', 'PS,S', 'PS,S,M', 'PS,S,M,R')")
    parser.add_argument("--outdir", default="results_pmh_ti", help="Output directory")
    args = parser.parse_args()

    outdir = pmhti_main(args.tree, args.labels, args.primary, args.pattern_set, args.outdir)
    print("Done.")
"""
def main():
    parser = argparse.ArgumentParser(
        description=(
            "PMH-TI: Parsimonious Migration History with Tree Inference. "
            "Jointly infers clone tree structure and migration history "
            "under pattern sets {PS}, {PS,S}, {PS,S,M}, {PS,S,M,R}."
        )
    )
    parser.add_argument("--tree", required=True, help="Unresolved clone tree file")
    parser.add_argument("--labels", required=True, help="Leaf tissue labeling file")
    parser.add_argument("--primary", required=True, help="Primary anatomical site")
    parser.add_argument(
        "--pattern-set",
        required=True,
        help="One of: 'PS', 'PS,S', 'PS,S,M', 'PS,S,M,R'",
    )

    args = parser.parse_args()

    # ----------------------------
    # Load inputs
    # ----------------------------
    tree = read_clone_tree(args.tree)
    leaf_labeling = read_leaf_labeling(args.labels)
    sites = SiteIndex(set(leaf_labeling.values()), primary_site_label=args.primary)

    # ----------------------------
    # Solve PMH-TI
    # Returns: list of (tree, labeling) pairs
    # ----------------------------
    solutions, mu_star, phi_star, sigma_star = solve_pmh_ti_with_pattern_set(
        tree, sites, leaf_labeling, args.pattern_set
    )

    normalized = args.pattern_set.replace(" ", "")
    num_mp = len(solutions)

    # ----------------------------
    # Compute stats for each solution
    # ----------------------------
    solution_stats = []
    for idx, (T, labeling) in enumerate(solutions):
        mu, phi, sigma = compute_migration_stats(T, sites, labeling)
        solution_stats.append((idx, T, labeling, mu, phi, sigma))

    # ----------------------------
    # Select optimal solutions
    # (PMH-TI is already ILP-optimized, but keep logic consistent)
    # ----------------------------
    opt_solutions = [
        (idx, T, lab, mu, phi, sigma)
        for (idx, T, lab, mu, phi, sigma) in solution_stats
        if (mu, phi, sigma) == (mu_star, phi_star, sigma_star)
    ]
    num_opt = len(opt_solutions)

    # ----------------------------
    # Output directory
    # ----------------------------
    base_dir = Path("results")
    pattern_dir = base_dir / f"{args.primary}_{normalized}"
    pattern_dir.mkdir(parents=True, exist_ok=True)

    # ----------------------------
    # Summary file
    # ----------------------------
    with open(pattern_dir / "summary.txt", "w") as f:
        f.write(f"pattern_set\t{normalized}\n")
        f.write(f"primary_site\t{args.primary}\n")
        f.write(f"num_sites\t{sites.m}\n")
        f.write(f"num_vertices\t{len(tree.vertices)}\n")
        f.write(f"mu_opt\t{mu_star}\n")
        f.write(f"phi_opt\t{phi_star}\n")
        f.write(f"sigma_opt\t{sigma_star}\n")
        f.write(f"n_mp_solutions\t{num_mp}\n")
        f.write(f"n_opt_solutions\t{num_opt}\n")

    # ----------------------------
    # All MP labelings
    # ----------------------------
    with open(pattern_dir / "labelings_all_mp.txt", "w") as f:
        for idx, T, lab, mu, phi, sigma in solution_stats:
            f.write(f"# solution {idx}\tmu={mu}\tphi={phi}\tsigma={sigma}\n")
            for i in sorted(lab.keys(), key=lambda x: T.vertices[x]):
                f.write(f"{T.vertices[i]}\t{sites.labels[lab[i]]}\n")
            f.write("\n")

    # ----------------------------
    # Optimal labelings
    # ----------------------------
    with open(pattern_dir / "labelings_opt.txt", "w") as f:
        for idx, T, lab, mu, phi, sigma in opt_solutions:
            f.write(f"# solution {idx}\tmu={mu}\tphi={phi}\tsigma={sigma}\n")
            for i in sorted(lab.keys(), key=lambda x: T.vertices[x]):
                f.write(f"{T.vertices[i]}\t{sites.labels[lab[i]]}\n")
            f.write("\n")

    # ----------------------------
    # Site graphs: all MP
    # ----------------------------
    sg_all_dir = pattern_dir / "site_graphs_all_mp"
    sg_all_dir.mkdir(exist_ok=True)
    for idx, T, lab, mu, phi, sigma in solution_stats:
        edges = site_graph_from_labeling(T, lab)
        write_site_graph(sg_all_dir / f"site_graph_{idx}.txt", sites, edges)

    # ----------------------------
    # Site graphs: optimal
    # ----------------------------
    sg_opt_dir = pattern_dir / "site_graphs_opt"
    sg_opt_dir.mkdir(exist_ok=True)
    for idx, T, lab, mu, phi, sigma in opt_solutions:
        edges = site_graph_from_labeling(T, lab)
        write_site_graph(sg_opt_dir / f"site_graph_{idx}.txt", sites, edges)

    # ----------------------------
    # Terminal output
    # ----------------------------
    print(f"P = {normalized}, primary = {args.primary}")
    print(f"mu* = {mu_star}, phi* = {phi_star}, sigma* = {sigma_star}")
    print(
        f"#MP solutions = {num_mp}, #optimal solutions = {num_opt}"
    )
    print(f"Results written to: {pattern_dir}")


if __name__ == "__main__":
    main()
