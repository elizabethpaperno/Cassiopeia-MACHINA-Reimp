"""
pmh-ti.py

Parsimonious Migration History - Tree Inference (PMH-TI)

Builds on the pmh.py implementation (Sankoff DP and ILP modes) to jointly infer:
  - frequency matrix F (descendant counts per site per clone node),
  - a resolved clone tree T (when the provided tree is ambiguous),
  - a vertex labeling l of T,
that minimizes migration number mu and then comigration number phi.

Strategy:
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

import pmh as pmh_og

MAX_POLYTOMY = 5  

def find_polytomies(tree):
    """Return list of vertex indices that are polytomies (degree > 2)."""
    return [i for i in range(len(tree.vertices)) if len(tree.children[i]) > 2]

def make_edges_from_tree(tree):
    """Return list of (parent_label, child_label) pairs from a pmh.CloneTree."""
    return [(tree.vertices[u], tree.vertices[v]) for u in range(len(tree.vertices)) for v in tree.children[u]]

def build_clonetree_from_edges(edges):
    """Construct pmh.CloneTree from an edges list of labels."""
    return pmh_og.CloneTree(edges)

def compute_site_support(leaf_labeling):
    """
    Count how many leaves map to each anatomical site.

    Returns:
        dict: site_label -> count
    """
    support = {}
    for site in leaf_labeling.values():
        support[site] = support.get(site, 0) + 1
    return support


def compute_support_penalty(site_graph_edges, site_support):
    """
    Penalize biologically implausible migrations where a site
    with lower support seeds a site with higher support.

    Params:
        site_graph_edges: list of (src_site, dst_site, weight)
        site_support: dict mapping site -> number of leaves

    Returns:
        int penalty score (lower is better)
    """
    penalty = 0
    for src, dst, w in site_graph_edges:
        src_support = site_support.get(src, 0)
        dst_support = site_support.get(dst, 0)

        if src_support < dst_support:
            penalty += w

    return penalty

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
        return [[(parent_label, ch) for ch in child_labels]]

    if k > MAX_POLYTOMY:
        raise ValueError(f"Polytomy size {k} > MAX_POLYTOMY ({MAX_POLYTOMY}); refusing to refine automatically.")

    refinements = []
    counter = {"c": 0}

    def recurse(active, edges_acc):
        if len(active) == 1:
            final_edges = edges_acc + [(parent_label, active[0])]
            refinements.append(final_edges)
            return
        n = len(active)
        for i in range(n):
            for j in range(i+1, n):
                a = active[i]
                b = active[j]
                new_label = f"{internal_label_prefix}{counter['c']}"
                counter['c'] += 1
                new_edges = edges_acc + [(new_label, a), (new_label, b)]
                new_active = [x for idx, x in enumerate(active) if idx not in (i, j)]
                new_active.append(new_label)
                recurse(new_active, new_edges)

    recurse(list(child_labels), [])
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
    descendant_counts = {v: {s: 0 for s in sites.labels} for v in tree.vertices}

    children = tree.children
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
        for leaf_idx in desc_leaves:
            s_idx = labeling[leaf_idx]
            site_label = sites.labels[s_idx]
            descendant_counts[vlabel][site_label] += 1

    return descendant_counts

def evaluate_candidate_tree(edges, leaf_labeling_map, sites, pattern_set):
    """
    Given edges (list of (parent_label, child_label)), build a CloneTree, convert
    leaf_labeling_map (label->site_label) into the required map, run solve_pmh_with_pattern_set,
    and return (tree_obj, labelings, mu_star, phi_star, sigma_star).
    """
    tree = build_clonetree_from_edges(edges)

    labelings, mu_star, phi_star, sigma_star = pmh_og.solve_pmh_with_pattern_set(tree, sites, leaf_labeling_map, pattern_set)
    return tree, labelings, mu_star, phi_star, sigma_star

def compute_hub_penalty(site_graph, primary_site):
    """
    Penalize non-primary sites that act as major migration hubs.
    """
    out_degree = {}
    for (s, t), w in site_graph.items():
        out_degree[s] = out_degree.get(s, 0) + w

    penalty = 0
    for site, deg in out_degree.items():
        if site != primary_site and deg >= 2:
            penalty += deg
    return penalty

def pmhti_main(tree_path, labels_path, primary_label, pattern_set, outdir):
    tree_orig = pmh_og.read_clone_tree(tree_path)
    leaf_labeling_map = pmh_og.read_leaf_labeling(labels_path)
    site_support = compute_site_support(leaf_labeling_map)
    all_sites = set(leaf_labeling_map.values())
    sites = pmh_og.SiteIndex(all_sites, primary_label)

    original_edges = make_edges_from_tree(tree_orig)

    candidates = []
    candidates.append(("original", original_edges))

    polytomies = find_polytomies(tree_orig)
    if polytomies:
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
            for idx, subedges in enumerate(subtree_refinements):
                new_edges = apply_refinement_to_edges(original_edges, parent_label, child_labels, subedges)
                candidates.append((f"refinement_{parent_label}_{idx}", new_edges))
    else:
        print("No polytomies found; using original tree only.")

    best = None  # tuple: (mu, phi, candidate_name, tree_obj, labeling, sigma, obj)
    results = []
    for name, edges in candidates:
        try:
            tree_cand, labelings, mu_star, phi_star, sigma_star = evaluate_candidate_tree(edges, leaf_labeling_map, sites, pattern_set)
        except Exception as e:
            print(f"Candidate {name} raised an error during PMH solving: {e}")
            continue
        for lab in labelings:
            """
            results.append((name, tree_cand, lab, mu_star, phi_star, sigma_star))
            key = (mu_star, phi_star)
            if best is None or key < (best[0], best[1]):
                best = (mu_star, phi_star, name, tree_cand, lab, sigma_star)
            """
            sg = pmh_og.site_graph_from_labeling(tree_cand, lab)
            hub_penalty = compute_hub_penalty(sg, primary_label)

            results.append(
                (name, tree_cand, lab, mu_star, phi_star, sigma_star, hub_penalty)
            )

            key = (mu_star, phi_star, sigma_star, hub_penalty)

            if best is None or key < (best[0], best[1], best[5], best[6]):
                best = (
                    mu_star,
                    phi_star,
                    name,
                    tree_cand,
                    lab,
                    sigma_star,
                    hub_penalty,
                )

    if best is None:
        raise RuntimeError("No feasible solution found for any candidate tree.")

    #mu_opt, phi_opt, cand_name, best_tree, best_lab, sigma_opt = best
    mu_opt, phi_opt, cand_name, best_tree, best_lab, sigma_opt, hub_penalty_opt = best
    print(f"Best candidate: {cand_name} with mu={mu_opt}, phi={phi_opt}, sigma={sigma_opt}")

    F = build_frequency_matrix(best_tree, best_lab, sites)

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

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

    lab_path = outdir / "labeling_pmh_ti.txt"
    with open(lab_path, "w") as f:
        f.write(f"# candidate {cand_name}, mu={mu_opt}, phi={phi_opt}, sigma={sigma_opt}\n")
        for i in sorted(best_lab.keys(), key=lambda idx: best_tree.vertices[idx]):
            v = best_tree.vertices[i]
            site_index = best_lab[i]
            f.write(f"{v}\t{sites.labels[site_index]}\n")

    F_path = outdir / "frequency_matrix_F.txt"
    with open(F_path, "w") as f:
        header = ["vertex"] + sites.labels
        f.write("\t".join(header) + "\n")
        for vlabel, site_counts in F.items():
            row = [vlabel] + [str(site_counts[s]) for s in sites.labels]
            f.write("\t".join(row) + "\n")

    sg = pmh_og.site_graph_from_labeling(best_tree, best_lab)
    sg_path = outdir / "site_graph_pmh_ti.txt"
    pmh_og.write_site_graph(sg_path, sites, sg)

    print(f"Wrote PMH-TI outputs to {outdir}")
    return outdir

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

if __name__ == "__main__":
    main()
