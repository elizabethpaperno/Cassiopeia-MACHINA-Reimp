#!/usr/bin/env python3
"""
draw_site_graph_gv.py

Render a site graph using Graphviz (dot) in a layered / hierarchical style,
with gradient arrows half-colored by source site and half by target site.
"""

import argparse
import os
from typing import Dict, List, Tuple
from graphviz import Digraph


# ----------------------------------------------------------------------
# I/O helpers
# ----------------------------------------------------------------------

def read_site_graph(path: str) -> Tuple[List[str], List[Tuple[str, str, int]]]:
    edges: List[Tuple[str, str, int]] = []
    sites = set()

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) not in (2, 3):
                raise ValueError(f"Invalid site-graph line: {line!r}")
            src, dst = parts[0], parts[1]
            w = int(parts[2]) if len(parts) == 3 else 1
            edges.append((src, dst, w))
            sites.add(src)
            sites.add(dst)

    return sorted(sites), edges


def read_colormap(path: str) -> Dict[str, str]:
    cmap = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            site, color = parts[0], parts[1]
            cmap[site] = color
    return cmap


# ----------------------------------------------------------------------
# Drawing with Graphviz
# ----------------------------------------------------------------------

def draw_site_graph_gv(site_graph_path: str,
                       colormap_path: str,
                       output_path: str,
                       primary: str | None = None) -> None:
    # Read inputs
    sites, edges = read_site_graph(site_graph_path)
    cmap = read_colormap(colormap_path)

    # Determine output format from extension
    base, ext = os.path.splitext(output_path)
    fmt = ext.lstrip(".") or "png"

    dot = Digraph(format=fmt)

    # Global graph attributes for a layered, straight-edge layout
    dot.graph_attr.update(
        rankdir="TB",      # top to bottom
        splines="spline",  
        nodesep="0.45",
        ranksep="0.75",
    )

    if primary:
        dot.graph_attr["root"] = primary

    # Nodes: rectangular boxes with colored borders
    for s in sites:
        border_color = cmap.get(s, "black")
        dot.node(
            s,
            label=s,
            shape="box",
            style="filled",
            fillcolor="white",
            color=border_color,
            penwidth="2",
        )

    # Edges: one arrow per unit weight, gradient-colored from src to dst
    for src, dst, w in edges:
        color_src = cmap.get(src, "black")
        color_dst = cmap.get(dst, "black")

        for _ in range(w):
            dot.edge(
                src,
                dst,
                color=f"{color_src};0.5:{color_dst}",  # gradient tail→head
                gradientangle="0",
                arrowsize="0.9",
                penwidth="2",
                minlen="1",        # mild spacing for layering
                # no label, no extra decorations
            )

    # Render
    dot.render(base, cleanup=True)


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Draw layered site graph with gradient arrows using Graphviz."
    )
    ap.add_argument("--site-graph", required=True,
                    help="Site graph text file (SRC DST WEIGHT)")
    ap.add_argument("--colors", required=True,
                    help="Colormap file (SITE COLOR)")
    ap.add_argument("--output", required=True,
                    help="Output image path (e.g. output.png)")
    ap.add_argument("--primary", default=None,
                    help="Primary site (placed at the top)")

    args = ap.parse_args()

    draw_site_graph_gv(
        site_graph_path=args.site_graph,
        colormap_path=args.colors,
        output_path=args.output,
        primary=args.primary,
    )


if __name__ == "__main__":
    main()
