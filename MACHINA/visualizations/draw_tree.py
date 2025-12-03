"""Draw a site graph using Graphviz."""

import argparse
import os
from graphviz import Digraph


def read_site_graph(path):
    """Read site graph text file into site list and weighted edges"""
    edges = []
    sites = set()

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) not in (2, 3):
                raise ValueError("Invalid site-graph line: {!r}".format(line))
            src, dst = parts[0], parts[1]
            w = int(parts[2]) if len(parts) == 3 else 1
            edges.append((src, dst, w))
            sites.add(src)
            sites.add(dst)

    return sorted(sites), edges


def read_colormap(path):
    """Read site to color mapping from text file"""
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


def draw_site_graph_gv(site_graph_path, colormap_path, output_path, primary=None):
    """Draw the site graph as a layered Graphviz figure"""
    sites, edges = read_site_graph(site_graph_path)
    cmap = read_colormap(colormap_path)

    base, ext = os.path.splitext(output_path)
    fmt = ext.lstrip(".") or "png"

    dot = Digraph(format=fmt)

    dot.graph_attr.update(
        rankdir="TB",
        splines="spline",
        nodesep="0.45",
        ranksep="0.75",
    )

    if primary:
        dot.graph_attr["root"] = primary

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

    for src, dst, w in edges:
        color_src = cmap.get(src, "black")
        color_dst = cmap.get(dst, "black")

        for _ in range(w):
            dot.edge(
                src,
                dst,
                color="{};0.5:{}".format(color_src, color_dst),
                gradientangle="0",
                arrowsize="0.9",
                penwidth="2",
                minlen="1",
            )

    dot.render(base, cleanup=True)


def main():
    """Parse command line arguments and draw the site graph."""
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
