"""Draw PMH-TR site graphs using Graphviz (based on PMH visualization script)."""

import argparse
import os
import sys
from pathlib import Path

# Set Graphviz executable path for Windows
if sys.platform == "win32":
    os.environ["PATH"] += os.pathsep + r"C:\Program Files\Graphviz\bin"

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


def generate_colormap(sites):
    """Generate a default colormap for sites using a color palette"""
    colors = [
        "#FF6B6B",  # Red
        "#4ECDC4",  # Teal
        "#45B7D1",  # Blue
        "#FFA07A",  # Light Salmon
        "#98D8C8",  # Mint
        "#F7DC6F",  # Yellow
        "#BB8FCE",  # Purple
        "#85C1E2",  # Light Blue
        "#F8B88B",  # Peach
        "#ABEBC6",  # Green
    ]
    
    colormap = {}
    for i, site in enumerate(sites):
        colormap[site] = colors[i % len(colors)]
    return colormap


def draw_site_graph_gv(site_graph_path, output_path, colormap=None, primary=None):
    """Draw the site graph as a layered Graphviz figure"""
    sites, edges = read_site_graph(site_graph_path)
    
    if not colormap:
        colormap = generate_colormap(sites)

    base, ext = os.path.splitext(output_path)
    fmt = ext.lstrip(".") or "png"

    dot = Digraph(format=fmt)

    dot.graph_attr.update(
        rankdir="TB",
        splines="spline",
        nodesep="0.45",
        ranksep="0.75",
        bgcolor="white",
    )

    if primary:
        dot.graph_attr["root"] = primary

    for s in sites:
        fill_color = colormap.get(s, "#E8E8E8")
        dot.node(
            s,
            label=s,
            shape="box",
            style="filled",
            fillcolor=fill_color,
            color="black",
            penwidth="2.5",
            fontsize="12",
            fontname="helvetica",
        )

    for src, dst, w in edges:
        color_src = colormap.get(src, "#E8E8E8")
        color_dst = colormap.get(dst, "#E8E8E8")

        for _ in range(w):
            dot.edge(
                src,
                dst,
                color="{};0.5:{}".format(color_src, color_dst),
                gradientangle="0",
                arrowsize="1.0",
                penwidth="2.5",
                minlen="1",
                fontsize="11",
                fontcolor="darkred",
            )
        
        # Add weight label if count > 1
        if w > 1:
            dot.edge(
                src,
                dst,
                label=str(w),
                fontsize="11",
                fontcolor="darkred",
                style="invis",
            )

    dot.render(base, cleanup=True)
    print(f"  ✓ Generated: {output_path}")


def process_results_directory(results_dir):
    """Process all site graphs in results directory structure"""
    results_path = Path(results_dir)
    
    if not results_path.exists():
        print(f"Error: Results directory not found: {results_dir}")
        sys.exit(1)
    
    # Find all site_graph*.txt files
    site_graph_files = list(results_path.rglob("site_graph_*.txt"))
    
    if not site_graph_files:
        print(f"No site graphs found in {results_dir}")
        return
    
    print(f"Found {len(site_graph_files)} site graph(s) to process...\n")
    
    for txt_file in sorted(site_graph_files):
        # Create output PNG in same directory as txt file
        png_file = txt_file.with_suffix(".png")
        
        # Get pattern set from parent directory
        pattern_dir = txt_file.parent.parent.name
        
        print(f"Processing {pattern_dir}/{txt_file.stem}...")
        
        try:
            draw_site_graph_gv(str(txt_file), str(png_file))
        except Exception as e:
            print(f"  ✗ Error: {e}")
    
    print("\n✅ All site graphs processed!")


def main():
    """Parse command line arguments and draw site graphs."""
    ap = argparse.ArgumentParser(
        description="Draw PMH-TR site graphs with colored nodes and gradient arrows using Graphviz."
    )
    ap.add_argument(
        "--results-dir",
        type=str,
        default=None,
        help="Results directory to process all site graphs (e.g., pmh_tr_analysis/results)"
    )
    ap.add_argument(
        "--site-graph",
        type=str,
        default=None,
        help="Individual site graph text file (SRC DST WEIGHT)"
    )
    ap.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output image path (e.g. output.png)"
    )

    args = ap.parse_args()

    if args.results_dir:
        # Process entire results directory
        process_results_directory(args.results_dir)
    elif args.site_graph and args.output:
        # Process single file
        draw_site_graph_gv(args.site_graph, args.output)
    else:
        ap.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
