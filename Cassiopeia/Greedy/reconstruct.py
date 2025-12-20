import argparse
from ete3 import Tree, TreeStyle, NodeStyle, faces, AttrFace

def layout(node):
    if node.is_leaf():
        name_face = AttrFace("name", fsize=10)
        faces.add_face_to_node(name_face, node, column=0)

def main():
    parser = argparse.ArgumentParser(description="Visualize Newick tree")
    parser.add_argument("newick_file", help="Path to Newick file")
    args = parser.parse_args()

    with open(args.newick_file, "r") as f:
        newick_str = f.read().strip()

    tree = Tree(newick_str, format=1)

    nstyle = NodeStyle()
    nstyle["fgcolor"] = "black"
    nstyle["size"] = 6
    nstyle["vt_line_width"] = 1
    nstyle["hz_line_width"] = 1

    for n in tree.traverse():
        n.set_style(nstyle)
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = layout
    ts.scale = 20
    ts.title.add_face(faces.TextFace("Greedy Tree Reconstruction", fsize=10), column=0)

    tree.show(tree_style=ts)

if __name__ == "__main__":
    main()
