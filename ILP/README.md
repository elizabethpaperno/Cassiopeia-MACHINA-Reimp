# Cassiopeia ILP

Folder for the ILP Cassiopeia implementation


- **cassiopeia_ilp.py**  
  Main entry point, parses input character matrix, builds potential graph, runs ILP solver, outputs lineage tree

- **potential_graph.py**  
  createst the potential graph over observed + inferred ancestral states

- **steiner_ilp.py**  
  solves Steiner tree that ILP uses - selects min-cost connecting root to all observed states

- **tree_utils.py**  
  functions for ILP process - funcs for parsing mutation states, handling missing data, computing distances, and for creating Newick tree


