# Cassiopeia Greedy

Folder for the Cassiopeia-Greedy implementation

### Input and Running the Algorithm

The command takes in an input tab-delimited file with index=cell and column state (we'll refer to this input input.txt). Running the following command executes the greedy algorithm on input.txt and prints the result to the terminal as well as storing it in tree.newick.

```
python cassiopeia_greedy.py input.txt

```

### Visualizing Output Tree

```
python reconstruct.py tree.newick
```

This visualizes the tree that is outputted by the Cassiopeia-Greedy algorithm with the input of input.txt.
