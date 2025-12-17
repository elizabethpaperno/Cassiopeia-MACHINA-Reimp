# Cassiopeia ILP

Folder for the Hybrid Cassiopeia implementation

### Newick Tree
The command takes in an input tab-delimited file with index=cell and column state. Use `-t` to indicate the ILP threshold, which is the maximum group size to switch from greedy to ILP. Output will be saved at `-o`. Use `--debug` to get print output for greedy splits and ILP group sizes.
```
python -m Cassiopeia.Hybrid.cassiopeia_hybrid cell_state.txt -t 10 --debug -o output.nwk
```

### Calculating parsimony cost
```
python -m Cassiopeia.compute_parsimony tree.nwk
```