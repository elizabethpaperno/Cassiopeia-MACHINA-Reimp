# Cassiopeia MACHINA Reimplementation

## How to Clone, Install, and Run

1. Create and activate an environment

```
python3 -m venv <<name>>
cd <<name>>
. bin/activate
```

2. Clone the repository

```
git clone git@github.com:elizabethpaperno/Cassiopeia-MACHINA-Reimp.git
```

3. Navigate to the root directory

```
cd Cassiopeia-MACHINA-Reimp
```

4. Install Requirements

```
pip install -r requirements.txt
```

## Cassiopeia

This repository contains code and data for running **Greedy, ILP, Hybrid** algorithms and visualizing cell lineage trees.

### `data/`

Contains data from the Allen Institute DREAM Challenge, which features publicly available lineage tracing datasets.

The data can be found at https://www.synapse.org/Synapse:syn20692755/files/ under the "data" folder. We used input from the "Subchallenge1" folder and compared it to the corresponding output in the "groundTruth_train" folder.

Each row represents:

- cell: a unique cell identifier
- state: sequence of integers representing mutations at different sites

### `Greedy/`

Contains main script for running Cassiopeia-Greedy

### Input and Running the Algorithm

The command takes in an input tab-delimited file with index=cell and column state (we'll refer to this input as input.txt, make sure it is in the Greedy folder). Running the following command executes the greedy algorithm on input.txt and prints the result to the terminal as well as storing it in tree.newick.

```
python cassiopeia_greedy.py input.txt

```

### Visualizing Output Tree

```
python reconstruct.py tree.newick
```

This visualizes the tree that is outputted by the Cassiopeia-Greedy algorithm with the input of input.txt.

### `ILP/`
Contains the main script for running Cassippeia-ILP code as well as other helper scripts to reconstruct lineage states with globally optimal parsimony.
### Running Cassiopeia-ILP:
The ILP solver takes in a tab-delimited file with each row corresponding to a cell and the mutation state string for barcode sites.  
Usage (where XXX corresopnds to the dataset you want to run on).
```
python -m Cassiopeia.ILP.cassiopeia_ilp \
  Cassiopeia/data/train/sub1_train_XXX.txt \
  -o Cassiopeia/experiment_results/general_full_dataset/ilp/ilp_sub1_train_XXX.nwk
```
Optional flags:
- (--debug): will write additional diagnostics to a debug.txt file for more information on the script run
- (--time): prints the total runtime.

The ILP script outputs and stores a lineage tree in .nwk Newick format.



### `Hybrid/`

Contains main script for running Cassiopeia-Hybrid.

### Newick Tree

The command takes in an input tab-delimited file with index=cell and column state. Use `-t` to indicate the ILP threshold, which is the maximum group size to switch from greedy to ILP. Output will be saved at `-o`. Use `--debug` to get print output for greedy splits and ILP group sizes.

```
python -m Cassiopeia.Hybrid.cassiopeia_hybrid cell_state.txt -t 10 --debug -o output.nwk
```

### Calculating parsimony cost

```
python -m Cassiopeia.compute_parsimony tree.nwk
```

### `visualizations/`

Contains script for visualizing Newick Tree.

## MACHINA

This repository contains code and data for running **PMH, PMH-TI, PMH-TR** analyses and visualizing inferred migration graphs.

### `data/`

Contains subdirectories named after the study or paper the data originates from (e.g. `mcpherson_2016`).  
Each subdirectory includes:

- `.tree` files: clone tree edge lists (`parent child` per line)
- `.labeling` files: mappings from leaf nodes to anatomical sites (`leaf site` per line)
- `.colormap` files: mappings from anatomical sites to colors (used for visualization)

---

### `visualizations/`

Contains scripts for visualizing PMH, PMH-TI, and PMH-TR outputs.

- `draw_migration_graph.py`:  
  Inputs:
  - `--site_graph`: a migration (site) graph file (e.g. produced by `pmh.py`) in the format:
    ```
    src_site  dst_site  weight
    ```
  - `--colors`: a colormap file (`site color`)
  - `--output`: file where visualization should be saved
    this script outputs a Graphviz visualization of the inferred migration (parsimony) history.

---

## PMH

### `PMH/pmh.py`

Main script for running Parsimonious Migration History.

The script:

- reads a clone tree and leaf labeling
- runs PMH under a specified pattern set
- outputs all maximum parsimony labelings and optimal solutions
- constructs a migration graphs

#### Inputs

- `--tree`: path to a `.tree` file
- `--labels`: path to a `.labeling` file
- `--primary`: primary anatomical site
- `--pattern-set`: one of `PS`, `PS,S`, `PS,S,M`, `PS,S,M,R`

---

### `PMH/results/`

Automatically created directory where all PMH outputs are stored.

For each run, a subdirectory with name {primary site}\_{pattern}\_{input file name} is created containing:

- `summary.txt`
- `labelings_all_mp.txt`
- `labelings_opt.txt`
- `site_graphs_all_mp/`
- `site_graphs_opt/`

---

## PMH Example Runs

### Run PMH

From PMH directory:

```bash
python pmh.py \
  --tree ../data/mcpherson_2016/patient1.tree \
  --labels ../data/mcpherson_2016/patient1.labeling \
  --primary ROv \
  --pattern-set "PS,S"
```

### Now Get Associated Visualization

cd to MACHINA directory:

```
cd ../
```

From MACHINA directory run:

```bash
python visualizations/draw_migration_graph.py \
  --site-graph PMH/results/ROv_PS,S_patient1/site_graphs_opt/site_graph_0.txt \
  --colors data/mcpherson_2016/patient1.colormap\
  --output PMH/results/ROv_PS,S_patient1/site_graphs_opt/site_graph_0.png \
  --primary ROv
```

---

## PMH-TI

### `PMH/pmh-ti.py`

Main script for running Parsimonious Migration History with Tree Inference.

The script:

- reads an unresolved clone tree and leaf labeling,
- identifies polytomies in the clone tree,
- generates candidate binary refinements for resolved polytomies,
- runs PMH on each candidate tree under a specified pattern set,
- selects the candidate and vertex labeling that minimize the objective,
- constructs the migration site graph and frequency matrix for the optimal solution.

#### Inputs

- `--tree`: path to a `.tree` file
- `--labels`: path to a `.labeling` file
- `--primary`: primary anatomical site
- `--pattern-set`: one of `PS`, `PS,S`, `PS,S,M`, `PS,S,M,R`
- `--outdir`: output directory

---

### `PMH/results_pmh_ti/`

Automatically created directory where all PMH-TI outputs are stored.

Each run produces the following files:

- `summary_pmh_ti.txt`
- `labelings_pmh_ti.txt`
- `site_graphs_pmh_ti.txt/`
- `frequency_matrix_F.txt/`

---

## PMH-TI Example Runs

### Run PMH-TI

From PMH-TI directory:

```bash
python pmh-ti.py \
  --tree ../data/hoadley_2016/A7.tree \
  --labels ../data/hoadley_2016/A7.labeling \
  --primary breast \
  --pattern-set PS,S,M,R
```

### Now Get Associated Visualization

cd to MACHINA directory:

```
cd ../
```

From MACHINA directory run:

```bash
python visualizations/draw_migration_graph.py \
  --site-graph PMH/results_pmh_ti/site_graph_pmh_ti.txt \
  --colors data/hoadley_2016/A7.colormap \
  --output ../figures/A7_breast_site_graph.png \
  --primary breast
```

---

## PMH-TR

### `PMH-TR/pmh_tr.py`

Main script for running Parsimonious Migration History with Tree Resolution (PMH-TR).

- Loads a clone tree and leaf site labels, fixing the primary site at the root.
- Resolves polytomies implicitly via vertex site labeling in an ILP-based PMH-TR formulation.
- Supports multiple migration pattern sets and selects solutions minimizing μ, then φ, then σ.
- Outputs optimal labelings, migration graphs, and summary results.

#### Inputs

- `--tree`: path to a `.tree` file (may contain polytomies)
- `--labels`: path to a `.labeling` file
- `--primary`: primary anatomical site
- `--pattern-set`: one of `PS`, `PS,S`, `PS,S,M`, `PS,S,M,R`
- `-o` / `--output`: output directory

#### Outputs

For each pattern set, PMH-TR creates a subdirectory containing:

- `summary.txt`: metadata and optimal statistics (mu, phi, sigma)
- `labelings_opt.txt`: optimal vertex labelings
- `site_graphs_opt/`: optimal migration (site) graphs

A top-level `result.txt` file provides a summary of all pattern sets in tabular format.

---

### `PMH-TR/results/`

Automatically created directory where all PMH-TR outputs are stored.

For each run, subdirectories are organized by pattern set, containing results from that optimization.

---

## PMH-TR Example Runs

### Run PMH-TR

From PMH-TR directory:

```bash
python pmh_tr.py \
  --tree ../data/mcpherson_2016/patient1.tree \
  --labels ../data/mcpherson_2016/patient1.labeling \
  --primary ROv \
  --pattern-set "PS,S"
```
