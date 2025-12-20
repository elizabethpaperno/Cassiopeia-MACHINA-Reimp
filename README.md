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

## Cassiopea 

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

- `draw_tree.py`:  
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
- reads a clone tree and leaf labeling,
- runs PMH under a specified pattern set,
- outputs all maximum parsimony labelings and optimal solutions,
- constructs a migration graphs.

#### Inputs
- `--tree`: path to a `.tree` file
- `--labels`: path to a `.labeling` file
- `--primary`: primary anatomical site
- `--pattern-set`: one of `PS`, `PS,S`, `PS,S,M`, `PS,S,M,R`

---

### `PMH/results/`
Automatically created directory where all PMH outputs are stored.

For each run, a subdirectory with name {primary site}_{pattern}_{tree/labeling file name} is created containing:
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
python visualizations/draw_tree.py \
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
From PMH directory:
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
python visualizations/draw_tree.py \
  --site-graph PMH/results_pmh_ti/site_graph_pmh_ti.txt \
  --colors data/hoadley_2016/A7.colormap \
  --output ../figures/A7_breast_site_graph.png \
  --primary breast
```

## PMH-TR
