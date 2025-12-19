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
Each subdirectory typically includes:
- `.tree` files: clone tree edge lists (`parent child` per line)
- `.labeling` files: mappings from leaf nodes to anatomical sites (`leaf site` per line)
- `.colormap` files: mappings from anatomical sites to colors (used for visualization)

---

### `visualizations/`
Contains scripts for visualizing PMH, PMH-TI, and PMH-TR outputs.

- `draw_tree.py`:  
  Given:
  - a migration (site) graph file (e.g. produced by `pmh.py`) in the format:
    ```
    src_site  dst_site  weight
    ```
  - a colormap file (`site color`)
  
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

For each run, a subdirectory is created containing:
- `summary.txt`
- `labelings_all_mp.txt`
- `labelings_opt.txt`
- `site_graphs_all_mp/`
- `site_graphs_opt/`

---

## Example Runs

### Run PMH
From PMH directory:
```bash
python pmh.py \
  --tree ../data/mcpherson_2016/patient1.tree \
  --labels ../data/mcpherson_2016/patient1.labeling \
  --primary ROv \
  --pattern-set "PS,S"

### Now Get Associated Visualization
From 
```bash 
python visualizations/draw_tree.py \
  --site-graph PMH/results/ROv_PS,S_patient1/site_graphs_opt/site_graph_0.txt \
  --colors data/mcpherson_2016/patient1.colormap\
  --output PMH/results/ROv_PS,S_patient1/site_graphs_opt/site_graph_0.png \
  --primary ROv
