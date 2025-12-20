#!/usr/bin/env bash
set -euo pipefail

ROOT=$(pwd)
PROJECT_ROOT=$(dirname "$ROOT")   # parent of Cassiopeia
GREEDY_DIR="$ROOT/experiment_results/general_full_dataset/greedy"
PARSIMONY_SCRIPT="$ROOT/compute_parsimony.py"

echo "Greedy parsimony values"
echo "======================="

for tree in "$GREEDY_DIR"/greedy_*.nwk; do
    base=$(basename "$tree")
    dataset=${base#greedy_}
    dataset=${dataset%.nwk}

    echo
    echo "Dataset: $dataset"
    PYTHONPATH="$PROJECT_ROOT" python "$PARSIMONY_SCRIPT" "$tree"
done

