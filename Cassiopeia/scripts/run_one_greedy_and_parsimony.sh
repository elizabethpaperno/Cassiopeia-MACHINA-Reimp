#!/usr/bin/env bash
set -euo pipefail

# usage: ./scripts/run_one_greedy_and_parsimony.sh 12
X="$1"

IN="data/train/sub1_train_${X}.txt"
OUT_TREE="experiment_results/general_full_dataset/greedy/greedy_sub1_train_${X}.nwk"

mkdir -p "$(dirname "$OUT_TREE")"

echo "Running Greedy on $IN"
python Greedy/cassiopeia_greedy.py "$IN" -o "$OUT_TREE" >/dev/null

echo "Saved tree: $OUT_TREE"
echo "First 200 chars:"
head -c 200 "$OUT_TREE"; echo
echo

echo "Parsimony (sum first branch length per edge):"
python compute_parsimony_greedy.py "$OUT_TREE"

