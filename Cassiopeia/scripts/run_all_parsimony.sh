#!/usr/bin/env bash
set -euo pipefail

ROOT="Cassiopeia/experiment_results"
PARSER="Cassiopeia/compute_parsimony.py"

echo "=============================="
echo "Parsimony results"
echo "=============================="
echo

run_parsimony_dir () {
  local label="$1"
  local dir="$2"

  if [ ! -d "$dir" ]; then
    return 0
  fi

  # Skip if no .nwk files
  shopt -s nullglob
  local files=("$dir"/*.nwk)
  shopt -u nullglob
  if [ ${#files[@]} -eq 0 ]; then
    return 0
  fi

  echo "--- $label ---"
  PYTHONPATH=. python "$PARSER" "${files[@]}"
  echo
}

echo "=== General full dataset ==="
echo
run_parsimony_dir "greedy" "$ROOT/general_full_dataset/greedy"
run_parsimony_dir "hybrid" "$ROOT/general_full_dataset/hybrid"
run_parsimony_dir "ilp"    "$ROOT/general_full_dataset/ilp"

echo "=== Scalability experiment ==="
echo
run_parsimony_dir "greedy" "$ROOT/scalability_experiment/greedy"
run_parsimony_dir "hybrid" "$ROOT/scalability_experiment/hybrid"
run_parsimony_dir "ilp"    "$ROOT/scalability_experiment/ilp"

