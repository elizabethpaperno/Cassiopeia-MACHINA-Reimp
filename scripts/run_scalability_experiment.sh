#!/usr/bin/env bash
set -euo pipefail

# Run from repo root: (the directory that contains Cassiopeia/)
# Creates trees + time logs for Greedy / ILP / Hybrid on sub1_train_44_k{10,20,30,32}.txt

ROOT="$(pwd)"
DATA_DIR="$ROOT/Cassiopeia/experiment_results/scalability_experiment/data"

OUT_BASE="$ROOT/Cassiopeia/experiment_results/scalability_experiment"
GREEDY_OUT="$OUT_BASE/greedy"
ILP_OUT="$OUT_BASE/ilp"
HYBRID_OUT="$OUT_BASE/hybrid"
LOG_DIR="$OUT_BASE/time_logs"

mkdir -p "$GREEDY_OUT" "$ILP_OUT" "$HYBRID_OUT" "$LOG_DIR"

# Make sure python can import Cassiopeia/ as a package-ish path
export PYTHONPATH="$ROOT"

# Settings (match what you've been using)
ILP_TIME_LIMIT=1200     # seconds
HYBRID_THRESHOLD=5      # your requested -t 5

# Inputs
KS=(10 20 30 32)

timestamp() { date "+%Y-%m-%d %H:%M:%S"; }

run_greedy () {
  local in_file="$1"
  local tag="$2"
  local out_tree="$GREEDY_OUT/greedy_${tag}.nwk"
  local time_file="$LOG_DIR/time_greedy_${tag}.txt"
  local stdout_file="$LOG_DIR/time_greedy_${tag}.stdout.txt"

  echo "[$(timestamp)] Greedy: $tag"
  : > "$time_file"
  : > "$stdout_file"

  # Greedy script writes "tree.newick" by default; we copy it immediately to avoid overwriting.
  /usr/bin/time -p python "$ROOT/Cassiopeia/Greedy/cassiopeia_greedy.py" "$in_file" \
    > "$stdout_file" 2> >(tee "$time_file" >&2)

  if [ ! -f "$ROOT/tree.newick" ]; then
    echo "ERROR: Greedy did not produce $ROOT/tree.newick for $tag" >&2
    exit 1
  fi

  cp "$ROOT/tree.newick" "$out_tree"
}

run_ilp () {
  local in_file="$1"
  local tag="$2"
  local out_tree="$ILP_OUT/ilp_${tag}.nwk"
  local time_file="$LOG_DIR/time_ilp_${tag}.txt"
  local stdout_file="$LOG_DIR/time_ilp_${tag}.stdout.txt"

  echo "[$(timestamp)] ILP: $tag (limit ${ILP_TIME_LIMIT}s)"
  : > "$time_file"
  : > "$stdout_file"

  /usr/bin/time -p python -m Cassiopeia.ILP.cassiopeia_ilp "$in_file" \
    -o "$out_tree" --debug --time-limit "$ILP_TIME_LIMIT" \
    > "$stdout_file" 2> >(tee "$time_file" >&2)
}

run_hybrid () {
  local in_file="$1"
  local tag="$2"
  local out_tree="$HYBRID_OUT/hybrid_${tag}.nwk"
  local time_file="$LOG_DIR/time_hybrid_${tag}.txt"
  local stdout_file="$LOG_DIR/time_hybrid_${tag}.stdout.txt"

  echo "[$(timestamp)] Hybrid: $tag (threshold ${HYBRID_THRESHOLD})"
  : > "$time_file"
  : > "$stdout_file"

  /usr/bin/time -p python "$ROOT/Cassiopeia/Hybrid/cassiopeia_hybrid.py" "$in_file" \
    -o "$out_tree" -t "$HYBRID_THRESHOLD" \
    > "$stdout_file" 2> >(tee "$time_file" >&2)
}

echo "[$(timestamp)] === START scalability: sub1_train_44_k{10,20,30,32} ==="

for k in "${KS[@]}"; do
  in_file="$DATA_DIR/sub1_train_44_k${k}.txt"
  tag="sub1_train_44_k${k}"

  if [ ! -f "$in_file" ]; then
    echo "Skipping missing input: $in_file"
    continue
  fi

  run_greedy "$in_file" "$tag"
  run_ilp "$in_file" "$tag"
  run_hybrid "$in_file" "$tag"
done

echo "[$(timestamp)] === DONE ==="
echo "Trees saved under:"
echo "  $GREEDY_OUT"
echo "  $ILP_OUT"
echo "  $HYBRID_OUT"
echo "Timing logs saved under:"
echo "  $LOG_DIR"

