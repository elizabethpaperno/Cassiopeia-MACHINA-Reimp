#!/usr/bin/env bash
set -euo pipefail

# Run from: Cassiopeia-MACHINA-Reimp/
# Uses threshold t = 5 for all runs

PYTHON_BIN=python
HYBRID_SCRIPT=Cassiopeia/Hybrid/cassiopeia_hybrid.py
THRESHOLD=5

GEN_DATA_DIR=Cassiopeia/data/train
GEN_OUT_DIR=Cassiopeia/experiment_results/general_full_dataset
GEN_HYBRID_DIR=${GEN_OUT_DIR}/hybrid
GEN_TIME_DIR=${GEN_OUT_DIR}/time_logs

SCALE_DATA_DIR=Cassiopeia/experiment_results/scalability_experiment/data
SCALE_OUT_DIR=Cassiopeia/experiment_results/scalability_experiment
SCALE_HYBRID_DIR=${SCALE_OUT_DIR}/hybrid
SCALE_TIME_DIR=${SCALE_OUT_DIR}/time_logs

mkdir -p "${GEN_HYBRID_DIR}" "${GEN_TIME_DIR}"
mkdir -p "${SCALE_HYBRID_DIR}" "${SCALE_TIME_DIR}"

echo "=== Running HYBRID on general_full_dataset ==="

for txt in ${GEN_DATA_DIR}/sub1_train_*.txt; do
    base=$(basename "${txt}" .txt)

    out_nwk="${GEN_HYBRID_DIR}/hybrid_${base}.nwk"
    time_log="${GEN_TIME_DIR}/time_hybrid_${base}.txt"

    if [[ -f "${out_nwk}" ]]; then
        echo "SKIP ${base} (already exists)"
        continue
    fi

    echo "Running HYBRID on ${base}"

    PYTHONPATH=. /usr/bin/time -p \
        ${PYTHON_BIN} ${HYBRID_SCRIPT} "${txt}" \
        -o "${out_nwk}" \
        -t ${THRESHOLD} \
        2> >(tee "${time_log}" >&2)
done

echo "=== Running HYBRID on scalability_experiment ==="

for txt in ${SCALE_DATA_DIR}/sub1_train_7_k*.txt; do
    base=$(basename "${txt}" .txt)

    out_nwk="${SCALE_HYBRID_DIR}/hybrid_${base}.nwk"
    time_log="${SCALE_TIME_DIR}/time_hybrid_${base}.txt"

    if [[ -f "${out_nwk}" ]]; then
        echo "SKIP ${base} (already exists)"
        continue
    fi

    echo "Running HYBRID on ${base}"

    PYTHONPATH=. /usr/bin/time -p \
        ${PYTHON_BIN} ${HYBRID_SCRIPT} "${txt}" \
        -o "${out_nwk}" \
        -t ${THRESHOLD} \
        2> >(tee "${time_log}" >&2)
done

echo "=== DONE: All hybrid trees generated ==="

