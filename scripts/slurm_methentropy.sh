#!/bin/bash
#SBATCH --job-name=methentropy
#SBATCH --output=methentropy_%A_%a.out
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=0-2
# Edit the paths below. One array task per sample.

set -euo pipefail
module purge
ml miniconda
conda activate /path/to/your/env          # env where `pip install -e .` was run

IN_DIR=/path/to/rlm_output                 # full RLM entropy files (23 columns)
OUT_DIR=/path/to/methentropy_output
mkdir -p "$OUT_DIR"

SAMPLES=(mef_p_1 mef_p0 mef_p5)
S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

methentropy calibrate -i "$IN_DIR/${S}_entropy.bed" --windows 2000 \
    > "$OUT_DIR/${S}_calibration.json"

methentropy run -i "$IN_DIR/${S}_entropy.bed" -o "$OUT_DIR/${S}_methentropy.tsv.gz" \
    --nsims 1000 --seed 1 --min-coverage 10 --max-coverage 1000 --chunksize 2000000
