#!/bin/bash
##------------------------------------
## pySCENIC step 2/3: GRNBoost2 -> cisTarget -> AUCell
## Author: Peng Yang
##
## Run from the repository root:  bash TF_activity/run_pyscenic.sh
## Requires pySCENIC 0.12.1 (see README for the pinned conda environment) and
## the hg38 cisTarget v10 databases in Data/cisTarget_db:
##   hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
##   hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
##   motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
## Optional: PYSCENIC_ENV=/path/to/conda/env bash TF_activity/run_pyscenic.sh
##------------------------------------
set -euo pipefail

WORK="TF_activity/pyscenic"
DB="$(cd Data/cisTarget_db && pwd)"     # absolute path, resolved before cd
NCORES=4

## dask: disable spill/pause/terminate (avoids the distributed spill.py assertion)
export DASK_DISTRIBUTED__WORKER__MEMORY__TARGET=False
export DASK_DISTRIBUTED__WORKER__MEMORY__SPILL=False
export DASK_DISTRIBUTED__WORKER__MEMORY__PAUSE=False
export DASK_DISTRIBUTED__WORKER__MEMORY__TERMINATE=False

## ---- pySCENIC command (conda environment first, if given) ----
if [ -n "${PYSCENIC_ENV:-}" ] && [ -x "$PYSCENIC_ENV/bin/pyscenic" ]; then
  export PATH="$PYSCENIC_ENV/bin:$PATH"
  PYSCENIC="$PYSCENIC_ENV/bin/pyscenic"
  PY="$PYSCENIC_ENV/bin/python"
elif command -v pyscenic >/dev/null 2>&1; then
  PYSCENIC="pyscenic"
  PY="python3"
else
  echo "ERROR: pyscenic not found (set PYSCENIC_ENV or add pyscenic to PATH)."
  exit 1
fi

cd "$WORK"

## ---- 1. loom file from the exported sparse matrix ----
if [ ! -f expr.loom ]; then
"$PY" - <<'PY'
import scipy.io, numpy as np, loompy as lp
m = scipy.io.mmread("counts.mtx").tocsc()          # genes x cells
genes = [l.strip() for l in open("genes.txt")]
cells = [l.strip() for l in open("barcodes.txt")]
lp.create("expr.loom", m, {"Gene": np.array(genes)}, {"CellID": np.array(cells)})
print("loom built:", m.shape[0], "genes x", m.shape[1], "cells")
PY
fi

D5="$DB/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
D10="$DB/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
ANN="$DB/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

## ---- 2. co-expression modules (GRNBoost2) ----
## re-use adj.tsv only if it is non-empty and newer than expr.loom
if [ -s adj.tsv ] && [ adj.tsv -nt expr.loom ]; then
  echo "adj.tsv exists - skipping GRN step."
else
  rm -f adj.tsv reg.csv auc.csv
  $PYSCENIC grn expr.loom tfs.txt -o adj.tsv --num_workers "$NCORES" --seed 777
  [ -s adj.tsv ] || { echo "ERROR: empty adj.tsv"; exit 1; }
fi

## ---- 3. regulons (cisTarget motif enrichment) ----
if [ -s reg.csv ] && [ reg.csv -nt adj.tsv ]; then
  echo "reg.csv exists - skipping ctx step."
else
  rm -f reg.csv auc.csv
  $PYSCENIC ctx adj.tsv "$D5" "$D10" \
    --annotations_fname "$ANN" \
    --expression_mtx_fname expr.loom \
    --mode "dask_multiprocessing" \
    --output reg.csv --num_workers "$NCORES"
  [ -s reg.csv ] || { echo "ERROR: empty reg.csv"; exit 1; }
fi

## ---- 4. regulon activity per cell (AUCell) ----
$PYSCENIC aucell expr.loom reg.csv --output auc.csv --num_workers "$NCORES"

echo "Done. Regulon activity written to $WORK/auc.csv"
