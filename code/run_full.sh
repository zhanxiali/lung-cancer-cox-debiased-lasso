#!/usr/bin/env bash
# Full analytic settings as used for the published results: 30 projection
# multipliers on a log grid from 0.001 to 5, five-fold cross-validation.
#
# Runtime scales with the number of multipliers, the number of folds, and p^2.
# On the simulated dataset (n = 400, p ~ 85) expect roughly one to two hours on
# a single core. On the study data (n = 979-2,322, p = 56-79) the per-subgroup
# analyses took several hours each.
#
# The study data are not distributed with this capsule; see the Data
# Availability Statement in the paper. This script therefore runs the full
# settings on the simulated dataset.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p ../results
Rscript demo_pipeline.R n_multi=30 n_cv=5 out=../results 2>&1 | tee ../results/full_log.txt
echo
echo "Outputs:"
ls -1 ../results
