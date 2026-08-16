#!/usr/bin/env bash
# Reduced-scale run on the simulated dataset. Completes in a few minutes.
# Used to verify that the environment is configured correctly and that the
# pipeline executes end to end.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p ../results
Rscript demo_pipeline.R n_multi=8 n_cv=3 out=../results 2>&1 | tee ../results/demo_log.txt
echo
echo "Outputs:"
ls -1 ../results
