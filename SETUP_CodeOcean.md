# Setting up the Code Ocean capsule

Code Ocean runs a capsule as three fixed locations: `/code` holds this repository,
`/data` holds read-only inputs, and `/results` collects whatever the run writes.
Pressing **Reproducible Run** executes `/code/run`.

Follow the steps below in order. The whole process takes about half an hour, most
of which is waiting for the environment to build.

---

## 1. Get the files into the capsule

**If the repository is already on GitHub** — in the capsule, open the **Code** pane,
choose **Import from Git**, and paste the repository URL. Everything lands under
`/code` with the directory structure intact.

**If it is not on GitHub yet** — drag the whole folder onto the Code pane. Keep the
subdirectories (`core/`, `code/`, `per_cohort/`, `meta/`, `sensitivity/`,
`figures/`, `example/`); do not flatten them.

Check afterwards that `/code/run` exists at the top level, not inside a
subdirectory. Code Ocean will not find it otherwise.

## 2. Configure the environment

Open the **Environment** pane.

- Base image: **R**, version **4.3** or later.
- Under **Packages**, add: `quadprog`, `Rcpp`, `RcppArmadillo`, `glmnet`,
  `survival`, `dplyr`, `ggplot2`, `metafor`.
- Under **System (apt-get)**, add: `build-essential`, `gfortran`, `libblas-dev`,
  `liblapack-dev`. These are required to compile `core/sim_univLib.cpp`; without
  them the run fails at `sourceCpp()`.

If you prefer to install packages by script instead of through the GUI, paste the
contents of `environment/postInstall` into the **postInstall** box.

Click **Build** and wait. The first build takes ten to twenty minutes; later runs
reuse the cached image.

## 3. Decide what goes in `/data`

The individual-level genotype and clinical data are restricted and must not be
uploaded. Leave `/data` empty. `run` detects this and uses the simulated dataset
in `example/`, which has the same structure.

If you ever want to run the capsule privately on real data, upload the files as
`/data/genotypes.raw` and `/data/clinical.csv`; `run` picks them up automatically.
Do not publish the capsule in that state.

## 4. Do a test run

Press **Reproducible Run**. Expect roughly:

```
Inputs: simulated dataset in example/
Mode: demo (multipliers=8, folds=3)
Compiling C++ routines
Reading simulated data
  n = 400, SNPs = 77, events = 240
  design matrix: n = 400, p = 85
Penalized Cox fit
  lambda.min = 0.0xxxx, nonzero coefficients = xx
Tuning the projection multiplier (8 values, 3 folds)
  ...
Final debiased estimation
```

`/results` should then contain `demo_debiased_results.csv`,
`demo_session_info.txt`, and `run_log.txt`.

To run the published settings instead, add `RUN_MODE=full` as an environment
variable in the Environment pane and run again. This takes considerably longer.

## 5. Write the capsule metadata

In the capsule **Metadata** pane:

- **Title** — the manuscript title.
- **Description** — state plainly what the capsule does and does not do, for
  example:

  > Analysis code for the Cox debiased lasso pipeline used to evaluate lung
  > function–associated variants against overall survival in non–small cell lung
  > cancer. The capsule executes the full pipeline on a simulated dataset of
  > identical structure. Individual-level genotype and clinical data are governed
  > by institutional data use agreements and are not distributed; see the Data
  > Availability Statement in the paper. The capsule documents the analytic
  > implementation rather than regenerating the reported estimates.

- **License** — MIT, matching `LICENSE`.

## 6. Link the capsule to the manuscript

Return to the submission system and use the Code Ocean link in the decision
letter. Linking from there associates the capsule with the manuscript record and
assigns it a DOI on publication. Then add that DOI to the Code Availability
Statement in the paper.

---

## If the run fails

**`sourceCpp` error, compiler not found** — the apt packages in step 2 were not
added. Add them and rebuild.

**`could not find function "neg_loglik_functions_cpp_ext"`** — `sim_univLib.cpp`
compiled but exported nothing, usually because the file did not upload completely.
Check its size in the Code pane.

**`File not found: /code/core/sim_univLib.cpp`** — the repository was flattened on
upload, or `run` sits inside a subdirectory. Restore the structure from step 1.

**`solve.QP` reports a non-positive-definite matrix** — the simulated dataset is
too small for the number of covariates. Regenerate it larger:
`python example/make_example_data.py --n 800`.

**Run exceeds the time limit** — stay in demo mode, or lower the settings further
by editing `N_MULTI` and `N_CV` in `run`.
