# Lung function–associated genetic variants and overall survival in NSCLC

Analysis code accompanying:

> Li Z, Zhang Z, Zhao Y, Su L, Li Y, Christiani DC. Lung Function–Associated Genetic Variants and
> Overall Survival in Non–Small Cell Lung Cancer: A Multi-Platform Analysis. *Cancer Epidemiology,
> Biomarkers & Prevention* (in press).

This repository contains the code used to produce every result in the paper. It contains **no
individual-level data**: the genotype and clinical data are governed by institutional data use
agreements (see *Data availability* below).

---

## What the pipeline does

A prespecified panel of 77 lung function GWAS–identified SNPs is evaluated against overall
survival in three genotyping subgroups of the Boston Lung Cancer Survival Cohort. Within each
subgroup, a Cox debiased lasso jointly estimates conditional per-allele log-hazard ratios for all
candidate SNPs together with clinical covariates and ancestry principal components. The debiased
estimates are then harmonized to a common effect allele and combined across subgroups by
inverse-variance–weighted fixed-effect meta-analysis.

```
per-subgroup:   genotypes + covariates
                      |
                L1-penalized Cox fit  (glmnet, lambda.min)
                      |
                score-based debiasing (quadratic programming, C++ backend)
                      |
                beta, SE, P per covariate
                      |
across subgroups:  allele harmonization  ->  IVW fixed-effect meta-analysis
                      |
                BH-FDR  ->  Tier 1 / Tier 2 evidence hierarchy
```

## Repository layout

```
core/
  applydebiaslasso.R                    core debiased lasso function
  sim_univLib.cpp                       C++ routines (score and score-squared matrices)

per_cohort/
  HSPH_cox_debiased_lasso_analysis.R    subgroup analyses; each writes
  Oncoarray_cox_debiased_lasso_analysis.R   <cohort>_cox_debiased_results_all.csv
  MGH_cox_debiased_lasso_analysis.R

meta/
  Meta_analysis_three_cohorts.R         allele harmonization + IVW meta-analysis + FDR
  Meta_sensitivity_NO_TREATMENT.R       sensitivity analysis excluding treatment covariates

sensitivity/
  Cohort_debiased_lasso_lambda1se.R     full pipeline rerun under lambda.1se  (Table S7)
  Meta_lambda1se_comparison.R           lambda.min vs lambda.1se comparison   (Table S7)
  PRS_LOCO_analysis.R                   leave-one-cohort-out PRS             (Table S8, Fig. S3)
  Supplementary_Analyses.R              EAF concordance, Schoenfeld residuals (Tables S3, S6; Fig. S1)

figures/
  Generate_Figure1_Workflow.py          Figure 1
  Regenerate_Figures_2_3.R              Figures 2 and 3
  Generate_Figure4_I2.py                Figure 4

example/
  example_genotypes.raw                 small synthetic dataset so the pipeline runs end to end
  example_clinical.csv                  (synthetic; not derived from any participant)
```

## Requirements

R (>= 4.3) with `glmnet`, `survival`, `quadprog`, `Rcpp`, `RcppArmadillo`, `dplyr`, `ggplot2`,
`metafor`; PLINK v1.9 for genotype extraction and allele alignment. Python (>= 3.10) with
`matplotlib` for Figures 1 and 4. A working C++ toolchain is required to compile
`core/sim_univLib.cpp` via `Rcpp::sourceCpp()`.

Exact package versions used for the published results are recorded in `sessionInfo.txt`.

## How to run

No paths are hard-coded. Each script prompts for its inputs interactively, so the pipeline can be
run from any directory:

```r
source("per_cohort/HSPH_cox_debiased_lasso_analysis.R")
```

You will be asked to select, in order, the C++ library, the aligned PLINK `.raw` genotype file, and
the cleaned clinical CSV. The synthetic files in `example/` can be used to verify that the
environment is set up correctly before running on real data.

Typical order: per-cohort analyses -> `meta/Meta_analysis_three_cohorts.R` -> the scripts in
`sensitivity/` -> the scripts in `figures/`.

## Input format

Genotypes are read from PLINK `.raw` files (`--recode A`) after effect alleles have been aligned
across platforms with `--a1-allele` using the HSPH subgroup as reference. Clinical files are CSVs
keyed on `IID` and must contain `OS_month`, `DEAD`, `AGE`, `SEX`, `smoksort`, `early_late`,
`surgery`, `chemotx`, `RADS`, and `pc1`–`pc3`.

Note on allele coding: for a subset of variants the per-subgroup `.raw` coding is opposite to the
harmonized reference coding. `meta/Meta_analysis_three_cohorts.R` detects and corrects this before
pooling; the per-subgroup result files retain each subgroup's own coding. Downstream scripts that
combine subgroups apply the same correction (see Supplementary Table S3 and Figure 1).

## Data availability

Individual-level genotype and clinical data are not included here and cannot be redistributed.
Genotype data for the HSPH subgroup are available through dbGaP (accession `phs000093.v2.p2`) and
OncoArray genotype data through the ILCCO-OncoArray deposit (accession `phs001273.v3.p2`). Data for
the MGH low-pass whole-genome sequencing subgroup contain protected health information and are
available only under an institutional data use agreement; requests should be directed to the
corresponding author. Summary-level meta-analysis results for all 65 meta-analyzed SNPs are
provided as Supplementary Table S5 of the paper.

## Citation

If you use this code, please cite the paper above and this archived release:

> [Author list]. Analysis code for "Lung Function–Associated Genetic Variants and Overall Survival
> in Non–Small Cell Lung Cancer". Zenodo; 2026. doi:10.5281/zenodo.XXXXXXX

## License

MIT (see `LICENSE`). The license covers the code only; it confers no rights to any data.

## Contact

David C. Christiani, Department of Environmental Health, Harvard T.H. Chan School of Public Health
(dchris@hsph.harvard.edu).

## Running without prompts

Every script selects its inputs at run time rather than from a hard-coded path.
In an interactive R session a file dialog appears; otherwise the path is read from
the console. To run unattended, set the corresponding environment variable:

| Variable      | Used by                                   | Points to                              |
|---------------|-------------------------------------------|----------------------------------------|
| `SIM_UNIVLIB` | core, per-cohort scripts                  | `core/sim_univLib.cpp`                 |
| `GENO_RAW`    | per-cohort scripts                        | aligned PLINK `.raw` genotype file     |
| `CLIN_CSV`    | per-cohort scripts                        | cleaned clinical CSV                   |
| `SNP_POS`     | per-cohort scripts                        | candidate SNP position file            |
| `RES_HSPH`    | `meta/Meta_analysis_three_cohorts.R`      | HSPH debiased results CSV              |
| `RES_ONCO`    | `meta/Meta_analysis_three_cohorts.R`      | OncoArray debiased results CSV         |
| `RES_MGH`     | `meta/Meta_analysis_three_cohorts.R`      | MGH debiased results CSV               |

```bash
export SIM_UNIVLIB=$PWD/core/sim_univLib.cpp
export GENO_RAW=$PWD/example/example_genotypes.raw
export CLIN_CSV=$PWD/example/example_clinical.csv
Rscript per_cohort/HSPH_cox_debiased_lasso_analysis.R
```

## Before publishing the repository

Run `bash scrub_check.sh` from the repository root. It searches every script for
working directories, cluster paths, hostnames, user accounts, and the IRB protocol
number, and reports the file and line of anything it finds. The check must print
`Clean` before the repository is made public.
