# Compute environment

Base image: R 4.3 (Code Ocean "R" base image, Ubuntu 22.04).

System packages required for compiling the C++ backend via `Rcpp::sourceCpp()`:

- build-essential
- gfortran
- libblas-dev, liblapack-dev

R packages are installed by `postInstall`. PLINK v1.9 is not required inside the
capsule: genotype extraction and allele alignment are upstream steps performed on
the restricted study data, and the capsule reads the resulting `.raw` file directly.
