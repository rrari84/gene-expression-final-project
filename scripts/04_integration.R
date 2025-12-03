#!/usr/bin/env Rscript

# 00_setup_environment.R
# -----------------------
# Sets up a reproducible R environment for the project using renv.
# - Should be run from the project root directory (where renv.lock will live).
# - Assumes you're already on the HPC with:
#     module purge
#     module load R/4.4.1
#
# Usage (from shell, in project root):
#   Rscript scripts/00_setup_environment.R
#
# After this finishes, all core packages should be installed and loadable.

message("=== 00_setup_environment.R: Starting environment setup ===")

## 1. Install renv if needed ----
if (!requireNamespace("renv", quietly = TRUE)) {
  message("[renv] Not found. Installing renv from CRAN...")
  install.packages("renv", repos = "https://cloud.r-project.org")
} else {
  message("[renv] Already installed.")
}

library(renv)

## 2. Initialize or restore renv project ----
if (!file.exists("renv.lock")) {
  message("[renv] No renv.lock found. Initializing a fresh renv project...")
  # bare = TRUE avoids snapshotting a bunch of random global packages
  renv::init(bare = TRUE)
} else {
  message("[renv] renv.lock found. Restoring packages from lockfile...")
  renv::restore()
}

# Make sure the project library is active
renv::activate()

## 3. Define package lists ----

# Core CRAN packages
cran_packages <- c(
  "tidyverse",    # dplyr, ggplot2, etc.
  "data.table",
  "here",         # convenient paths
  "readr",
  "survival",
  "survminer",
  "ggrepel",
  "patchwork"     # combining ggplots
)

# Core Bioconductor packages
bioc_packages <- c(
  "GEOquery",         # download GEO datasets
  "limma",            # microarray / voom DE
  "edgeR",            # count-based DE (optional)
  "DESeq2",           # RNA-seq DE
  "sva",              # batch correction (if needed)
  "clusterProfiler",  # GO/KEGG enrichment
  "org.Hs.eg.db",     # human gene annotations
  "ReactomePA",       # Reactome pathway enrichment
  "oligo",            # microarray CEL preprocessing (RMA, etc.)
  "Biobase",          # ExpressionSet / pData / exprs
  "AnnotationDbi"     # for mapping probes -> gene symbols (platform-specific)
  # NOTE: platform-specific pd-* / *db packages can be added later once we know the array
)

## 4. Install CRAN packages (inside renv) ----

message("=== Installing CRAN packages (if missing) ===")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("[CRAN] Installing %s ...", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    message(sprintf("[CRAN] %s already installed.", pkg))
  }
}

## 5. Install Bioconductor packages (inside renv) ----

message("=== Installing Bioconductor packages (if missing) ===")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("[BiocManager] Not found. Installing from CRAN...")
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

library(BiocManager)

BiocManager::install(
  bioc_packages,
  ask    = FALSE,  # don't prompt
  update = FALSE   # avoid updating everything globally
)

## 6. Verify that all packages load ----

message("=== Verifying that all required packages can be loaded ===")

all_pkgs <- c(cran_packages, bioc_packages)

failed <- character(0)

for (pkg in all_pkgs) {
  message(sprintf("[check] Loading %s ...", pkg))
  ok <- require(pkg, character.only = TRUE, quietly = TRUE)
  if (!ok) {
    warning(sprintf("!!! Failed to load package: %s", pkg))
    failed <- c(failed, pkg)
  }
}

if (length(failed) > 0) {
  message("=== WARNING ===")
  message("The following packages failed to load:")
  print(failed)
  message("Please check the installation messages above and try again.")
} else {
  message("All required packages loaded successfully.")
}

## 7. Snapshot environment ----

message("=== Taking renv snapshot ===")
renv::snapshot(prompt = FALSE)

message("=== 00_setup_environment.R: Finished successfully ===")
