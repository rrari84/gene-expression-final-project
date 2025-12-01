# Project Gameplan 

## 1. Big Picture (What are we trying to learn?)

We're using **public gene expression data** from breast cancer and normal breast tissue to answer three main questions:

1. **Which genes behave differently in breast cancer vs normal?**  
   → "Differentially expressed genes" (DEGs).

2. **What biological processes and pathways are most disturbed?**  
   → Cell cycle? DNA repair? Immune response? Hormones?

3. **Are any of these genes linked to patient outcomes or good drug targets?**  
   → Survival analysis + network analysis + (maybe) simple drug repurposing.

The whole project is about turning large tables of numbers into a **biological story** about how breast cancer cells differ from normal cells.

---

## 2. The datasets (What data are we using?)

We will mainly use:

- **GEO breast cancer datasets**  
  - Each dataset has:
    - Gene expression (microarray or RNA-seq)
    - Sample labels: tumor vs normal (or different clinical groups)

- **(Optional/if time) TCGA BRCA dataset**  
  - Gene expression + clinical data (overall survival, etc.)

We don't need to collect any new samples. We're re-analyzing **existing public data**.

---

## 3. The pipeline (Step-by-step in plain English)

### Step 1 – Get and clean the data

**Goal:** Turn raw GEO datasets into clean expression tables we can analyze.

What we do:

- Download expression data + sample information from GEO (in R).
- Make sure:
  - All samples have clear labels (e.g., tumor vs normal).
  - Gene names are standardized (e.g., HGNC symbols).
- Log-transform / normalize as needed so values are on a comparable scale.

**Biology translation:**  
We're organizing the data so that each row is a gene and each column is a patient/sample, with a clear label of "cancer" or "normal".

---

### Step 2 – Find differentially expressed genes (DEGs) per dataset

**Goal:** For each dataset, find genes that are significantly **up** or **down** in cancer vs normal.

What we do:

- For each dataset separately:
  - Compare cancer samples vs normal samples.
  - For each gene, compute:
    - How much higher or lower it is in cancer (**log2 fold change**)
    - A **p-value** and **adjusted p-value (FDR)** for significance.
- Choose "DEGs" based on thresholds (e.g., |log2FC| > 1 and FDR < 0.05).
- Make **volcano plots** to visualize:
  - X-axis: change in expression
  - Y-axis: significance

**Biology translation:**  
We're answering: "Which genes are turned up or down in breast cancer?" separately in each study.

---

### Step 3 – Combine multiple datasets (robust DEGs)

**Goal:** Find genes that are **consistently** different across several datasets, not just in one.

What we do:

- Compare DEG lists from each dataset.
- Define "robust DEGs" as:
  - Genes that are significantly up in cancer in **several datasets**, or
  - Genes that are significantly down in **several datasets**.
- Optionally use simple meta-analysis (combining p-values) to rank genes.

**Biology translation:**  
We're looking for genes whose abnormal expression **keeps showing up across many experiments**, making them more trustworthy as real cancer-related genes.

---

### Step 4 – Functional enrichment (GO / pathways)

**Goal:** Turn a long list of genes into a **biological story**.

What we do:

- Take the robust **up-regulated** genes and run:
  - GO Biological Process enrichment
  - Pathway enrichment (KEGG, Reactome, etc.)
- Do the same for **down-regulated** genes.
- Identify top enriched terms (e.g. "cell cycle", "DNA repair", "immune response").

**Biology translation:**  
Instead of saying "we saw 200 weird genes," we can say things like:
> "Breast cancer samples show up-regulation of cell cycle and proliferation pathways and down-regulation of immune response pathways."

---

### Step 5 – Protein–protein interaction (PPI) network and hub genes

**Goal:** See how our genes interact as proteins and identify **key players (hubs)**.

What we do:

- Take robust DEGs and input them to **STRING** to build a PPI network.
- Load the network into **Cytoscape** to:
  - Visualize interaction clusters (modules).
  - Compute network measures (degree, centrality) to find **hub genes**.

**Biology translation:**  
We're asking:
- "How do the proteins made by these genes talk to each other?"
- "Which ones are at the center of the conversation?"
  
Hub genes that are strongly up-regulated in cancer are interesting as:
- Possible biomarkers
- Possible therapeutic targets

---

### Step 6 – Survival / prognosis analysis (if data available)

**Goal:** See if gene expression is linked to how patients actually do.

What we do (likely using TCGA BRCA):

- Pick a small set of important genes (e.g., top hubs).
- For each gene:
  - Divide patients into "high expression" vs "low expression" groups.
  - Make Kaplan–Meier survival curves.
  - Use a log-rank test to see if survival differs.

**Biology translation:**  
We're asking:
> "Do patients whose tumors have high expression of this gene live longer or shorter than those with low expression?"

If yes, the gene might be a **prognostic marker**.

---

### Step 7 – (Stretch goal) Simple drug repurposing angle

**Goal:** See if any known drugs might push gene expression back toward "normal."

What we do (if there's time):

- Use our up/down gene signature to query a drug perturbation resource (e.g., Connectivity Map / LINCS).
- Look for compounds whose effects on gene expression are **opposite** to the breast cancer pattern.

**Biology translation:**  
We're asking:
> "Are there existing drugs that make cells move in the opposite direction of the cancer signature?"  

Even a rough list of candidate drugs is useful for discussion.

---

## 4. How the HPC fits in (short version)

We use the university's HPC as a **big shared machine** to:

- Store large expression datasets and intermediate results.
- Run heavier scripts (e.g., DEG analysis across all datasets) without overloading laptops.
- Keep one consistent R environment for the whole group.

Practical setup:

- Shared project folder under `/projects/...`
- Everyone runs:
  ```bash
  module purge
  module load R/4.4.1
  ```
- Then uses the same R scripts / notebooks from the repo.

---

## 5. What each of us can focus on (bio-friendly roles)

Even without heavy coding, there are important biology-focused tasks:

### Pathway & literature leads

- Interpret enriched GO/pathway results:
  - Which pathways make sense for breast cancer?
  - Which are surprising or new?
- Look up our top genes in the literature.

### Network biology lead

- Explore PPI networks in Cytoscape.
- Identify interesting clusters:
  - Cell cycle cluster?
  - Immune/chemokine cluster?
- Connect clusters to known cancer hallmarks.

### Clinical relevance lead

- Interpret survival plots (if we do TCGA analysis).
- Ask:
  - Are these genes already known prognostic markers?
  - Are they associated with specific subtypes?

### Writing & storytelling

- Turn the technical results into:
  - Clear figures with readable labels.
  - A narrative that flows from:
    - DEGs → pathways → networks → possible clinical/therapeutic angles.

---

## 6. Quick glossary

| Term | Definition |
|------|------------|
| **Gene expression** | How much a gene is "on" or "off" (mRNA level). |
| **DEG (Differentially Expressed Gene)** | A gene with significantly different expression between groups (e.g., tumor vs normal). |
| **FDR (False Discovery Rate)** | p-value adjusted for many tests; helps control for false positives when testing thousands of genes. |
| **GO (Gene Ontology)** | Standardized terms for biological processes, molecular functions, and cellular components. |
| **Pathway** | A set of genes/proteins that work together in a specific process. |
| **PPI (Protein–Protein Interaction) network** | A map showing which proteins physically or functionally interact. |
| **Hub gene** | A protein with many connections in the PPI network; often biologically important. |
| **Kaplan–Meier curve** | Plot showing the proportion of patients surviving over time in different groups. |
| **Connectivity / drug repurposing analysis** | Matching our gene signature to drugs that reverse or mimic it. |

---

## 7. One-sentence summary of the project

We are combining multiple public breast cancer gene expression datasets to find robust, biologically meaningful changes in gene activity, understand which pathways and networks are disrupted, and (if data allows) connect these changes to patient survival and potential therapeutic strategies.

---

## 8. R scripts & acceptance criteria

This section explains what each R file in the `scripts/` folder should do, and how we'll know it's "done and working."

---

### 8.1 `00_setup_environment.R`

#### Purpose (plain English)
Set up the R environment for the whole project: install required packages and initialize renv so everyone uses the same versions.

#### What it should do

1. Install `renv` if not already installed.
2. Initialize `renv` in the project (or restore from an existing lockfile).
3. Install all required CRAN and Bioconductor packages, for example:
   - `tidyverse`, `data.table`
   - `GEOquery`
   - `limma`, `edgeR`, `DESeq2`
   - `clusterProfiler`, `org.Hs.eg.db`, `ReactomePA`
   - `survival`, `survminer`

#### Inputs
- None (other than an internet connection and R).

#### Outputs
- `renv.lock` file in the project root.
- Packages installed in the local renv library.

#### Acceptance criteria
- ✅ Script runs from start to finish without errors.
- ✅ After running, the following succeed in R without "package not found":
  - `library(DESeq2)`, `library(limma)`, `library(GEOquery)`, `library(clusterProfiler)`, `library(survival)`
- ✅ `renv::status()` reports no major problems.

---

### 8.2 `01_download_and_prepare_data.R`

#### Purpose (plain English)
Download breast cancer datasets from GEO and turn them into clean expression + sample annotation tables.

#### What it should do

1. Read a simple config (e.g., a vector or CSV) listing GEO IDs to use.
2. For each GEO dataset:
   - Download the data using `GEOquery`.
   - Extract:
     - Expression matrix
     - Sample metadata (e.g., tumor vs normal labels)
   - Standardize gene symbols (where possible).
3. Save processed objects to `data/processed/` (e.g., as `.rds`).

#### Inputs
- List of GEO IDs (could be hard-coded in the script or read from `docs/datasets.csv`).

#### Outputs
For each dataset, for example:
- `data/processed/<GSE_ID>_expr.rds` – gene × sample expression matrix.
- `data/processed/<GSE_ID>_pheno.rds` – sample metadata with at least:
  - Sample ID
  - Group label (e.g., tumor / normal)

#### Acceptance criteria
- ✅ For each GEO ID in the config, both `.rds` files exist in `data/processed/`.
- ✅ Each expression matrix:
  - Has >0 genes and >0 samples.
  - Column names match sample IDs in the phenotype table.
- ✅ Each phenotype table:
  - Has a clear grouping variable (e.g., a column named `group` with values like `tumor` / `normal`).

---

### 8.3 `02_run_deg_per_dataset.R`

#### Purpose (plain English)
For each dataset, find genes that are significantly up- or down-regulated in cancer vs normal.

#### What it should do

1. Loop over all processed datasets created by `01_download_and_prepare_data.R`.
2. For each dataset:
   - Create a design matrix (e.g., `~ group`).
   - Run `limma` (or `DESeq2`) to compute:
     - log2 fold change
     - p-values
     - adjusted p-values (FDR)
   - Add gene symbols as a column.
3. Save results per dataset and optionally combine into one big table.

#### Inputs
- `data/processed/<GSE_ID>_expr.rds`
- `data/processed/<GSE_ID>_pheno.rds`

#### Outputs
For each dataset:
- `results/deg/<GSE_ID>_deg_results.csv` (or `.rds`), with columns like:
  - `gene`, `log2FC`, `pvalue`, `padj`, `dataset`, `direction`

Optional combined file:
- `results/deg/all_datasets_deg_results.csv`

#### Acceptance criteria
- ✅ Script completes without errors.
- ✅ For each dataset:
  - Output file exists and has at least these columns:
    - `gene`, `log2FC`, `pvalue`, `padj`
  - There are some genes with `padj < 0.05`.
  - A quick check in R (or Excel) shows both positive and negative log2FCs.

---

### 8.4 `03_integrate_deg_across_datasets.R`

#### Purpose (plain English)
Combine DEG results from multiple datasets to identify robust genes that are consistently up- or down-regulated.

#### What it should do

1. Read all per-dataset DEG tables from `results/deg/`.
2. For each gene:
   - Count in how many datasets it is significant (e.g., `padj < 0.05`) and direction-consistent (all up or all down).
3. Define "robust up-DEGs" and "robust down-DEGs" based on a threshold, e.g.:
   - Up in ≥ N datasets.
   - Down in ≥ N datasets.
4. Optionally compute simple meta-analysis p-values (e.g., Fisher's method).

#### Inputs
- `results/deg/*_deg_results.csv` (from previous step).

#### Outputs
- `results/integration/robust_up_genes.csv`
- `results/integration/robust_down_genes.csv`

Optional:
- `results/integration/gene_support_summary.csv` – counts per gene of how many datasets support it.

#### Acceptance criteria
- ✅ Script completes without errors.
- ✅ `robust_up_genes.csv` and `robust_down_genes.csv` exist and:
  - Each has at least `gene` and `n_datasets_supporting` columns.
  - At least a modest number of up and/or down genes are present (not empty).
  - For a random subset of genes, you can cross-check that they are indeed DE in multiple per-dataset files.

---

### 8.5 `04_run_enrichment.R`

#### Purpose (plain English)
Run GO and pathway enrichment on the robust DEG lists to identify key biological processes and pathways.

#### What it should do

1. Load:
   - `results/integration/robust_up_genes.csv`
   - `results/integration/robust_down_genes.csv`
2. Use `clusterProfiler` (and `org.Hs.eg.db`, `ReactomePA`) to run:
   - GO Biological Process enrichment.
   - KEGG / Reactome pathway enrichment (if desired).
3. Produce:
   - Tables of enriched GO terms / pathways.
   - A few summary plots (barplots, dotplots).

#### Inputs
- Robust DEG lists from Step 8.4.

#### Outputs

**Tables**, e.g.:
- `results/enrichment/robust_up_go_bp.csv`
- `results/enrichment/robust_down_go_bp.csv`
- Optionally, pathway tables for KEGG/Reactome.

**Figures**, e.g.:
- `results/figures/enrichment_up_go_bp_barplot.png`
- `results/figures/enrichment_down_go_bp_barplot.png`

#### Acceptance criteria
- ✅ Script completes without errors.
- ✅ At least some GO terms/pathways have:
  - Adjusted p-values (`p.adjust`) below a chosen cutoff (e.g., 0.05).
  - Enriched terms make biological sense for breast cancer (e.g., cell cycle, DNA repair, proliferation, etc.).
- ✅ Plots open and are readable (axes labeled, term names visible).

---

### 8.6 `05_prepare_ppi_gene_lists.R`

#### Purpose (plain English)
Prepare clean gene lists from robust DEGs for upload to STRING / Cytoscape for PPI network analysis.

#### What it should do

1. Take robust up/down gene lists:
   - Optionally filter to top N genes (e.g., top 200 by |log2FC| or support).
2. Create text files with one gene symbol per line, in the format that STRING expects.
3. Optionally export a combined list and simple annotation table with:
   - Gene
   - Direction (up/down)
   - Support counts, etc.

#### Inputs
- `results/integration/robust_up_genes.csv`
- `results/integration/robust_down_genes.csv`

#### Outputs
- `results/ppi/genes_up_for_string.txt`
- `results/ppi/genes_down_for_string.txt`

Optional:
- `results/ppi/ppi_gene_annotations.csv`

#### Acceptance criteria
- ✅ Script completes without errors.
- ✅ Gene list files exist and:
  - Contain sensible human gene symbols.
  - Can be pasted/uploaded into the STRING web interface without errors.
- ✅ Optional check: counts of genes match your expectations (e.g., 100–500 genes, not 5 or 50,000).

---

### 8.7 `06_run_survival_analysis.R` (if using TCGA BRCA)

#### Purpose (plain English)
Test whether expression of selected genes is associated with patient survival.

#### What it should do

1. Load a TCGA BRCA expression + clinical dataset (or any dataset with survival info).
2. Choose a small set of genes of interest (e.g., hub genes from PPI network).
3. For each gene:
   - Split patients into "high" vs "low" expression groups (e.g., by median).
   - Fit survival curves (Kaplan–Meier) for the two groups.
   - Perform log-rank tests to see if survival is different.
4. Save:
   - Summary statistics (HR, p-values)
   - Kaplan–Meier plots.

#### Inputs
- `data/processed/tcga_brca_expr.rds` (or similar)
- `data/processed/tcga_brca_clinical.rds`
- A list of genes to test (hard-coded or read from a file).

#### Outputs
- `results/survival/survival_results.csv` – one row per gene with:
  - `gene`, `HR`, `log-rank p-value`, etc.

**Plots** such as:
- `results/figures/survival_<GENE>.png`

#### Acceptance criteria
- ✅ Script completes without errors.
- ✅ For at least some genes, KM plots show:
  - Two distinct curves (high vs low expression).
  - Reasonable number of patients in each group.
- ✅ At least a few genes have survival p-values that are not completely random (e.g., some near significance, even if not all < 0.05).
- ✅ The curves and stats are interpretable for the biology students (axes labeled, legends clear).

---

### 8.8 Summary

If all scripts above:

1. ✅ Run without crashing,
2. ✅ Produce the described output files, and
3. ✅ Look sensible on basic inspection,

then our pipeline is functionally complete, and the biology-focused members can concentrate on:

- Interpreting DEGs and pathways,
- Exploring PPI networks and hub genes,
- Interpreting survival plots,
- Connecting everything back to breast cancer biology and clinical questions.

---
