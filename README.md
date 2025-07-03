# Drug Repurposing in Type 2 Diabetes using Multi-Trait Bayesian Gene Set Analysis

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

This repository contains code and data used in the study:  
**Leveraging Genetic Correlations to Prioritize Drug Groups for Repurposing in Type 2 Diabetes**

## Study Overview
Type 2 diabetes (T2D) is a complex polygenic disorder with limited success in therapeutic translation despite numerous genome-wide association studies (GWAS). In this study, we developed and applied a **Bayesian Linear Regression (BLR) gene set framework** to integrate GWAS summary statistics with curated drug-gene interaction data from the **Drug-Gene Interaction Database (DGIdb)**.

We prioritized druggable gene sets, evaluated at the **ATC 4th level**, by calculating **posterior inclusion probabilities (PIP)** for each set. Our findings:
- Validated the model by identifying strong genetic associations for known antidiabetic agents
- Highlighted additional candidate drug classes with genetic relevance to T2D, including:
  - Carboxamide derivatives
  - Fibrates
  - Uric acid inhibitors
  - Immunomodulatory and antineoplastic agents  
- Identified **key T2D-associated genes** such as *PPARG*, *KCNQ1*, *TNF*, and *GCK*
- Suggested **bezafibrate**, a PPAR pan-agonist, as a potential repurposing candidate

## Requirements
To run our analyses, you need to install the following R packages from GitHub. You can do so using the devtools package:
```r
library(devtools)
devtools::install_github("psoerensen/gact")
devtools::install_github("psoerensen/qgg")
```

For details on how to use the `gact` and `qgg` packages, please refer to their respective GitHub documentation and associated scientific publications. You can find links to these resources and more at [https://pdrohde.github.io/](https://pdrohde.github.io/).

## gact Workflow Scripts

This repository contains seven R scripts that together form a complete workflow for building and analyzing a GWAS database using the `gact` framework. The scripts are organized to facilitate reproducible research, from data ingestion to gene- and pathway-level interpretation.

### Overview of Scripts

#### 1. Install and Setup
Installs the `gact` and `qgg` R packages, downloads a versioned gact database of GWAS summary statistics and annotations, and sets up the infrastructure for downstream analyses. Includes commands to explore the content and structure of the database.

#### 2. Prepare Genotype Summary Data
Generates genotype summary objects (`Glist`) using 1000 Genomes (1000G) data. This includes full and ancestry-specific datasets (e.g., EUR). Variants are filtered and matched to those in the gact database. Filtered `Glist` objects are saved for later use. Similar procedures can be used for other ancestries.

#### 3. Compute Sparse LD Matrices
Computes sparse linkage disequilibrium (LD) matrices using high-quality variants from the 1000G reference panel. Applies filtering criteria (e.g., MAF, missingness, indels, HWE) and saves updated `Glist` and an annotated marker table with LD scores. Fast computation is recommended using OpenBLAS, MKL, or similar libraries.

#### 4. Ingest GWAS Summary Statistics
Harmonizes and ingests GWAS summary statistics for nine traits into the gact database. Performs allele matching and quality control, and annotates each dataset with relevant metadata (e.g., trait type, ancestry, sample size, publication reference). Uses `updateStatDB()` for integration.

#### 5. Gene-Level Analysis (VEGAS)
Performs gene-level association analysis using the VEGAS method for all nine traits. For each study, GWAS summary statistics are aligned with the reference panel and tested for enrichment across genes (±40kb upstream / 10kb downstream). Results are saved per trait.

#### 6. Multi-Trait Gene-Set Analysis (MAGMA)
Performs multi-trait Bayesian gene-set enrichment analysis using the `magma()` function. Gene-level Z-scores (from VEGAS) are combined across traits and analyzed using ATC level 4 drug gene sets from the gact database to identify enriched pathways.

#### 7. Disease-Gene Enrichment (Hypergeometric Test)
Runs hypergeometric tests to identify ATC level 4 drug classes that are significantly enriched for diabetes-related gene sets. Tests are performed using multiple sources of gene-disease associations (curated knowledge, text mining, experimental data).

## Citation
If you use this code or refer to our results, please cite:
<div style="text-indent: -30px; padding-left: 30px;">
<p>Hjelholt AJ, Gholipourshahraki T, Bai Z, Shrestha M, Kjølby M, Sørensen P, <b><span class="my-name">Rohde PD</span></b> (2025). <b>Leveraging Genetic Correlations to Prioritize Drug Groups for Repurposing in Type 2 Diabetes</b>, <em>medRxiv</em>,<a href="https://doi.org/10.1101/2025.06.13.25329590"> link</a> </p>
</div>

## Contact
For questions or feedback, please contact <font color="#E97132">palledr(at)hst.aau.dk</font> .

## License
This work is licensed under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.  
You are free to share and adapt the material with appropriate credit.

[View license](https://creativecommons.org/licenses/by/4.0/)