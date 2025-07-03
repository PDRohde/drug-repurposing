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
```markdown
```r
library(devtools)
devtools::install_github("psoerensen/gact")
devtools::install_github("psoerensen/qgg")
```

For details on how to use the `gact` and `qgg` packages, please refer to their respective GitHub documentation and associated scientific publications. You can find links to these resources and more at [https://pdrohde.github.io/](https://pdrohde.github.io/).

## Citation
If you use this code or refer to our results, please cite:
<div style="text-indent: -30px; padding-left: 30px;">
<p>Hjelholt AJ, Gholipourshahraki T, Bai Z, Shrestha M, Kjølby M, Sørensen P <b><span class="my-name">Rohde PD</span></b> (2025). <b>Leveraging Genetic Correlations to Prioritize Drug Groups for Repurposing in Type 2 Diabetes</b>, <em>medRxiv</em>,<a href="https://doi.org/10.1101/2025.06.13.25329590"> link</a> </p>
</div>

## Contact
For questions or feedback, please contact <font color="#E97132">palledr(at)hst.aau.dk</font> .

## License
This work is licensed under the **Creative Commons Attribution 4.0 International (CC BY 4.0)** license.  
You are free to share and adapt the material with appropriate credit.

[View license](https://creativecommons.org/licenses/by/4.0/)