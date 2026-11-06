# Does machine learning beat a transparent evidence rule? Prospective validation of human longevity gene prioritization

R code and data to reproduce the analyses of:

Cava C, Castiglioni I. *Does machine learning beat a transparent evidence rule? Prospective validation of human longevity gene prioritization* (manuscript in preparation, 2026).

The 213 candidate genes and the evidence rule come from our previous study, whose code is included here unchanged (scripts 00–21) and is also available at https://github.com/claudiacava/longevity-mirna-prioritization:

Cava C, Castiglioni I (2026) Integrative prioritization of human longevity genes and their miRNA regulators. *Netw Model Anal Health Inform Bioinform* 15:233. https://doi.org/10.1007/s13721-026-00908-y

## Design in brief

At a cut-off of 1 June 2017 (release of LongevityMap Build 3), every method ranks the 213 candidate genes using only information available before that date:

- **Rule**: at least two significant gene–study units, at least two population labels, and at least 60% of units significant (21 genes).
- **M1**: ridge logistic regression (primary) and balanced random forest trained on the same LongevityMap evidence by back-testing (pseudo-cut-offs 2013 and 2011).
- **M2**: two-step positive-unlabelled models trained on GO (2017), MSigDB 6.0, BioGRID 3.4.149 and GTEx v8, with GenAge genes as positives and the 213 candidates held out.
- **Baselines and combinations**: publications before the cut-off, number of GO terms, number of significant units; average rank of the rule with M1 or M2.

The 21 top-ranked genes of each method are then checked against evidence published after the cut-off: Open Genes (O1) and the GWAS Catalog, P ≤ 5 × 10⁻⁸ within 500 kb (O2).

## Requirements

- R ≥ 4.3 (analyses run with R 4.5.1)
- CRAN: readr, dplyr, tidyr, stringr, purrr, tibble, ggplot2, scales, patchwork, broom, jsonlite, openxlsx, glmnet (4.1-9), ranger (0.17.0), logistf (1.26.1)
- Bioconductor: AnnotationDbi, org.Hs.eg.db (3.21.0), GO.db, clusterProfiler, enrichplot, multiMiR

## Usage

Open R in this folder and run:

```r
source("run_all.R")       # previous study: candidate genes, rule, gene coordinates (scripts 00-21)
source("run_temporal.R")  # this study (scripts 30-36), including the label-permutation control
```

Large third-party files (full GWAS Catalog, NCBI gene2pubmed, BioGRID, GTEx, MSigDB) are not stored in the repository; scripts 30 and 33 download them and record source URL, date, size and MD5 checksum. The GWAS Catalog associations for the selected traits and the Open Genes snapshot used in the article (2 November 2026) are included, so the outcomes can be reproduced exactly.

## Scripts of this study

| Script | Description | Output |
|---|---|---|
| `R/30_download_temporal.R` | Dated snapshots of GWAS Catalog, gene2pubmed, Open Genes; publication dates from Europe PMC | `data/raw/temporal/`, `data/derived/` |
| `R/30b_fix_doi_map.R` | Corrects three Open Genes DOIs recorded with typographical errors | `data/derived/` |
| `R/31_temporal_validation.R` | Outcomes O1 and O2 after the cut-off; research attention; sensitivity analyses of the rule | `results/30_temporal/` |
| `R/32_ml_feasibility.R` | Training labels for back-testing (pseudo-cut-offs 2013 and 2011) | `results/32_ml_feasibility/` |
| `R/33_download_ml_features.R` | Feature sources for the M2 models | `data/raw/ml/` |
| `R/34_build_ml_features.R` | Feature matrices and gene universe | `data/derived/ml/`, `results/34_ml_features/` |
| `R/35_ml_models.R` | M1 and M2 models; scores of the candidate genes (with `ML_PERMUTE <- TRUE`: permuted-label control) | `results/35_ml_models/`, `results/35_ml_models_permuted/` |
| `R/36_rule_vs_ml.R` | Rule vs ML on the post-cut-off outcomes (the only script that reads outcomes and scores together) | `results/36_rule_vs_ml/`, `figures/` |

## Version history

- **v1.0** (November 2026): version used in the manuscript. Includes two corrections made before submission: three mistyped Open Genes DOIs are mapped to their PMIDs (script 30b), and the GWAS window is matched on official gene symbols (scripts 31–32).

## Licence

Code: MIT licence (see `LICENSE`). Data obtained from LongevityMap, Open Genes, the GWAS Catalog, NCBI, GenAge, Gene Ontology, MSigDB, BioGRID, GTEx and multiMiR remain subject to the terms of the original resources.

## Contact

Claudia Cava, University School for Advanced Studies IUSS Pavia
