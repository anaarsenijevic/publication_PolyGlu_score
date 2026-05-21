# publication_PolyGlu_score

This repository contains scripts for the analysis of the Polyglutamylation (PolyGlu) score in cancer.

## TTLL/CCP polyglutamylation signature in cancer

This repository contains all scripts used to perform the bioinformatic analyses presented in:

**Arsenijevic et al.**  
*Transcriptional dysregulation of the microtubule polyglutamylation/deglutamylation enzyme network is a hallmark of human cancer*  
*(Title subject to final revision)*

---

## Overview

The scripts reproduce the main analyses of the study, from TCGA data processing to machine learning and external validation. The scripts are organised by figures.

### Main analyses
 
1. Expression profiling of TTLL and CCP (AGBL) genes across TCGA cancers  
2. Principal component analysis (PCA) of tumour vs normal samples  
3. Random forest classification of tumour vs normal samples  
4. Backward feature elimination and model reduction  
5. Segmented regression analysis to identify optimal gene panel size  
6. External validation in colorectal, kidney, and lung cancer datasets  
7. Comparison with cancer-specific classifiers  
8. Random gene panel and label permutation controls  
9. Scoring of cell lines using DepMap data  
10. Association analyses with TCGA clinical categories  
11. Survival analyses  
12. Scoring of stages of cell transformation  

---

## Requirements

R 4.4.2

**Main packages:**

- DESeq2 (1.46.0)  
- tximport (1.34.0)  
- tidyverse (2.0.0)  
- dplyr (1.1.4)  
- ggplot2 (3.5.2)  
- tibble (3.3.0)  
- tidyr (1.3.1)  
- purrr (1.0.4)  
- readr (2.1.5)  
- stringr (1.5.1)  
- forcats (1.0.0)  
- biomaRt (2.62.1)  
- AnnotationDbi (1.68.0)  
- org.Hs.eg.db (3.20.0)  
- GenomicFeatures (1.58.0)  
- SummarizedExperiment (1.36.0)  
- randomForest (4.7-1.2)  
- caret (7.0-1)  
- pROC (1.18.5)  
- segmented (2.2-1)  
- pheatmap (1.0.13)  
- patchwork (1.3.2)  
- ggpubr (0.6.2)  
- survival (3.7-0)  
- survminer (0.5.1)  
- TCGAbiolinks (2.34.1)
- MSigDBR (26.1.0)
- UCSCXenaTools (1.7.0)

---

## Notes

- Results may vary slightly due to random forest stochasticity  

---

## Citation

If you use this code, please cite:

Arsenijevic et al.  
*Transcriptional dysregulation of the microtubule polyglutamylation/deglutamylation enzyme network is a hallmark of human cancer*  
*(Title subject to final revision)*
