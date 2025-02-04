septum_lateral
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![DOI](https://zenodo.org/badge/445640720.svg)](https://zenodo.org/badge/latestdoi/445640720)

## Overview

Welcome to the septum_lateral project! These experiments generated
several interactive websites for you to browse and download.

In our previous work,
[www.nature.com/articles/s41386-022-01487-y](https://www.nature.com/articles/s41386-022-01487-y),
we demonstrated that expression of tropomyosin kinase receptor B (TrkB)
in LS neurons is required for social novelty recognition. To better
understand how TrkB signaling may control social behavior, we locally
knocked down TrkB in the LS and used bulk RNA-sequencing to identify
changes in gene expression between 4 mice with TrkB knockdown in the LS
(cre) and 4 control mice with viral expression of GFP (GFP).
Additionally, we generated molecular profiles for LS cell types using
single nucleus RNA-sequencing (snRNA-seq) on N = 4 samples (generated
from n=4 males, n=4 females, each sample contained pooled LS dissections
of 2 mice of the same sex). Once these datasets were generated, our
final goal was to intersect them to identify whether DEGs identified
from the TrkB knockdown bulk RNA-seq data are enriched in the cell types
identified using our snRNA-seq data.

For more detailed information about study design and experimental
results, please refer to our manuscript
[www.biorxiv.org/content/10.1101/2023.06.29.547069v2](https://www.biorxiv.org/content/10.1101/2023.06.29.547069v2).
This work was performed by Keri Martinowich, Stephanie Cerceo Page, and
Leonardo Collado-Torres teams at the Lieber Institute for Brain
Development.

This project involves the GitHub repository
[github.com/LieberInstitute/septum_lateral](https://github.com/LieberInstitute/septum_lateral)

## Study Design

<img src="http://research.libd.org/septum_lateral/img/study_overview.png" width="1000px" align="left" />

**Study design to identify molecular changes induced by TrkB knockdown
in newly identified LS cell-types using snRNA-seq**. (**A**) Schematic
of viral strategy using cre-mediated recombination to locally knockdown
TrkB expression in the LS (n=4 TrkB KD, n=4 TrkB Control). (**B**)
Schematic of experimental design for snRNA-seq of mouse LS tissue.
Tissues from 2 mice of the same sex were pooled together for each
individual sample for a total of N=4 samples (generated from n=4 male,4
female mice). (**C**) Volcano plot of differentially expressed genes in
bulk RNA-seq of LS tissue samples comparing control versus local TrkB
knockdown. (**D**) Uniform manifold approximation and projection (UMAP)
of identified cell types, with nuclei counts per clusters. (**E**)
Enrichment analysis of all, positive, and negative DEGs from the TrkB
knockdown dataset across broad cellular clusters. (**F**) Heatmap of DE
genes from the TrkB knockdown dataset enriched in lateral septal and
microglial broad clusters across broad clusters. DE genes from the TrkB
knockdown dataset are divided into genes unique to the LS, plasticity
genes, neurodevelopmental, and microglia genes.

## Interactive websites:

All of these interactive websites are powered by open source software,
namely:

- 👀 [`iSEE`](https://doi.org/10.12688%2Ff1000research.14966.1)

We provide the following interactive websites, organized by dataset with
software labeled by emojis:

- snRNA-seq (n = 4)
  - 👀
    [snRNAseq_lateral_septum](https://libd.shinyapps.io/snRNAseq_lateral_septum/):
    main website for the snRNA-seq data from this study.
- bulk RNA-seq (n = 8)
  - 👀
    [bulkseq_lateral_septum](https://libd.shinyapps.io/bulkseq_lateral_septum/):
    similar to `snRNAseq_lateral_septum`, but for the bulk RNA-seq data.

## Contact

We value public questions, as they allow other users to learn from the
answers. If you have any questions, please ask them at
[LieberInstitute/septum_lateral/issues](https://github.com/LieberInstitute/septum_lateral/issues)
and refrain from emailing us. Thank you again for your interest in our
work!

## Citing our work

Please cite this [manuscript](https://doi.org/10.1101/TODO) if you use
data from this project.

> TODO

Below is the citation in [`BibTeX`](http://www.bibtex.org/) format.

    @article {TODO
    }

## Data Access

We highly value open data sharing and believe that doing so accelerates
science.

### Processed Data

TODO (by Leo)

### Raw data

The source data is publicly available from the Globus endpoint
`jhpce#septum_lateral`, which is also listed at
<http://research.libd.org/globus>.

## Internal

- JHPCE path: `/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse`
- Files at `snRNAseq_mouse` are mostly organized following
  <https://github.com/LieberInstitute/template_project>.
- Slack channel:
  [`libd_lateralseptum_snrna-seq`](https://jhu-genomics.slack.com/archives/C01KSLPCGAC)
