library("here")
library("sessioninfo")
library("dplyr")
library("purrr")
library("spatialLIBD")
library("ComplexHeatmap")
library("circlize")



############## Load objects for gene_set_enrichment_plot_complex ##############

## Load gene_set_enrichment_objects.rda which includes the gene list (FDR < 0.05 and FDR < 0.01) and modeling results
load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "gene_set_enrichment_objects.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   modeling_results_broad
#   modeling_results
#   gene_list_FDR05
#   gene_list_FDR01

## Load gene_set_enrihcment_result_tables.rda. The FDR in the file names refers to the gene set filtering, not the fdr_cut used in gene_set_enrichment which was 0.1.
load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "gene_set_enrichment_1vsAll_result_tables.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   enrichTab_broad_glFDR01_ctFDR05
#   enrichTab_broad_glFDR05_ctFDR05
#   enrichTab_glFDR01_ctFDR05
#   enrichTab_glFDR05_ctFDR05

load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "gene_set_enrichment_1vs1_result_tables.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   prwiseTab_LS_ct_broad_glFDR01_ctFDR05
#   prwiseTab_LS_ct_broad_glFDR05_ctFDR05
#   prwiseTab_ct_LS_broad_glFDR01_ctFDR05
#   prwiseTab_ct_LS_broad_glFDR05_ctFDR05
#   prwiseTab_glFDR01_ctFDR05
#   prwiseTab_glFDR05_ctFDR05

lobstr::obj_size(enrichTab_glFDR05_ctFDR05)
# 4.10 kB
lobstr::obj_size(enrichTab_broad_glFDR05_ctFDR05)
# 5.88 kB
lobstr::obj_size(prwiseTab_glFDR05_ctFDR05)
# 23.25 kB
lobstr::obj_size(prwiseTab_broad_glFDR05_ctFDR05)
# 5.70 kB

###############################################################################



#### Function to plot with functions in gene_set_enrichment_plot_complex.R ####

source(
    here(
        "snRNAseq_mouse",
        "code",
        "04_clinical_set_enrichment",
        "gene_set_enrichment_plot_complex.R"
    )
)

use_gsepc <- function(modeling_results, model_type, gene_list, enrichTab, path_to_plot) {
    gene_enrichment_count <- get_gene_enrichment_count(model_results = modeling_results, model_type = model_type, bayes_anno = NULL)
    gene_list_count <- get_gene_list_count(gene_list)

    gse_plot <- gene_set_enrichment_plot_complex(
        enrichment = enrichTab,
        gene_count_col = log10(gene_list_count),
        gene_count_row = gene_enrichment_count,
        anno_title_col = "DE Genes\n(log10)",
        anno_title_row = "Cluster\nGenes"
    )

    pdf(path_to_plot, height = 4, width = 6)
    print(gse_plot)
    dev.off()
}

###############################################################################



############################# Plots for broad data ############################

#### 1vsAll ####
use_gsepc(
    modeling_results = modeling_results_broad,
    model_type = "enrichment",
    gene_list = gene_list_FDR05,
    enrichTab = enrichTab_broad_glFDR05_ctFDR05,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "GSEA-1vsAll_broad-glFDR05_ctFDR05.pdf"
    )
)

use_gsepc(
    modeling_results = modeling_results_broad,
    model_type = "enrichment",
    gene_list = gene_list_FDR01,
    enrichTab = enrichTab_broad_glFDR01_ctFDR05,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "GSEA-1vsAll_broad-glFDR01_ctFDR05.pdf"
    )
)

#### 1vs1 ####
modeling_results_broad$pairwise <- modeling_results_broad$'pairwise_LS-ct'
use_gsepc(
    modeling_results = modeling_results_broad,
    model_type = "pairwise",
    gene_list = gene_list_FDR05,
    enrichTab = prwiseTab_LS_ct_broad_glFDR05_ctFDR05,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "GSEA-1vs1_broad_LS_ct-glFDR05_ctFDR05.pdf"
    )
)

use_gsepc(
    modeling_results = modeling_results_broad,
    model_type = "pairwise",
    gene_list = gene_list_FDR01,
    enrichTab = prwiseTab_LS_ct_broad_glFDR01_ctFDR05,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "GSEA-1vs1_broad_LS_ct-glFDR01_ctFDR05.pdf"
    )
)

modeling_results_broad$pairwise <- modeling_results_broad$'pairwise_ct-LS'
use_gsepc(
    modeling_results = modeling_results_broad,
    model_type = "pairwise",
    gene_list = gene_list_FDR05,
    enrichTab =  prwiseTab_ct_LS_broad_glFDR05_ctFDR05,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "GSEA-1vs1_broad_ct_LS-glFDR05_ctFDR05.pdf"
    )
)

use_gsepc(
    modeling_results = modeling_results_broad,
    model_type = "pairwise",
    gene_list = gene_list_FDR01,
    enrichTab = prwiseTab_ct_LS_broad_glFDR01_ctFDR05,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "GSEA-1vs1_broad_ct_LS-glFDR01_ctFDR05.pdf"
    )
)

###############################################################################



########################## Plots for LS and Sept data #########################

## 1vsAll
use_gsepc(
    modeling_results = modeling_results,
    model_type = "enrichment",
    gene_list = gene_list_FDR05,
    enrichTab = enrichTab_glFDR05_ctFDR05,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "GSEA-1vsAll_8clusts-glFDR05_ctFDR05.pdf"
    )
)

use_gsepc(
    modeling_results = modeling_results,
    model_type = "enrichment",
    gene_list = gene_list_FDR01,
    enrichTab = enrichTab_glFDR01_ctFDR05,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "GSEA-1vsAll_8clusts-glFDR01_ctFDR05.pdf"
    )
)

## 1vs1
# use_gsepc(
#     modeling_results = modeling_results,
#     model_type = "pairwise",
#     gene_list = gene_list_FDR05,
#     enrichTab = prwiseTab_glFDR05_ctFDR05,
#     path_to_plot = here(
#         "snRNAseq_mouse",
#         "plots",
#         "04_clinical_set_enrichment",
#         "GSEA-1vs1_SeptLS-glFDR05_ctFDR05.pdf"
#     )
# )

# use_gsepc(
#     modeling_results = modeling_results,
#     model_type = "pairwise",
#     gene_list = gene_list_FDR01,
#     enrichTab = prwiseTab_glFDR01_ctFDR05,
#     path_to_plot = here(
#         "snRNAseq_mouse",
#         "plots",
#         "04_clinical_set_enrichment",
#         "GSEA-1vs1_SeptLS-glFDR01_ctFDR05.pdf"
#     )
# )

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

##############################################################################
