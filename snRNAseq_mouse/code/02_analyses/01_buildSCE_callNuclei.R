### LC snRNA-seq analysis
### Build SCE from raw count matrices; perform nuclei calling
### qrsh -l bluejay,mf=52G,h_vmem=56G 
### Adapted from https://github.com/lmweber/locus-c/blob/main/code/analyses_sn/01_buildSCE_callNuclei.R
### Initiated MNT 08Mar2022

library(SingleCellExperiment)
library(DropletUtils)
library(rtracklayer)
library(BiocParallel)
library(jaffelab)
library(here)
library(sessioninfo)

here()
  #[1] "/dcs04/lieber/marmaypag/pilotLS_LIBD1070"

### Create raw SCE ====
# Basic sample info
sample.info <- read.table(here("snRNAseq_mouse", "raw_data","sample_info",
                               "sample_info_simplified.tsv"))

sample.info$path <- file.path(
  here("snRNAseq_mouse","processed_data", "cellranger"),
  sample.info$V1,
  "outs",
  "raw_feature_bc_matrix"
)
stopifnot(all(file.exists(sample.info$path)))

## Build basic SCE (will add more subject-level colData later, once obtained)
Sys.time()
    # [1] "2022-03-09 09:01:55 EST"
sce <- read10xCounts(
  samples = sample.info$path,
  sample.names = paste0(sample.info$V3,"_LS"),
  type = "sparse",
  col.names = TRUE
)
Sys.time()
    # [1] "2022-03-09 09:08:40 EST"


## Read in the gene information from the annotation GTF file
# (following Leo's method in https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/build_basic_sce.R)
gtf <-
  rtracklayer::import(
    "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-mm10-2020-A/genes/genes.gtf"
  )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[ , c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SCE object
rowRanges(sce) <- gtf[match_genes]

## Inspect object
sce
    # class: SingleCellExperiment 
    # dim: 32285 9521311 
    # metadata(1): Samples
    # assays(1): counts
    # rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ...
    # ENSMUSG00000095019 ENSMUSG00000095041
    # rowData names(6): source type ... gene_name gene_type
    # colnames(9521311): 1_AAACCCAAGAAACACT-1 1_AAACCCAAGAAACCAT-1 ...
    # 4_TTTGTTGTCTTTGATC-1 4_TTTGTTGTCTTTGCGC-1
    # colData names(2): Sample Barcode
    # reducedDimNames(0):
    #   mainExpName: NULL
    # altExpNames(0):

table(sce$Sample)
    # ms5478_ms5480_LS ms5479_ms5481_LS ms5483_ms5484_LS ms5485_ms5487_LS 
    # 2310153          2382021          2555720          2273417 

sce.ls <- sce

sample.idx <- splitit(sce.ls$Sample)
e.out.custom <- list()

Sys.time()
# 
for(i in names(sample.idx)){
  
  # Re-define the `lower=` param by identifying the 'second knee point'
  newknee <- barcodeRanks(counts(sce.ls[ ,sample.idx[[i]]]),
                          fit.bounds=c(10,1e3))
  cat(paste0("'Second knee point' for ", i," is: ",
             metadata(newknee)$knee,"\n"))

  # Set `lower = newknee + 100` (to capture the 'plateau' mode of low UMI totals)
  cat(paste0("Simulating empty drops for: ", i,"... \n"))
  set.seed(109)
  e.out.custom[[i]] <- emptyDrops(counts(sce.ls[ ,sample.idx[[i]]]), niters=20000,
                                  lower = metadata(newknee)$knee + 100,
                                  BPPARAM=BiocParallel::MulticoreParam(2))
  cat(paste0("\n\t...Simulations complete. \n\t", Sys.time(), "\n\n\n"))
}
Sys.time()
# 'Second knee point' for ms5478_ms5480_LS is: 430
# 'Second knee point' for ms5479_ms5481_LS is: 583
# 'Second knee point' for ms5483_ms5484_LS is: 560
# 'Second knee point' for ms5485_ms5487_LS is: 399


## For reference, plot the barcode rank plots
dir.create(here("snRNAseq_mouse","plots"))
for(i in names(sample.idx)){
  # Re-define the `lower=` param by identifying the 'second knee point'
  newknee <- barcodeRanks(counts(sce.ls[ ,sample.idx[[i]]]),
                          fit.bounds=c(10,1e3))
  cat(paste0("'Second knee point' for ", i," is: ",
             metadata(newknee)$knee,"\n"))
  
  # Plot the barcode rank plot with +100 and +300
  png(here("snRNAseq_mouse","plots",
           paste0("barcodeRankPlot_",i,"_dropletUtils-w-fitbounds_plus100-300UMIs.png")))
  plot(newknee$rank, newknee$total,
       log="xy", xlab="Barcode Rank", ylab="Total UMIs",
       main=paste0("Barcode rank plot for: ",i,"\n( `fit.bounds=c(10,1e3)` )"),
       cex.axis=0.8, cex.lab=1.2, las=1)
    o <- order(newknee$rank)
    lines(newknee$rank[o], newknee$fitted[o], col="red")
    # k_2 from above
    abline(h=metadata(newknee)$knee, col="darkblue", lty=2)
    # k_2 + 100:
    abline(h=metadata(newknee)$knee + 100, col="dodgerblue", lty=2)
    # k_2 + 300:
    abline(h=metadata(newknee)$knee + 300, col="darkgreen", lty=2)
    legend("bottomleft", lty=2, col=c("darkblue", "dodgerblue", "darkgreen"),
           legend=c("knee_2", "knee_2+100", "knee_2+300"))
    dev.off()
    
}


## emptyDrops() results ===
for(i in 1:length(e.out.custom)){
  print(names(e.out.custom)[[i]])
  print(table(Signif = e.out.custom[[i]]$FDR <= 0.001,
              Limited = e.out.custom[[i]]$Limited))
}

    #[1] "ms5478_ms5480_LS"
    #       Limited
    # Signif  FALSE  TRUE
    #   FALSE 36701     0
    #   TRUE     27  6676

    # [1] "ms5479_ms5481_LS"
    #       Limited
    # Signif  FALSE  TRUE
    #   FALSE 33846     0
    #   TRUE     21  6855

    # [1] "ms5483_ms5484_LS"
    #       Limited
    # Signif  FALSE  TRUE
    #   FALSE 47113     0
    #   TRUE    742 17098   - even with a sample-defined `lower=` threshold,
    #                         this is way higher than should be... (2F_C_LS)

    # [1] "ms5485_ms5487_LS"
    #       Limited
    # Signif  FALSE  TRUE
    #   FALSE 32454     0
    #   TRUE     30  7786


    # and btw these sample names are not exactly the same order:
    table(ss(colnames(sce.ls),"_",1), sce.ls$Sample)
        #   ms5478_ms5480_LS ms5479_ms5481_LS ms5483_ms5484_LS ms5485_ms5487_LS
        # 1          2310153                0                0                0
        # 2                0                0          2555720                0
        # 3                0          2382021                0                0
        # 4                0                0                0          2273417
    sample.info
        #        V1             V2            V3      V4
        # 1 1M_C_LS Lateral_Septum ms5478_ms5480 1M-C-LS
        # 2 2F_C_LS Lateral_Septum ms5483_ms5484 2F-C-LS
        # 3 3M_C_LS Lateral_Septum ms5479_ms5481 3M-C-LS
        # 4 4F_C_LS Lateral_Septum ms5485_ms5487 4F-C-LS

    
## Save this with the raw SCE for interactive downstream analyses ===
dir.create(here("snRNAseq_mouse","processed_data","SCE"))
README.custom <- "Object 'e.out.custom' is the output of emptyDrops() with lower= set to the quantified
  'second knee point' (+100) to better model empty/debris-containing droplets"
save(sce.ls,
     e.out.custom, README.custom,
     file=here("snRNAseq_mouse","processed_data","SCE", "sce_raw_LS.rda"))



## Stricter approach for sample 'ms5483_ms5484_LS'? Make refined barcodeRanks stats:
bcranks.2F <- barcodeRanks(counts(sce.ls[ ,sample.idx[["ms5483_ms5484_LS"]]]),
                           fit.bounds=c(1e3,1e5), lower=560)
                                                # setting lower to the above-ID'd '2nd knee'
                                                # bc otherwise is ID'ing ~the '2nd inflxn point'
metadata(bcranks.2F)
    # $knee
    # [1] 10524
    # 
    # $inflection
    # [1] 6491

table(bcranks.2F$total >= metadata(bcranks.2F)$inflection)  # 3983

# What about restricting to this and having a high scores emptyDrops?
table(rownames(bcranks.2F) == rownames(e.out.custom[["ms5483_ms5484_LS"]])) # all good

table(bcranks.2F$total >= metadata(bcranks.2F)$inflection &
        e.out.custom[["ms5483_ms5484_LS"]]$FDR <= 0.001)
    #  FALSE    TRUE 
    #2551738    3982   (about the same)

# For reference:

plot(bcranks.2F$rank, bcranks.2F$total,
     #col=ifelse(e.out.custom[["ms5483_ms5484_LS"]]$FDR <= 0.001,"darkblue","red"),
     log="xy", xlab="Barcode Rank", ylab="Total UMIs",
     main=paste0("Barcode rank plot for: 2F_ms5483_ms5484\n( `fit.bounds=c(10,1e3), lower=560` )"),
     cex.axis=0.8, cex.lab=1.2, las=1)
o <- order(bcranks.2F$rank)
lines(bcranks.2F$rank[o], bcranks.2F$fitted[o], col="red")
# k_1 from barcodeRanks
abline(h=metadata(bcranks.2F)$knee, col="darkblue", lty=2)
# k_2 previously ID'd, above:
abline(h=560, col="dodgerblue", lty=2)
abline(h=860, col="dodgerblue", lty=2) # k_2+100

# Inflection point to threshold on, with this lower quality sample:
abline(h=metadata(bcranks.2F)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("darkblue", "dodgerblue", "forestgreen"),
       legend=c("knee_1", "knee_2", "inflection"))








## Filter based on this revised nuclei-calling approach and save into another .rda ===
# BCs.keep <- unlist(
#   lapply(e.out.custom, function(x){ rownames(x)[which(x$FDR <= 0.001)] }
#   ))
# 
# sce.ls <- sce.ls[ ,BCs.keep]
# 
# save(sce.ls, file=here("snRNAseq_mouse","processed_data","SCE", "sce_working_LS.rda"))


## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
#
proc.time()
options(width = 120)
session_info()
