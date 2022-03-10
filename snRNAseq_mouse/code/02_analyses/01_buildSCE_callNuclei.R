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
    # [1] "2022-03-10 07:53:33 EST"
sce <- read10xCounts(
  samples = sample.info$path,
  sample.names = ss(sample.info$V1, "_", 1),
  type = "sparse",
  col.names = TRUE
)
Sys.time()
    # [1] "2022-03-10 07:56:43 EST"


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
    #     1M      2F      3M      4F 
    #2310153 2555720 2382021 2273417 

sce.ls <- sce
sample.idx <- splitit(sce.ls$Sample)

## First use barcodeRanks to ID the '2nd knee' and make some plots to diagnose ===
#dir.create(here("snRNAseq_mouse","plots"))
bcRanks.ls <- list()

Sys.time()
    # [1] "2022-03-10 08:09:53 EST"
for(i in names(sample.idx)){
  
  # Re-define the `lower=` param by identifying the 'second knee point'
  bcRanks.ls[[i]] <- barcodeRanks(counts(sce.ls[ ,sample.idx[[i]]]),
                          fit.bounds=c(10,1e3))
  cat(paste0("'Second knee point' for ", i," is: ",
             metadata(bcRanks.ls[[i]])$knee,"\n"))
  
  # Plot the barcode rank plot with +100 and +300
  png(here("snRNAseq_mouse","plots",
           paste0("barcodeRankPlot_",i,"_w-fitbounds_k2plus100-300UMIs.png")))
  plot(bcRanks.ls[[i]]$rank, bcRanks.ls[[i]]$total,
       log="xy", xlab="Barcode Rank", ylab="Total UMIs",
       main=paste0("Barcode rank plot for: ",i,"\n( `fit.bounds=c(10,1e3)` )"),
       cex.axis=0.8, cex.lab=1.2, las=1)
    o <- order(bcRanks.ls[[i]]$rank)
    lines(bcRanks.ls[[i]]$rank[o], bcRanks.ls[[i]]$fitted[o], col="red")
    # k_2 from above
    abline(h=metadata(bcRanks.ls[[i]])$knee, col="darkblue", lty=2)
    # k_2 + 100:
    abline(h=metadata(bcRanks.ls[[i]])$knee + 100, col="dodgerblue", lty=2)
    # k_2 + 300:
    abline(h=metadata(bcRanks.ls[[i]])$knee + 300, col="darkgreen", lty=2)
    legend("bottomleft", lty=2, col=c("darkblue", "dodgerblue", "darkgreen"),
           legend=c("knee_2", "knee_2+100", "knee_2+300"))
    dev.off()
    
}
Sys.time()
    # 'Second knee point' for 1M is: 430
    # 'Second knee point' for 2F is: 560
    # 'Second knee point' for 3M is: 583
    # 'Second knee point' for 4F is: 399

    # Diagnosis: +300 UMIs to the '2nd knee point' will be better than +100 that sufficed
    #            for the LC project.  However, it would be better to calculate the 2nd deriv
    #            minimum to adaptively identify where this '2nd plateau' starts




## NOW run emptyDrops ===
e.out.custom <- list()

Sys.time()
    # [1] "2022-03-10 08:46:33 EST"
for(i in names(sample.idx)){
  
  cat(paste0("'Second knee point' for ", i," is: ",
             metadata(bcRanks.ls[[i]])$knee,"\n"))
  cat(paste0("\tSetting `lower` for `empyDrops()` to = ",
             metadata(bcRanks.ls[[i]])$knee+300,"\n\n"))
  
  # Set `lower = bcRanks.ls[[i]] + 300` (to capture the 'plateau' mode of low UMI totals)
  cat(paste0("\tSimulating empty drops for: ", i,"... \n"))
  set.seed(109)
  e.out.custom[[i]] <- emptyDrops(counts(sce.ls[ ,sample.idx[[i]]]), niters=20000,
                                  lower = metadata(bcRanks.ls[[i]])$knee + 300
                                  )
  cat(paste0("\n\t...Simulations complete. \n\t", Sys.time(), "\n\n\n"))

  }
    # 


## emptyDrops() results ===
for(i in 1:length(e.out.custom)){
  print(names(e.out.custom)[[i]])
  print(table(Signif = e.out.custom[[i]]$FDR <= 0.001,
              Limited = e.out.custom[[i]]$Limited))
}

      # [1] "1M"
      #       Limited
      # Signif  FALSE TRUE
      #   FALSE  7918    0
      #   TRUE     67 6661

      # [1] "2F"
      #       Limited
      # Signif  FALSE  TRUE
      #   FALSE 23312     0
      #   TRUE   1033 18588   - this is still too high...

      # [1] "3M"
      #       Limited
      # Signif  FALSE TRUE
      #   FALSE  6172    0
      #   TRUE     54 6907

      # [1] "4F"
      #       Limited
      # Signif  FALSE TRUE
      #   FALSE  3128    0
      #   TRUE     70 7786

# For reference:
sample.info
    #        V1             V2            V3      V4
    # 1 1M_C_LS Lateral_Septum ms5478_ms5480 1M-C-LS
    # 2 2F_C_LS Lateral_Septum ms5483_ms5484 2F-C-LS
    # 3 3M_C_LS Lateral_Septum ms5479_ms5481 3M-C-LS
    # 4 4F_C_LS Lateral_Septum ms5485_ms5487 4F-C-LS

    
## Save this with the raw SCE for interactive downstream analyses ===
dir.create(here("snRNAseq_mouse","processed_data","SCE"))
README.custom <- "Object 'e.out.custom' is the output of emptyDrops() with lower= set to the quantified
  'second knee point' (+300) to better model empty/debris-containing droplets"
save(sce.ls,
     e.out.custom, README.custom,
     file=here("snRNAseq_mouse","processed_data","SCE", "sce_raw_LS.rda"))


## Btw without a 'lower=' param for barcodeRanks, it'll most likely pick the
 #    '2nd' inflection point
bcRanks.2F <- barcodeRanks(counts(sce.ls[ ,sample.idx[["2F"]]]),
                           fit.bounds=c(1e3,1e5),
                           lower=bcRanks.ls[["2F"]] + 300)

metadata(bcRanks.2F)
    #      knee inflection 
    #     10524       6491

# Re-plot this with that inflection point
i <- "2F"
png(here("snRNAseq_mouse","plots",
         paste0("barcodeRankPlot_2F_w-fitbounds_k2plus100-300UMIs.png")))
plot(bcRanks.ls[[i]]$rank, bcRanks.ls[[i]]$total,
     log="xy", xlab="Barcode Rank", ylab="Total UMIs",
     main=paste0("Barcode rank plot for: ",i,"\n( `fit.bounds=c(10,1e3)` )"),
     cex.axis=0.8, cex.lab=1.2, las=1)
o <- order(bcRanks.ls[[i]]$rank)
lines(bcRanks.ls[[i]]$rank[o], bcRanks.ls[[i]]$fitted[o], col="red")
# k_2 from above
abline(h=metadata(bcRanks.ls[[i]])$knee, col="darkblue", lty=2)
# k_2 + 100:
abline(h=metadata(bcRanks.ls[[i]])$knee + 100, col="dodgerblue", lty=2)
# k_2 + 300:
abline(h=metadata(bcRanks.ls[[i]])$knee + 300, col="darkgreen", lty=2)
# Re-computed inflection:
abline(h=metadata(bcRanks.2F)$inflection, col="darkred", lty=2)
legend("bottomleft", lty=2, col=c("darkblue", "dodgerblue", "darkgreen", "darkred"),
       legend=c("knee_2", "knee_2+100", "knee_2+300", "inflection"))
dev.off()


# What about restricting to this and having a high scores emptyDrops?
table(rownames(bcRanks.2F) == rownames(e.out.custom[["2F"]])) # all good

table(bcRanks.2F$total >= metadata(bcRanks.2F)$inflection &
        e.out.custom[["2F"]]$FDR <= 0.001)
    #  FALSE    TRUE
    #2551738    3982




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
