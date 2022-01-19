### MNT 29Nov2021 =======
  # Check Read 1 files for discrepant read lengths, as seen before in Tran-Maynard 2021 project:
  # Should be exactly 28bp = 16 [BC] + 12 [UMI]

library(ShortRead)
library(jaffelab)
library(sessioninfo)

FASTQ.dir <- "/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/raw_data/FASTQ/"

### Read in abridged sample info (MNT generated for processing through CR)
sampleInfo <- read.table("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/raw_data/sample_info/sample_info_simplified.tsv",
                            sep="\t", header=F)$V1

# Reference folder/sample name (not the edited FASTQ names)
R1files <- data.frame(
  sampleName = unlist(sapply(sampleInfo, function(x){
    rep(x,length(list.files(paste0(FASTQ.dir, x), pattern="R1")))}), use.names=F),
  
  R1 = unlist(sapply(sampleInfo, function(x){list.files(paste0(FASTQ.dir,x),
                                                   pattern="R1")}), use.names=F)
)
dim(R1files)  # 

for(i in 1:nrow(R1files)){
  cat(paste0("Checking R1 length distribution for: ", R1files[i,2], "\n"))
  temp.R1s <- readFastq(paste0(FASTQ.dir, R1files[i,1], "/", R1files[i,2]),
                        withIds=F)
  print(head(sread(temp.R1s), n=4))
  print(table(width(sread(temp.R1s))))
  rm(temp.R1s)
  cat("\n\n")
}


## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
