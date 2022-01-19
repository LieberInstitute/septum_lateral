# From `/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/raw_data`

```R
library(readxl)
sampleInfo <- as.data.frame(read_excel("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/raw_data/sample_info/snRNAseq_LateralSeptum_Master_Spreadsheet_20211118.xlsx",
			    sheet=1, col_names=T))

# Simplify - only need first few columns (ignore example 'DLPFC Mid' sample)
sampleInfo <- sampleInfo[2:5, 1:3]

# Clean up
colnames(sampleInfo) <- gsub(" ","_",colnames(sampleInfo))
colnames(sampleInfo) <-	gsub("#","Num",colnames(sampleInfo))
sampleInfo$Tissue <- gsub(" ","_",sampleInfo$Tissue)
sampleInfo$Animal_ID <- gsub("\\+","_", sampleInfo$Animal_ID)

# Add changed FASTQ file prefix... (only the files; the encasing folders are fine)
sampleInfo$filePrefix <- gsub("_","-", sampleInfo$Sample_Num)


write.table(sampleInfo, quote=F, col.names=F, row.names=F, sep="\t",
	    file="sample_info/sample_info_simplified.tsv")


```

# Organizing FASTQs

```bash
## 1M_C_LS
mkdir FASTQ/1M_C_LS
ln -s /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-12-19_SPag121621/1M_C_LS*/* FASTQ/1M_C_LS/

## 2F_C_LS
mkdir FASTQ/2F_C_LS
ln -s /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-12-19_SPag121621/2F_C_LS*/* FASTQ/2F_C_LS/

## 3M_C_LS
mkdir FASTQ/3M_C_LS
ln -s /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-12-19_SPag121621/3M_C_LS*/* FASTQ/3M_C_LS/

## 4F_C_LS
mkdir FASTQ/4F_C_LS
ln -s /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2021-12-19_SPag121621/4F_C_LS*/* FASTQ/4F_C_LS/


```
