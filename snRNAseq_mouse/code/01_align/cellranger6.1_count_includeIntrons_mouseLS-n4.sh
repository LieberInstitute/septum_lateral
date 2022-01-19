#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -pe local 4
#$ -N mouseLS_n4
#$ -o logs/mouseLS_n4.$TASK_ID.txt
#$ -e logs/mouseLS_n4.$TASK_ID.txt
#$ -m ea
#$ -M nnguye44@jhmi.edu
#$ -t 1-4
#$ -tc 2

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## List current modules for reproducibility
module list

## load CellRanger
module load cellranger/6.1.1

## Locate file
FASTQDIR=/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/raw_data/FASTQ
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' ${FASTQDIR}/../sample_info/sample_info_simplified.tsv | awk "NR==${SGE_TASK_ID}")

# Fixed FASTQ file sample prefix ("-" instead of "_")
FIXED=$(awk 'BEGIN {FS="\t"} {print $4}' ${FASTQDIR}/../sample_info/sample_info_simplified.tsv | awk "NR==${SGE_TASK_ID}")

echo "Processing sample ${SAMPLE}"
date

## Run CellRanger
cellranger count --id=${SAMPLE} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-mm10-2020-A \
    --fastqs=${FASTQDIR}/${SAMPLE} \
    --sample=${FIXED} \
    --jobmode=local \
    --localcores=4 \
    --localmem=20 \
    --include-introns \
    --r1-length=28

## Move output
echo "Moving data to new location"
date
mv ${SAMPLE} /dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/cellranger/

echo "**** Job ends ****"
date

