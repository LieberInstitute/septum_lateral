#!/bin/bash
#$ -cwd
#$ -l mem_free=40G,h_vmem=40G,h_fsize=100G
#$ -N create_pseudobulk_data
#$ -o logs/01_create_pseudobulk_data
#$ -e logs/01_create_pseudobulk_data
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript 01_create_pseudobulk_data

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.2
## available from http://research.libd.org/sgejobs/
