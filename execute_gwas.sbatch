#!/bin/bash

#SBATCH -p short   # queue name
#SBATCH -t 0-4:00       # hours:minutes runlimit after which job will be killed.
#SBATCH -n 4      # number of cores requested
#SBATCH --mem=24G # memory requested
#SBATCH -J gwas_analysis         # Job name
#SBATCH -o %j.out       # File to which standard out will be written
#SBATCH -e %j.err       # File to which standard err will be written
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Peter_Shen@hms.harvard.edu

# make dependencies
module load gcc/6.2.0
module load R/3.4.1

Rscript gwas.r -a "WTCCC/58C/Affx_gt_58C_Chiamo_22.tped.gz" -b "WTCCC/NBS/Affx_gt_NBS_Chiamo_22.tped.gz" -d "WTCCC/NBS/Affx_gt_NBS_Chiamo_22.tped.gz" -n 22 -s "WTCCC/58C/snps/snps_22"
