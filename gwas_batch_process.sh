#! /bin/bash

for i in $(seq -f "%02g" 1 22)
do
    sbatch -p short -n 6 -t 0-2:00 --mem=2G --job-name gwas-analysis -o %j.out -e %j.err --mail-type=ALL --mail-user=Peter_Shen@hms.harvard.edu --wrap="./gwas_single_chr.sh $i"
    sleep 1 # wait 1 second between each job submission
done
exit
