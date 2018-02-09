#! /bin/bash

for i in {1..2}
do
    sbatch -p short -n 6 -t 0-2:00 --mem=2G --job-name gwas-analysis -o %j.out -e %j.err --mail-type=ALL --mail-user=Peter_Shen@hms.harvard.edu --wrap="sh gwas_single_chr $i"
    sleep 1 # wait 1 second between each job submission
done
exit
