#! /bin/bash

#58C  BD  CAD  CD  HT  NBS  RA  T1D  T2D
#--mail-type=ALL --mail-user=Peter_Shen@hms.harvard.edu

declare -a arr=("BD" "CAD" "CD" "HT" "RA" "T1D")

for j in "${arr[@]}"
    do
    for i in $(seq -f "%02g" 1 22)
        do
            sbatch -p short -n 2 -t 0-4:00 --mem=16G --job-name $j-$i -o %j.out -e %j.err --wrap="./gwas_single_chr.sh $j $i"
    done
done
exit
