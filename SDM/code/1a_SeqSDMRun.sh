#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=80G
#SBATCH -t 24:00:00
#SBATCH -p normal_q
#SBATCH --account=usgs_rcs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=chloe9mo@vt.edu


echo "Date of job run:"
date
echo " "
cd /fastscratch/chloe9mo/PNW_fishes

module load containers/singularity
module list

## set variables
MotherFile="10_MaxEntSDM"
MotherFileR="${MotherFile}.R"
TodaysDate=$(date +%Y-%m-%d)
SandBox=SDM/SDMruns-$TodaysDate
mkdir /fastscratch/chloe9mo/PNW_fishes/$SandBox

echo "Beginning SDM runs..."

for index in $(seq 1 1 24); do
    cp /fastscratch/chloe9mo/PNW_fishes/SDM/code/$MotherFileR /fastscratch/chloe9mo/PNW_fishes/$SandBox/${MotherFile}_Species${index}.R
    sed -i "s/species.list\[1/species.list\[${index}/g" $SandBox/${MotherFile}_Species${index}.R
    echo "========================"
    echo "  >> Running Rscript:  ${SandBox}/${MotherFile}_Species${index}.R"
    singularity exec --bind=/fastscratch/chloe9mo/PNW_fishes:/data /projects/arcsingularity/ood-rstudio141717-geospatial_4.1.1.sif Rscript /data/$SandBox/${MotherFile}_Species${index}.R
    echo "========================"
    echo "  "
done

