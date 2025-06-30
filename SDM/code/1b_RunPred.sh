#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=40G
#SBATCH -t 12:00:00
#SBATCH -J pred_run
#SBATCH -p normal_q
#SBATCH --account=usgs_rcs
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=chloe9mo@vt.edu

cd /fastscratch/chloe9mo/PNW_fishes

module load containers/singularity
echo "Modules loaded:"
module list
echo " "
echo "============================="
echo "Running from:"
pwd
echo " "
echo "============================="


echo "Running background pt creation step..."
echo "============================="
singularity exec --bind=/fastscratch/chloe9mo/PNW_fishes:/data /projects/arcsingularity/ood-rstudio141717-geospatial_4.1.1.sif Rscript /data/SDM/code/11_SDM_Predictions.R
echo "============================="
echo " "




