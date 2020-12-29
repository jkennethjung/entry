#!/bin/sh
#SBATCH -t 24:00:00
#SBATCH -c 12 
module load MATLAB/2019a
rm -rf ../output/
mkdir ../output/

matlab -nodisplay -nodesktop -r "run('./analysis.m')"
