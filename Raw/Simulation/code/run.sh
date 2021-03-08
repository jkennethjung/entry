#!/bin/sh
#SBATCH -t 1-00:00:00
#SBATCH -c 24 
module load MATLAB/2019a
rm -rf ../output/
rm -rf ../temp/

mkdir ../temp/
mkdir ../output/
ln -s ../../../Data/Simulation/output/* ../temp/

matlab -nodisplay -nodesktop -r "run('./analysis.m')"
