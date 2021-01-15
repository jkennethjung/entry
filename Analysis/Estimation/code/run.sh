#!/bin/sh
#SBATCH -t 24:00:00
#SBATCH -c 36 
module load MATLAB/2019a
rm -rf ../output/
rm -rf ../temp/

mkdir ../temp/
mkdir ../output/
ln -s ../../../Data/Simulation/output/* ../temp/

matlab -nodisplay -nodesktop -r "run('./analysis.m')"
