#!/bin/sh
#SBATCH -t 24:00:00
module load MATLAB/2019a
rm -rf ../output/
rm -rf ../temp/

mkdir ../temp/
mkdir ../output/
ln -s ../../../Data/output/* ../temp/

matlab -nodisplay -nodesktop -r "run('./analysis.m')"
