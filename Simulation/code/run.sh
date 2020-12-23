#!/usr/bin/sh
rm -rf ../output/
mkdir ../output/

matlab -nodisplay -nodesktop -r "run('./analysis.m')"
