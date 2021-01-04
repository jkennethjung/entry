#!/bin/bash
rm -rf ../output/
rm -rf ../temp/

mkdir ../output/
mkdir ../temp/
ln -s ../../Raw/Simulation/output/* ../temp/

stata analysis.do 
