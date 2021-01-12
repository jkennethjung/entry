#!/bin/bash
rm -rf ../output/
rm -rf ../temp/

mkdir ../output/
mkdir ../temp/
ln -s ../../../Raw/CBP/*.csv ../temp/
ln -s ../../../Raw/Geography/*.csv ../temp/
ln -s ../../../Raw/QCEW/*.csv ../temp/
ln -s ../../../Raw/HUD/*.csv ../temp/

stata analysis.do 
