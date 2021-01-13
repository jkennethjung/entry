#!/bin/bash
rm -rf ../output/
rm -rf ../temp/

mkdir ../output/
mkdir ../temp/
ln -s ../../../Raw/Infogroup/*.dta ../temp/
ln -s ../../../Raw/CBP/*.csv ../temp/
ln -s ../../../Raw/Geography/*.csv ../temp/
ln -s ../../Costs/output/*.dta ../temp/

stata analysis.do 
