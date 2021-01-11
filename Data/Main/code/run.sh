#!/bin/bash
rm -rf ../output/
rm -rf ../temp/

mkdir ../output/
mkdir ../temp/
ln -s ../../../Raw/Infogroup/*.dta ../temp/
ln -s ../../../Raw/CBP/*.csv ../temp/
ln -s ../../../Raw/Geography/*.csv ../temp/

stata analysis.do 
