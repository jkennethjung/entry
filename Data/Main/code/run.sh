#!/bin/bash
rm -rf ../output/
rm -rf ../temp/

mkdir ../output/
mkdir ../temp/
ln -s ../../../Raw/Infogroup/* ../temp/
ln -s ../../../Raw/CBP/* ../temp/

stata analysis.do 
