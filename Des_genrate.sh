#!/bin/sh
#Change the filename to your input filename
filename=example_smiles.csv
echo "Calculation of monomer properties"
python 01Des_monomer.py $filename
wait
echo "Calculation of Mordred descriptors"
python 02Des_Mordred.py $filename
wait
echo "Calculation of MD descriptors"
python 03Des_MD.py $filename
wait
echo "Descriptor files merge"
python 04Des_merger.py

echo "Successful! Automatic generation of physical descriptors"