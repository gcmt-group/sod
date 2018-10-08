#!/bin/bash
# This script runs the spbe code (simple pair-based extrapolation) programe of the SOD package, starting from x=1.
# It requires the following input files: EQMATRIX, OUTSOD, OUTSOD1, OUTSOD2, ENERGIES0, ENERGIES1, ENERGIES2
# The following commands might be useful to collect the input files - comment out if you already have them.
cp ../../EQMATRIX ./
cp ../OUTSOD ./OUTSOD_original
invertOUTSOD
mv OUTSOD_inverted OUTSOD
cp ../../n11/OUTSOD ./OUTSOD_original
invertOUTSOD
mv OUTSOD_inverted OUTSOD1 
cp ../../n10/OUTSOD ./OUTSOD_original
invertOUTSOD
mv OUTSOD_inverted OUTSOD2 
cp ../../n12/ENERGIES ./ENERGIES0
cp ../../n11/ENERGIES ./ENERGIES1
cp ../../n10/ENERGIES ./ENERGIES2
rm ./OUTSOD_original
spbesod

