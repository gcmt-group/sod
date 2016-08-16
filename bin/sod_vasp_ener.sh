#!/bin/bash

rm -f ENERGIES 
ls -l c??/OUTCAR |awk '{print "enervasp ", $9}' > rungetener
chmod +x rungetener
./rungetener
rm rungetener

