#!/bin/bash

rm -f ENERGIES 
ls -l *.gout |awk '{print "energulp ", $8}' > rungetener
chmod +x rungetener
./rungetener
rm rungetener

