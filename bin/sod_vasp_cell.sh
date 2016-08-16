#!/bin/bash

rm -f cell.dat a.dat b.dat c.dat
ls -l c??/CONTCAR |awk '{print "cellvasp", $9}' > rungetcell
chmod +x rungetcell
./rungetcell
paste  a.dat b.dat c.dat > cell.dat
rm rungetcell a.dat b.dat c.dat 


# Cubic structure
echo "1" > DATA
awk '{print $1*$5*$9}' cell.dat >> DATA


