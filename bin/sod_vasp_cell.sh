#!/bin/bash

rm -f cell.dat a.dat b.dat c.dat

n_columns_ls=`ls -l |tail -1 |awk '{ FS = "|" } ; { print NF}'`
ls -l */CONTCAR |awk -v nc=$n_columns_ls '{print "cellvasp.sh", $nc}' > rungetcell
chmod +x rungetcell
./rungetcell
paste  a.dat b.dat c.dat > cell.dat
rm rungetcell a.dat b.dat c.dat 


# Cubic structure
echo "1" > DATA
awk '{print $1*$5*$9}' cell.dat >> DATA


