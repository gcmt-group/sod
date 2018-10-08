#!/bin/bash

n_columns_ls=`ls -l |tail -1 |awk '{ FS = "|" } ; { print NF}'`
ls -l *.gout |awk -v nc=$n_columns_ls '{print "cellgulp.sh", $nc}' > rungetcell
chmod +x rungetcell
./rungetcell
paste  a.dat b.dat c.dat alpha.dat beta.dat gamma.dat volume.dat > cell.dat
rm rungetcell a.dat b.dat c.dat alpha.dat beta.dat gamma.dat volume.dat

#Hexagonal structure
awk '{print $1*$2*sin((3.1415926)*$6/(180.0)), $7}' cell.dat > av.dat

