#!/bin/bash

ls -l *.mout |awk '{print "./getcell " $9}' > rungetcell
chmod +x rungetcell
rm a.dat b.dat c.dat
./rungetcell
paste  a.dat b.dat c.dat > cell.dat
awk '{print $1*$5, $1*$5*$9}' cell.dat > av.dat

