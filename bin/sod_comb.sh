#!/bin/bash
# This script generates the correct names for the input files as well as
# the executable that will run all the simulations

rm  matrix* coordinates.xyz  OUTSOD supercell
rm -r CALCS

clear

combsod

set FILER=`cat  filer`

if ($FILER > 0) then

mkdir CALCS
cd CALCS
mv fort.* .
cp OUTSOD .

# FILER=1 for GULP
if ($FILER == 1)  then
set extin=gin
set extout=gout
set program=gulp
endif

# FILER=2 for METADISE
if ($FILER == 2) then
set extin=min
set extout=mout
set program=metadise
endif

# FILER=11 for VASP
if ($FILER == 11)  then
set extin=vaspin
set extout=vaspout
set program=vasp
endif

ls fort.* > tmp1

sed s/fort.// tmp1 |awk '{printf ("%+4s\n",$1-100000)}' > tmp3

awk '{ if ($1<10) {print "0000"$1}}'                 tmp3 > tmp4
awk '{ if (($1>9)&&($1<100)) {print "000"$1}}'       tmp3 >>tmp4
awk '{ if (($1>99)&&($1<1000)) {print "00"$1}}'      tmp3 >>tmp4
awk '{ if (($1>999)&&($1<10000)) {print   "0"$1}}'   tmp3 >>tmp4
awk '{ if (($1>9999)&&($1<100000)) {print    $1}}'   tmp3 >>tmp4
awk '{ if ($1>99999) {print    "Error in file numbering, too many configurations (> SOD limit = 99999)"}}'  tmp3

awk -v extin=$extin '{print "c"$1"."extin}' tmp4 > tmp5

paste tmp1 tmp5 > tmp6
awk '{print "mv",$0}' tmp6 > tmp7
chmod +x tmp7
./tmp7
rm tmp*

# This is to create the script that is going to run the calculations, with appropriate extension
ls -l *$extin | awk -v extout=$extout -v program=$program '{print program,"<",$9,">",$9"."extout}' |sed s/$extin.$extout/$extout/  > job_sender
chmod +x job_sender

cd ..

endif

rm filer
rm coord* supercell matrix* operat* conf* indc* 

