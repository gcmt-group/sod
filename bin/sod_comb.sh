#!/bin/bash
# This script generates the correct names for the input files 
# It also generates a script that can be used to run all the simulations

rm -f OUTSOD SUPERCELL EQMATRIX OPERATORS cSGO 
rm -rf CALCS

clear

combsod

FILER=$(<filer)

if [ $FILER -ne 0 ]; then
  mkdir CALCS
  cd CALCS
  mv ../fort.* .
  cp ../OUTSOD .

  # FILER=1 for GULP
  if [ $FILER -eq 1 ];  then
    extin="gin"
    extout="gout"
    program="gulp"
  fi  

  # FILER=2 for METADISE
  if [ $FILER -eq 2 ]; then
    extin="min"
    extout="mout"
    program="metadise"
  fi 

  # FILER=11 for VASP
  if [ $FILER -eq 11 ]; then
    extin="vasp"
    extout="vout"
    program="vasp"
  fi

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

 n_columns_ls=`ls -l |tail -1 |awk '{ FS = "|" } ; { print NF}'`
 ls -l *$extin |awk -v nc=$n_columns_ls -v extout=$extout -v program=$program '{print program,"<",$nc,">",$nc"."extout}' |sed s/$extin.$extout/$extout/  > job_sender
  chmod +x job_sender

fi

cd ..
rm -f filer coord* 
