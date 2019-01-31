#!/bin/bash

rm -f tmp* errors.txt

###### Do one calculation for each sgo file
for i in $( ls *.bcs |sed s/.bcs// ) ; do
echo Space group: $i


###### Read the cif file
sed s/\,/\ /g $i.bcs | awk '{print $2,$3,$4}' > tmp1
#awk '{print $2,$3,$4}' $i.bcs > tmp1

noperators=`grep -c "" tmp1 `
echo $noperators operators read.
echo 

awk -v spgr=$i -v noperators=$noperators   -f $ROOTSOD/bin/bcs2sgo.awk tmp1 > tmp2

cp tmp2 tmpf
mv tmpf $i.sgo
echo "0" >> $i.sgo

grep ERROR $i.sgo >> errors.txt

done


###### Error reporting
exist_errors=`grep -c ERROR errors.txt`
if [ $exist_errors == 0 ] ; then

 rm errors.txt
 echo
 echo SGO file successfully created.
 echo

else

 echo $exist_errors
 echo "ERROR in conversion script. Operator not read. Edit the bcs2sgo.awk script to solve this format problem."
 echo List of ERRORS:
 cat errors.txt

fi

rm tmp1 tmp2

