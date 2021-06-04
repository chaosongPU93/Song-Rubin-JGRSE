#!/bin/bash

# generate the input phase file 'chat00p' for hypoinverse scripts written by
# John Armbruster in folder /home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssi
# from a summary offset file combining detections of all days in one fam, result 
# in one file. for fam 002

locdir="$MHOME/Seisbasics/hypoinverse/lzbrst"
stnm1="TWK"
stnm2="LZB"
stnm3="MGC"

sps=40

cd $locdir

## +-1 sample in 40 sps
offfile="offset_13fam_locres_hf1spl"

phafile="chat00p.offset_13fam_locres_hf1spl"
rm -f $phafile
touch $phafile
IFS=$'\n'
for line in `cat $offfile`
do
    #echo $line
    cont12=`echo $line | awk '{ print $1 }'`    # offset for control fam
    cont13=`echo $line | awk '{ print $2 }'`    # offset for control fam
    off12=`echo $line | awk '{ print $3 }'`     # in sample
    off13=`echo $line | awk '{ print $4 }'`     # in sample
    #echo $off12 $off13
    arr2=`echo "scale=4; ($cont12 - $off12) / $sps * 1000 + 30000" | bc | cut -d'.' -f1`
    arr3=`echo "scale=4; ($cont13 - $off13) / $sps * 1000 + 30000" | bc | cut -d'.' -f1`
    evtline1="${stnm1}  s 0 0301010101 30000"
    evtline2="${stnm2}  s 0 0301010101 ${arr2}"
    evtline3="${stnm3}  s 0 0301010101 ${arr3}"
    endline="03010101012000048 3000123 3000-3500"
    echo $evtline1 >> $phafile
    echo $evtline1 >> $phafile
    echo $evtline2 >> $phafile
    echo $evtline2 >> $phafile
    echo $evtline3 >> $phafile
    echo $evtline3 >> $phafile
    echo $endline >> $phafile
done
echo "loc resolution +-1 sample in 40 sps input files for hypoinverse is ready"

cp $phafile chat00p
./blow.iterchao

evtfile="evtloc.13fam_locres_hf1spl"

cat chat09s | awk '{ printf "%.6f %.6f %.4f \n", -($6 + $7 /60), $4 + $5 /60, $8}' > ${evtfile}
echo "Inversion of +-1 sample in 40 sps of all fams is done!"


## +-2 sample in 40 sps
offfile="offset_13fam_locres_hf2spl"

phafile="chat00p.offset_13fam_locres_hf2spl"
rm -f $phafile
touch $phafile
IFS=$'\n'
for line in `cat $offfile`
do
    #echo $line
    cont12=`echo $line | awk '{ print $1 }'`    # offset for control fam
    cont13=`echo $line | awk '{ print $2 }'`    # offset for control fam
    off12=`echo $line | awk '{ print $3 }'`     # in sample
    off13=`echo $line | awk '{ print $4 }'`     # in sample
    #echo $off12 $off13
    arr2=`echo "scale=4; ($cont12 - $off12) / $sps * 1000 + 30000" | bc | cut -d'.' -f1`
    arr3=`echo "scale=4; ($cont13 - $off13) / $sps * 1000 + 30000" | bc | cut -d'.' -f1`
    evtline1="${stnm1}  s 0 0301010101 30000"
    evtline2="${stnm2}  s 0 0301010101 ${arr2}"
    evtline3="${stnm3}  s 0 0301010101 ${arr3}"
    endline="03010101012000048 3000123 3000-3500"
    echo $evtline1 >> $phafile
    echo $evtline1 >> $phafile
    echo $evtline2 >> $phafile
    echo $evtline2 >> $phafile
    echo $evtline3 >> $phafile
    echo $evtline3 >> $phafile
    echo $endline >> $phafile
done
echo "loc resolution +-2 sample in 40 sps input files for hypoinverse is ready"

cp $phafile chat00p
./blow.iterchao

evtfile="evtloc.13fam_locres_hf2spl"

cat chat09s | awk '{ printf "%.6f %.6f %.4f \n", -($6 + $7 /60), $4 + $5 /60, $8}' > ${evtfile}
echo "Inversion of +-2 sample in 40 sps of all fams is done!"






