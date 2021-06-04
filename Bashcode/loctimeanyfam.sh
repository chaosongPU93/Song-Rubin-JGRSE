#!/bin/bash

# prepare the files for hypoinverse, locate, clean files and generate the event location file for plotting purpose
# this locate the summary phase file 'chat00p.allfam.time' from single phase files like 'chat00p.???.time' from script gentime???.sh
# this deals with fam as requested in the 'famlist'
# RUN this ONLY after results of new fam are obtained

# different from other loctimeXXX, this one reads the fams needed to combined from a list file

rstdir="$ALLAN/data-no-resp/LZBtrio/MAPS"
locdir="$MHOME/Seisbasics/hypoinverse/lzbrst"

# HF parameters
hihf="6.5"
lohf="1.25"
wlenhf="4"

# LF parameters
lolf="0.5"
hilf="1.25"
#wlenlf="24"
#wlenlf="12"  # several choices
#wlenlf="20"
wlenlf="16"
lofflf="4"
ccminlf="0.35"


## LF
cd $rstdir
phafmerge="chat00p.requestfam.lf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
rm -f $phafmerge

for fam in `cat $locdir/famlist`
do
   phafsingle="chat00p.${fam}.lf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
   #echo $phafsingle
   cat $phafsingle >> $phafmerge
done

offfmerge="off.lf.sort.requestfam"
rm -f $offfmerge

for fam in `cat $locdir/famlist`
do
    offfsingle="off.lf.sort${fam}"
    #echo $offfsingle
    cat $offfsingle >> $offfmerge
done
echo "Merge of unique events of all days at lf of all fams is done!"

cp $rstdir/$phafmerge $locdir
cd $locdir

cp $phafmerge chat00p
./blow.iterchao

evtfile="evtloc.requestfam.lf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"

cat chat09s | awk '{ printf "%.6f %.6f %.4f \n", -($6 + $7 /60), $4 + $5 /60, $8}' > ${evtfile}.temp
echo "Inversion of unique events of all days at lf of all fams is done!"

paste ${evtfile}.temp $rstdir/$offfmerge > ${evtfile}
rm ${evtfile}.temp



## HF
cd $rstdir
phafmerge="chat00p.requestfam.up.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
rm -f $phafmerge

for fam in `cat $locdir/famlist`
do
    phafsingle="chat00p.${fam}.up.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
    cat $phafsingle >> $phafmerge
done

offfmerge="off.up.hf.sort.requestfam"
rm -f $offfmerge

for fam in `cat $locdir/famlist`
do
    offfsingle="off.hf.sort${fam}"
    cat $offfsingle >> $offfmerge
done
echo "Merge of unique events of all days at hf of all fams is done!"

cp $rstdir/$phafmerge $locdir
cd $locdir

cp $phafmerge chat00p
./blow.iterchao

evtfile="evtloc.requestfam.up.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"


cat chat09s | awk '{ printf "%.6f %.6f %.4f \n", -($6 + $7 /60), $4 + $5 /60, $8}' > ${evtfile}.temp
echo "Inversion of unique events of all days at hf of all fams is done!"

paste ${evtfile}.temp $rstdir/$offfmerge > ${evtfile}
rm ${evtfile}.temp









