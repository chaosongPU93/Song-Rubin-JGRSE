#!/bin/bash

# prepare the files for hypoinverse, locate, clean files and generate the event location file for plotting purpose
# this locate the summary phase file 'chat00p.allfam.time' from single phase files like 'chat00p.???.time' from script gentime???.sh
# this deals with original 11 fams
# MAke sure that your disk only has results of those 11 fams when u RUN this script!!! Otherwise delete them first

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
phafmerge="chat00p.allfam.lf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
#phafmerge="unsort.chat00p.allfam.lf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
rm -f $phafmerge
for phafsingle in `ls chat00p.???.lf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}`
#for phafsingle in `ls unsort.chat00p.???.lf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}`
do
    cat $phafsingle >> $phafmerge
done
offfmerge="off.lf.sort.allfam"
#offfmerge="off.lf.allfam"
rm -f $offfmerge
for offfsingle in `ls off.lf.sort???`
#for offfsingle in `ls off.lf.???`
do
    cat $offfsingle >> $offfmerge
done
echo "Merge of unique events of all days at lf of all fams is done!"

cp $rstdir/$phafmerge $locdir
cd $locdir

cp $phafmerge chat00p
./blow.iterchao

evtfile="evtloc.allfam.lf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
#evtfile="unsort.evtloc.allfam.lf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"

cat chat09s | awk '{ printf "%.6f %.6f %.4f \n", -($6 + $7 /60), $4 + $5 /60, $8}' > ${evtfile}.temp
echo "Inversion of unique events of all days at lf of all fams is done!"

#cat $rstdir/off.lf.sort | awk '{print $5, $6, $7}' > time.temp
#paste ${evtfile}.temp time.temp > ${evtfile}
#rm time.temp ${evtfile}.temp

paste ${evtfile}.temp $rstdir/$offfmerge > ${evtfile}
rm ${evtfile}.temp



## HF
cd $rstdir
##phafmerge="chat00p.allfam.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
phafmerge="chat00p.allfam.up.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
#phafmerge="unsort.chat00p.allfam.up.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
rm -f $phafmerge
##for phafsingle in `ls chat00p.???.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}`
for phafsingle in `ls chat00p.???.up.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}`
#for phafsingle in `ls unsort.chat00p.???.up.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}`
do
    cat $phafsingle >> $phafmerge
done
##offfmerge="off.hf.sort.allfam"
offfmerge="off.up.hf.sort.allfam"
#offfmerge="off.up.hf.allfam"
rm -f $offfmerge
for offfsingle in `ls off.hf.sort???`
#for offfsingle in `ls off.hf.???`
do
    cat $offfsingle >> $offfmerge
done
echo "Merge of unique events of all days at hf of all fams is done!"

cp $rstdir/$phafmerge $locdir
cd $locdir

cp $phafmerge chat00p
./blow.iterchao

##evtfile="evtloc.allfam.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
#evtfile="evtloc.allfam.up.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
evtfile="unsort.evtloc.allfam.up.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"

cat chat09s | awk '{ printf "%.6f %.6f %.4f \n", -($6 + $7 /60), $4 + $5 /60, $8}' > ${evtfile}.temp
echo "Inversion of unique events of all days at hf of all fams is done!"

#cat $rstdir/off.lf.sort | awk '{print $5, $6, $7}' > time.temp
#paste ${evtfile}.temp time.temp > ${evtfile}
#rm time.temp ${evtfile}.temp

paste ${evtfile}.temp $rstdir/$offfmerge > ${evtfile}
rm ${evtfile}.temp









