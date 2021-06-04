#!/bin/bash

# prepare the files for hypoinverse, locate, clean files and generate the event location file for plotting purpose
# this locate the summary phase file from gen-hypo-input-alldays.sh
# this deals with original 11 fams plus a new fam 006
# RUN this ONLY after results of new fam are obtained

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


# LF
cd $rstdir
phafmerge="chat00p.all12fam.lf.count.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
rm -f $phafmerge
for phafsingle in `ls chat00p.???.lf.count.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}`
do
    cat $phafsingle >> $phafmerge
done
offfmerge="countori_lf.all12fam"
rm -f $offfmerge
for offfsingle in `ls countori_???.*${lofflf}*${ccminlf}*_${lolf}-${hilf}_${wlenlf}*`
do
    cat $offfsingle >> $offfmerge
done
echo "Merge of unique events of all days at lf of all fams is done!"

cp $rstdir/$phafmerge $locdir
cd $locdir

cp $phafmerge chat00p
./blow.iterchao

evtfile="evtloc.all12fam.lf.count.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"

cat chat09s | awk '{ printf "%.6f %.6f %.4f \n", -($6 + $7 /60), $4 + $5 /60, $8}' > ${evtfile}.temp
echo "Inversion of unique events of all days at lf of all fams is done!"

#cat $rstdir/$offfile | awk '{print $5}' > count.temp
#paste ${evtfile}.temp count.temp > ${evtfile}
#rm count.temp ${evtfile}.temp

paste ${evtfile}.temp $rstdir/$offfmerge > ${evtfile}
rm ${evtfile}.temp


# HF
cd $rstdir
#phafmerge="chat00p.allfam.hf.count.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
phafmerge="chat00p.all12fam.up.hf.count.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
rm -f $phafmerge
#for phafsingle in `ls chat00p.???.hf.count.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}`
for phafsingle in `ls chat00p.???.up.hf.count.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}`
do
    cat $phafsingle >> $phafmerge
done
#offfmerge="countori_hf.allfam"
offfmerge="countori_up.hf.all12fam"
rm -f $offfmerge
#for offfsingle in `ls countori_???.*_${lohf}-${hihf}_${wlenhf}*`
for offfsingle in `ls countori_???.up.*_${lohf}-${hihf}_${wlenhf}*`
do
    cat $offfsingle >> $offfmerge
done
echo "Merge of unique events of all days at hf of all fams is done!"

cp $rstdir/$phafmerge $locdir
cd $locdir

cp $phafmerge chat00p
./blow.iterchao

#evtfile="evtloc.allfam.hf.count.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
evtfile="evtloc.all12fam.up.hf.count.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"

cat chat09s | awk '{ printf "%.6f %.6f %.4f \n", -($6 + $7 /60), $4 + $5 /60, $8}' > ${evtfile}.temp
echo "Inversion of unique events of all days at hf of all fams is done!"

#cat $rstdir/$offfile | awk '{print $5}' > count.temp
#paste ${evtfile}.temp count.temp > ${evtfile}
#rm count.temp ${evtfile}.temp

paste ${evtfile}.temp $rstdir/$offfmerge > ${evtfile}
rm ${evtfile}.temp
