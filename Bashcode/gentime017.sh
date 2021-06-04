#!/bin/bash

# generate the input phase file 'chat00p' for hypoinverse scripts written by
# John Armbruster in folder /home/data2/chaosong/Seisbasics/hypoinverse/LOCpgsssi
# from a summary offset file combining detections of all days in one fam, result 
# in one file. for fam 017

rstdir="$ALLAN/data-no-resp/LZBtrio/MAPS"
stnm1="TWK"
stnm2="LZB"
stnm3="MGC"
fam="017"

# HF parameters
sps=40
hihf="6.5"
lohf="1.25"
wlenhf="4"

# LF parameters
hilf="1.25"
lolf="0.5"
#wlenlf="24"
#wlenlf="12"  # several choices
#wlenlf="20"
wlenlf="16"
lofflf="4"
ccminlf="0.35"

# HF
cd $rstdir
#offfile=`ls timeori_${fam}.*_${lohf}-${hihf}_${wlenhf}*`
offfile=`ls timeori_${fam}.up.*_${lohf}-${hihf}_${wlenhf}*`

### not sure why i sorted it, seems unnecessary, now testing if the same without sorting,
### change several lines accordingly, 2020/05/18
#cat $offfile | sort -t ' ' -nk 5,5 -nk 6,6 -nk 7,7 -nk 1,1 -nk 2,2 -o off.hf.sort${fam} 
cat $offfile > off.hf.${fam}

##phafile="chat00p.${fam}.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
#phafile="chat00p.${fam}.up.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
phafile="unsort.chat00p.${fam}.up.hf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
rm -f $phafile
touch $phafile
IFS=$'\n'
#for line in `cat off.hf.sort${fam}`
for line in `cat off.hf.${fam}`
do
    #echo $line
    off12=`echo $line | awk '{ print $3 }'`     # in sec
    off13=`echo $line | awk '{ print $4 }'`     # in sec
    #echo $off12 $off13
    arr2=`echo "scale=4; (-18 / $sps - $off12) * 1000 + 30000" | bc | cut -d'.' -f1`
    arr3=`echo "scale=4; (6 / $sps - $off13) * 1000 + 30000" | bc | cut -d'.' -f1`
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
echo "HF input files for hypoinverse is ready"

# LF
offfile=`ls timeori_${fam}.*${lofflf}*${ccminlf}*_${lolf}-${hilf}_${wlenlf}*`
#echo $offfile

#cat $offfile | sort -t ' ' -nk 5,5 -nk 6,6 -nk 7,7 -nk 1,1 -nk 2,2 -o off.lf.sort${fam}
cat $offfile > off.lf.${fam}

#phafile="chat00p.${fam}.lf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
phafile="unsort.chat00p.${fam}.lf.time.${wlenhf}_${wlenlf}.${lofflf}.${ccminlf}"
rm -f $phafile
touch $phafile
IFS=$'\n'
#for line in `cat off.lf.sort${fam}`
for line in `cat off.lf.${fam}`
do
    #echo $line
    off12=`echo $line | awk '{ print $3 }'`     # in sec
    off13=`echo $line | awk '{ print $4 }'`     # in sec
    #echo $off12 $off13
    arr2=`echo "scale=4; (-18 / $sps - $off12) * 1000 + 30000" | bc | cut -d'.' -f1`
    arr3=`echo "scale=4; (6 / $sps - $off13) * 1000 + 30000" | bc | cut -d'.' -f1`
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
echo "LF input files for hypoinverse is ready"


