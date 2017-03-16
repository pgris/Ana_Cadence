#!/bin/bash

fieldname=$1
fieldid=$2
season=$3
nevts=$4
sntype=$5
versionnum=$6
rolling=$7
nmerger=$8
percent=$9



thedir=/sps/lsst/data/dev/pgris/Make_Cadence/Sim_minion_1016

prefix=SuperNova
addit=''
rrolling=''
if [ "$rolling" == "1" ]
then
prefix=${prefix}_Rolling
addit=${nmerger}_${percent}_
rrolling=_Rolling
fi
output=List_${sntype}${rrolling}_${fieldname}_${fieldid}_${nevts}_${addit}season_${season}_minion.dat

#ls ${thedir}/SuperNova_${sntype}_${fieldname}_${fieldid}_*_${nevts}_season_${season}_*.pkl | cut -d '/' -f2 >& ${output}
ls ${thedir}/${prefix}_${sntype}_${fieldname}_${fieldid}_*_${nevts}_${addit}season_${season}_${versionnum}.pkl >& ${output}
