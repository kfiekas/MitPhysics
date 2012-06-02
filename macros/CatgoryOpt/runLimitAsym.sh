#! /bin/bash
which=$1
numCat=$2
catVal1=$3
catVal2=$4
catVal3=$5
catVal4=$6
catVal5=$7

if [ "$catVal1" == "" ] || [ $numCat -lt 2 ]
    then
    catVal1=-100
fi

if [ "$catVal2" == "" ] || [ $numCat -lt 3 ]
    then
    catVal2=-100
fi

if [ "$catVal3" == "" ] || [ $numCat -lt 4 ]
    then
    catVal3=-100
fi

if [ "$catVal4" == "" ] || [ $numCat -lt 5 ]
    then
    catVal4=-100
fi

if [ "$catVal5" == "" ] || [ $numCat -lt 6 ]
    then
    catVal5=-100
fi


sourceDir=/afs/cern.ch/user/m/mingyang/optimizeMVACats
thisDir=$(pwd)

# set up correct ROOT version
# ming: do i need to do set up ROOT every time?
cd /afs/cern.ch/user/m/mingyang/CMSSW_4_2_3_patch3/src
eval `scram ru -sh`
cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.28.00b/x86_64-slc5-gcc43-opt/root
. bin/thisroot.sh
cd $thisDir
cp $sourceDir/doFit_asymptotic.C ./

#root -l -q -b $thisDir"/doFit_asymptotic.C+("$which","$numCat","$catVal1","$catVal2","$catVal3","$catVal4","$catVal5")"
root -l -q -b $thisDir"/doFit_asymptotic.C+("$which","$numCat","$catVal1","$catVal2","$catVal3","$catVal4","$catVal5")"