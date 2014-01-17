#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Install the present bambu release.
#
#                                                                              C.Paus (Aug 15, 2010)
#---------------------------------------------------------------------------------------------------
# base checks whether we are ready for install
if [ -z "$MIT_TAG" ] || [ -z "$MIT_VERS" ] || [ -z "$CMSSW_BASE" ]
then
  echo ' Environment is not prepared. Please first setup CMSSW etc.'
  exit 1;
else
  echo "Install BAMBU for MIT_TAG: $MIT_TAG (CMSSW: $CMS_VERS, MIT production: $MIT_VERS)"
  echo ' Environment is not prepared. Please first setup CMSSW etc.'
fi

# make sure there is a valid kerberos ticket
klist -s
if [ "$?" != 0 ]
then
  kinit -f
fi

# go to the proper environment directory
cd $CMSSW_BASE/src

# checkout the relevant code
git clone -b $MIT_TAG https://github.com/cpausmit/MitCommon.git
git clone -b $MIT_TAG https://github.com/cpausmit/MitAna.git
git clone -b $MIT_TAG https://github.com/cpausmit/MitPhysics.git
git clone -b $MIT_TAG https://github.com/cpausmit/MitPlots.git

git clone -b hggpaperV6 https://github.com/bendavid/GBRLikelihood.git HiggsAnalysis/GBRLikelihood

# Initialize the production software
source $CMSSW_BASE/src/MitAna/bin/mitana.sh

# Apply MIT specific patches
./MitPhysics/bin/setup.sh

# build the software
cd $CMSSW_BASE/src
scram build -j8

exit 0
