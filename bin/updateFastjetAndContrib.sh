#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Update fastjet and all relevant contributions to be able to fully compile MitPhysics.
# Here we make sure to cleanup old links in the $CMSSW_BASE/external area and update them
#
#
#                                                                   May 29, 2014 - V0 LDM
#---------------------------------------------------------------------------------------------------
function reconfigureScram {
  FASTJET_VAR=$1
  EXTERNAL=$2
  # cleanup old links
  rm -rf $CMSSW_BASE/external/$SCRAM_ARCH/*
  # backup scram configuration file 
  mv $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml \
     $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml-backup
  # add dummy location to enforce srcram setup refresh 
  DUMMYEXTERNAL=/home/dimatteo/cms/dummyexternal
  echo \
'
  <tool name="fastjet" version="xx-VERSION-xx">
    <info url="http://www.lpthe.jussieu.fr/~salam/fastjet/"/>
    <lib name="fastjetplugins"/>
    <lib name="fastjettools"/>
    <lib name="siscone"/>
    <lib name="siscone_spherical"/>
    <lib name="fastjet"/>
    <lib name="fastjetcontrib"/>
    <client>
      <environment name="FASTJET_BASE" default="xx-PATH-xx"/>
      <environment name="LIBDIR" default="$FASTJET_BASE/lib"/>
      <environment name="INCLUDE" default="$FASTJET_BASE/include"/>
    </client>
  </tool>
' | sed "s/xx-VERSION-xx/$FASTJET_VER/"  | sed "s#xx-PATH-xx#$DUMMYEXTERNAL#" \
> $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml

  # refresh scram
  cd $CMSSW_BASE/src
  scram setup fastjet

  # restore the correct config file and reconfigure scram one last time 
  mv $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml-backup \
     $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml
  scram setup fastjet
}

# test whether environment is clean and we can start
if [ -z "$CMSSW_BASE" ]
then
  echo ""
  echo " ERROR - cmssw is not setup. Please, setup release directory and checkout MitPhysics."
  echo ""
  exit 1
fi
if [ -d "$CMSSW_BASE/src/MitPhysics" ]
then
  echo ""
  echo " INFO - found MitPhysics location at: $CMSSW_BASE/src/MitPhysics"
  echo ""
else
  echo ""
  echo " ERROR - MitPhysics is not in your release. Please check it from GITHUB/CVS."
  echo ""
  exit 1
fi

# the external library location
EXTERNAL="/home/cmsprod/cms/external"
if [ -d "$EXTERNAL" ]
then
  echo ""
  echo " INFO - updating links to external libraries located at: $EXTERNAL"
  echo ""
  # use existing location just adjust scram configuration
  FASTJET_VER=`ls -1 $EXTERNAL | grep fastjet | tail -1 |cut -d '-' -f2`
  reconfigureScram $FASTJET_VER $EXTERNAL
  exit 0
else
  echo " ERROR - trying to update links to non-existing location $EXTERNAL ."
  echo ""
  exit 1
fi

exit 0
