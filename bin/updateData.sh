#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Update the data directory in the MitPhysics package. The data are to large to maintain in GITHUB
# therefore they are stored on a server from where they will be downloaded using this script.
#
# It is important that CMSSW is setup ($CMSSW_BASE/MitPhysics/data is the target). The target
# directory will be created if it does not exist and the content will be overwritten with the
# current server content.
#
#                                                                   Jan 12, 2014 - V0 Christoph Paus
#---------------------------------------------------------------------------------------------------
# define what we want to download
SERVER=t3serv001.mit.edu
LOCATION='~cmsprod'
TAR='MitPhysics_data.tgz'

# test whether environment is clean
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

# go to our target area
cd $CMSSW_BASE/src/

# is there already a version of this file?
if [ -e "$TAR" ]
then
  echo ""
  echo -n " WARNING - $TAR exists already. Overwrite it? [N/y] "
  read answer
  if [ "$answer" != "y" ]
  then 
    echo ""
    echo " YOUR ANSWER: $answer -> exit now without further action."
    echo ""
    exit 0
  else
    echo ""
  fi
fi

# now do the download
echo " INFO - download starting"
echo ""
wget "http://$SERVER/$LOCATION/$TAR" -O $TAR

# unpack it (no debug, could add 't' option)
echo " INFO - unpacking"
echo ""
tar fzx $TAR

# cleanup to avoid junk
echo " INFO - cleaning up"
echo ""
rm -rf $TAR

# exit clean
exit 0
