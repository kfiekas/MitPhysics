#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Install fastjet and all relevant contributions to be able to fully compile MitPhysics.
#
# This situation should be temporary but for now it is required to work with the contributions to
# fastjet and the most up to date versions.
#
#                                                                   Jan 13, 2014 - V0 Christoph Paus
#---------------------------------------------------------------------------------------------------
function configureScram {
  FASTJET_VAR=$1
  EXTERNAL=$2
  # add local fastjet external to scarm config
  mv $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml \
     $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml-last.$$
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
' | sed "s/xx-VERSION-xx/$FASTJET_VER/"  | sed "s#xx-PATH-xx#$EXTERNAL#" \
> $CMSSW_BASE/config/toolbox/$SCRAM_ARCH/tools/selected/fastjet.xml

  # commit the scram config changes
  cd $CMSSW_BASE/src
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

# the default location
EXTERNAL_DEF="/home/cmsprod/cms/external"
if [ -d "$EXTERNAL_DEF" ]
then
  EXTERNAL="$EXTERNAL_DEF"
  echo ""
  echo " INFO - found external location at: $EXTERNAL"
  echo ""
  # use existing location just adjust scram configuration
  FASTJET_VER=`ls -1 $EXTERNAL | grep fastjet | tail -1 |cut -d '-' -f2`
  configureScram $FASTJET_VER $EXTERNAL
  exit 0
else
  EXTERNAL="/home/$USER/cms/external"
  echo " INFO - default external location ($EXTERNAL_DEF) not found make own external at: $EXTERNAL"
  echo ""
fi

# Here the real work starts (fastjet, fjcontrib, tweaks)

mkdir -p $EXTERNAL

# install best fastjet version
#-----------------------------

# set all relevant variables
FASTJET_URL="http://fastjet.fr/repo"
FASTJET_TGZ="fastjet-3.0.6.tar.gz"
FASTJET_DIR=`echo $EXTERNAL/$FASTJET_TGZ | sed 's/.tar.gz//'`
FASTJET_VER=`echo $FASTJET_TGZ | sed 's/.tar.gz//' |cut -d'-' -f2`

# in the right location
cd $EXTERNAL

# now do the download
echo " INFO - download starting"
echo ""
wget "$FASTJET_URL/$FASTJET_TGZ" -O $FASTJET_TGZ

# unpack it (no debug, could add 't' option)
echo " INFO - unpacking"
echo ""
tar fzx $FASTJET_TGZ

# cleanup to avoid junk
echo " INFO - cleaning up"
echo ""
rm -rf $FASTJET_TGZ

# installing fastjet
echo " INFO - installing"
echo ""
cd `echo $FASTJET_TGZ | sed 's/.tar.gz//'`
./configure --prefix=$EXTERNAL
make
make check
make install


# install best fastjet contribution version
#------------------------------------------

# define all relevant variables
FJCONTRIB_URL="http://t3serv001.mit.edu/~cmsprod"
FJCONTRIB_TGZ="fjcontrib-1.011_nsub-2.0.0-rc3.tar.gz"
FJCONTRIB_DIR=`echo $EXTERNAL/$FJCONTRIB_TGZ | sed 's/.tar.gz//'`

# in the right location
cd $EXTERNAL

# now do the download
echo " INFO - download starting"
echo ""
wget "$FJCONTRIB_URL/$FJCONTRIB_TGZ" -O $FJCONTRIB_TGZ

# unpack it (no debug, could add 't' option)
echo " INFO - unpacking"
echo ""
tar fzx $FJCONTRIB_TGZ

# cleanup to avoid junk
echo " INFO - cleaning up"
echo ""
rm -rf $FJCONTRIB_TGZ

# installing fastjet
echo " INFO - installing"
echo ""
cd `echo $FJCONTRIB_TGZ | sed 's/.tar.gz//'`
./configure --fastjet-config=$FASTJET_DIR/fastjet-config --prefix=$EXTERNAL CXXFLAGS="-O3 -Wall -Woverloaded-virtual -g -fPIC -I$EXTERNAL/include"
make
make check
make install

# make shared libraries for the contributions
g++ -shared -fPIC -o $EXTERNAL/lib/libfastjetcontrib.so -Wl,-soname,libfastjetcontrib.so $FJCONTRIB_DIR/*/[A-Z]*.o 

# final adjustment to scram configuration
configureScram $FASTJET_VER $EXTERNAL

exit 0
