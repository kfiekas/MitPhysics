#!/bin/bash
#---------------------------------------------------------------------------------------------------
#
# Upload a new version of the data directory used in the MitPhysics package. The package will have
# to be carefully build we do not takje responsibility for that. The specified file will make a
# backup of the existing file (with todays date) and upload the new version. Then all other
# maintained copies will be updated
#
#                                                                   May 29, 2014 - V0 Christoph Paus
#---------------------------------------------------------------------------------------------------
# define what we want to download
SERVER=t3serv001.mit.edu
LOCATION='~cmsprod'
TAR='MitPhysics_data.tgz'

date=`date +%Y.%m.%d-%H:%M:%S`

# command line arguments
SOURCE="$1"

# backup existing version on the main server
ssh cmsprod@t3serv001.mit.edu mv public_html/$TAR public_html/${TAR}-$date
ssh cmsprod@t3serv001.mit.edu ls -lhrt public_html/${TAR}\*

# upload new version in proper location for web access
scp $SOURCE cmsprod@t3serv001.mit.edu:public_html/$TAR

# upload new version in proper location for Tier-2 and Tier-3 access
scp $SOURCE t3btch100.mit.edu:/mnt/hadoop/cms/store/user/paus/$TAR
scp $SOURCE se01.cmsaf.mit.edu:/mnt/hadoop/cms/store/user/paus/$TAR

# exit clean
exit 0
