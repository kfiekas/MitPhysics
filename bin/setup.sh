#!/bin/bash
#---------------------------------------------------------------------------------------------------
# Setup the MitPhysics package adjusting things that are needed for it to compile and run properly.
#
#                                                                   Jan 16, 2014 - V0 Christoph Paus
#---------------------------------------------------------------------------------------------------

# download the MitPhysics/data directory
$CMSSW_BASE/src/MitPhysics/bin/updateData.sh

# check for existing fastjet+contribution directory or install it
$CMSSW_BASE/src/MitPhysics/bin/installFastjetAndContrib.sh

exit 0
