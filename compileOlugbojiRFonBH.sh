#!/bin/csh -f
# Author: Tolulope Olugboji
# Updated on June 7, 2018 - script written to compile code on Pippin machine - use to compile @UofR
# script to Debug and Iteratively develop the Harmonic stack mode: -M option
#    Look at Code Blog to follow thoughts....
# Currently built to implement methods in JJParks 2010 analysis and legacy code rfmig_cboot.f
#
# From legacy code to version control in git. Repos. Clone the following repos to get all src files
#
# ----------------------------------------- REPO LIST
#  git clone https://tolugboji@bitbucket.org/tolugboji/seismorfcodes.git
#  git clone https://tolugboji@bitbucket.org/tolugboji/numrecipesc.git
#
# ------------------------------------------ END REPO LIST

# ----------------------------------------- CHANGE FOLLOWING MAKEFILE INCLUDES
#
# HOME =/home/jamespippin/Documents/1.UofR_Seismology/numrecipesc
# PWD = ${HOME}/code/
#
# ------------------------------------------ END CHANGE FILE INCLUDES




# Output Files returned on execution of script
#                     1. Azimuth_Rad.xyz     Azimuth_Trans.xyz
#                     2. OutputFName.log      : Stack & Bin Statistics
#                     3. OutputFNameEvent.txt : Event Statistics
#											  : Azimuth, EvLong,Evlat,Mag,Dist.

# Never Push this Code to the repo.
# Either I always put out two different


set rootDir = /scratch/tolugboj_lab/softwares/rfcodes
#sset rootDir = /home/jamespippin/Documents/1.UofR_Seismology/seismorfcodes
set Distrib = $rootDir/RFMigHarmonic_bl
set LocalMake = $Distrib/Makefile
set isRForHK = 1   # If 1 then Compile RecFunc, Otherwise Compile HK Stacking Code



if ($isRForHK == 1 ) then
make help -f $LocalMake
make RecFunc -f $LocalMake
else
make help -f $LocalMake
make RFVelStck -f $LocalMake

### Independent tests for now... Working on Sequential HkVk Stacks
set makeStatus = $status

if ( $makeStatus == 2 ) then
echo "Error in Make file revisit its execution"
else
echo  "Compilation Success... Running Basics"

endif

goto finish

endif

set makeStatus = $status

if ( $makeStatus == 2 ) then
echo "Error in Make file revisit its execution"
endif


cp $Distrib/RecFunc $rootDir/bin/RecFuncHarmonic_BlueHive

finish:



