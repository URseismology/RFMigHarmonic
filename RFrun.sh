#!/bin/csh -f
# script to run multiple RF analyses
# Currently built to test my code which re-implents JJParks 2000 analysis
# Arrays are all in the file Events.in
# Other issues to revisit: Output options, and azimuth options. see RecFunc.cpp
#

make >& /dev/null 
RecFunc -LCleanEvents.in -F1 -T80 -ORF_Azim.grid -B2 
