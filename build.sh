#!/bin/csh -f
# build.sh
# 
#
# Created by Tolulope Olugboji on 11/22/10.
# Copyright 2010 Yale University. All rights reserved.

if ($#argv == 0) then
	echo "No argumets, provide arguments please"
	exit
endif

if ( -e $1 ) then
	echo "Removing existing command files"
	rm $1
endif
	

g++ $1.cpp -o $1 -I/Users/tmo22/Seismology/nr_c304/code/  
