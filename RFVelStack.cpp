/**  RFVelStack.cpp 
 *
 *  @author  Created by Tolulope  Olugboji (tolumorayo@gmail.com) on 7/23/13.
 ** @section LICENSE
 *
 *  Copyright (c) 2013 Yale University. All rights reserved.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * @section DESCRIPTION
 *
 *	This code is the first code for interpreting RF data.
 *	It does a grid search on elastic parameters and depth that optimizes the
 *	migrated stacks - Epicentral migrated stacks.
 *					  Azimuthal migrated stacks.
 *					  Harmonic migrated stacks.
 *

 * To use we need - the timingWindow of pulse to migrate -
 *					the seismograms - event list used to compute the RFs
 *					phase to be used in migration.
					bounds of grid search to parameterize migration
 */

#include <iostream>
#include "error.h"
#include <cstdio>
#include <cctype>
#include "GrdSrchRouter.h"


using namespace std;

void usage_and_exit(){
	cout << "\n\n\n***ScatteredWaveAnalysis(Grid Search Sequential HK Velocity Analysis) Version 1.0. *** \nCreated by Tolulope Olugboji. \nCopyright 2013 Yale University. All rights reserved.  \nUpdated July 1, 2013. \n \n";
	cout << "Usage: RFVelStack  [options] \n \n";
	
	
	cout << "-L \t File name holding the [L]ist of SAC records. The leading [.RTZ] attached internally \n";
	cout << "-O \t File name where the [O]utput grd files should be stored. \n";
	cout << "-S \t File name holding the [S]earch parameters. H -K bounds. \n";
	cout << "-I [incDecim] \t incDecim, number to increase decimation and use bicubic spline interpolation before search. \n";
	cout << "-T \t [T]ime window in seconds. Default 50 secs.\n";
	cout << "-F \t Maximum [F]requency for computing the Receiver Functions. Default 1 Hz. \n";
	cout << "-E(bazmax/bazmin/epicmin/epicmax) \t Parameters for restricting stack to earthquakes in particular location.\n\t\t\t bazmax: left back azimuth \n\t\t\t bazmin: right back azimuth \n\t\t\t epicmin: minimum distance in deg. \n\t\t\t epicmax: maximum distance in degrees \n";
	cout << "-H \t Use [H]eader tag look up. IF time tag is NOT set, skip record.  \n";
	cout << "**new: -Z  \t activate this flag to use frequency domain migrtion in stack \n";
	cout << "-R \t Set 1 if LQT [R]otation is to be done.  \n";
	cout << "-V \t Increase verbosity \n" << endl;
	
	
	
	exit(1);
	
}

int main(int argc, char** argv){
	
	// check if all is running well
	bool err = false; bool getHelp = false;
	
	// Variables to hold option flags. 9 Options in Total. Group by Requisites
	bool getEvents=false, getOut = false, getSearchParams=false; // 3. Must Haves.
	bool getTimewin=false, getFreq=false, getEpic=false;		 // If not set default
	bool getHeaderTag=false, getLqtSwch=false,  getVerbose=false;      // If not set default.
	
	bool getInterpolate=false;   // Interpolate. Improves resolution - maybe?
	bool getFreqMigrateFlag = false; // Set flag to compute stack in migration domain?
	
	// Variables to support options: additional option arguments. 
	// Modelled after RecFunc but now RFVelStack... 
	char *eventsfn,  *outfn, *searchParamsfn;
	float timeWin, Fcutoff, epbazmax, epbazmin, epicmin, epicmax;
	float timeTag, lqtRotArg;
	int nDecim;   // factor to improve nGrid points.
	
	//Retrieve arguments and ensure array file well parsed.
	for (int i=1; i<argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
					//1.
				case 'L':
					eventsfn = &argv[i][2];
					getEvents = true;
					break;
					//2.
				case 'T':
					getTimewin = true;
					timeWin = atof(&argv[i][2]);
					break;
					//3.
				case 'F':
					Fcutoff = atof(&argv[i][2]);
					getFreq = true;
					break;
					//4.
				case 'E':
					if (sscanf(&argv[i][2], "%f/%f/%f/%f", 
							   &epbazmin, &epbazmax, &epicmin, &epicmax) == 4) {
						getEpic = true;
					}else {
						error("Invalid -E option. Using default Epic. setting");
					}
					break;
					//5.
				case 'O':
					outfn = &argv[i][2];
					getOut = true;
					break;
					//6.
				case 'H':
					getHeaderTag = true;
					timeTag = atoi(&argv[i][2]);
					break;
					//7.
				case 'R':
					getLqtSwch = true;
					lqtRotArg = atof(&argv[i][2]);
					break;
					//8.
				case 'V':
					getVerbose = true;
					break;
					//9.
				case 'S':
					searchParamsfn = &argv[i][2];
					getSearchParams = true;
					break;
				case 'I':
					nDecim = atoi(&argv[i][2]);
					getInterpolate = true;
					break;
				case 'Z':
					getFreqMigrateFlag = true;
					break;
				case 'h':
					getHelp = true;
					break;
				default:
					err=true;
					break;
			}
		}else {
			err = true;
		}
		
	}
	
	// These 3 Options are the minimum options needed to be set for code to run
	if (!getEvents || !getOut || !getSearchParams || getHelp) {
		err = true;
	}
	
	//These 3 Options  Data Resolution options
	if (!getTimewin || !getFreq || !getEpic) {
		// Update Message ...
		if (!getTimewin) timeWin = 50; if (!getFreq) Fcutoff = 1; 
		if (!getEpic){
			epbazmax = 360;
			epbazmin = 0;
			epicmax = 180;
			epicmin = 0;
		}
	}
	
	// These 3 Options are Data selection operators + Debug
	if (!getHeaderTag || !getLqtSwch) {
		timeTag = 1; lqtRotArg = 1;
	}
	 
	
	
	// Administrative info Here ... 
	if (err){
		
		usage_and_exit();
			
	} else {
		cout << "\n\n\n***ScatteredWaveAnalysis(RF - Sequential HK Analysis)*** Created by Tolulope Olugboji. \n************************** Copyright 2013 Yale University. All rights reserved. \n**Option: Code does sequential grid search for Upper crustal structure using Radial RFs only!: " << endl;

		
	}
	
	
	if (getInterpolate == true) {
		int ROUTERCODE = HKSTCKMULTIPLEINTERPOLATE;
		cout << "!****Working on the demo for interpolation and findLABDepth: \n usage: -I2 interpolates with twice the no. of points." << endl;
		
		GrdSrchRouter singleHKStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, timeTag, searchParamsfn, ROUTERCODE, epbazmin, epbazmax, epicmin, epicmax, nDecim, getVerbose);
	}else{
		//	int ROUTERCODE = HKSTCKSINGLE;
		int ROUTERCODE = HKSTCKMULTIPLE;
		// There are some modifications from normal RFrouter. No binsize, or migration params
		// No need to increment bins also ...
		GrdSrchRouter singleHKStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, timeTag, searchParamsfn, ROUTERCODE, epbazmin, epbazmax, epicmin, epicmax, getVerbose);
		
	}
	
	if (getFreqMigrateFlag) {
		int ROUTERCODE = HKSTCKRECURSEFREQ;
		// There are some modifications from normal RFrouter. No binsize, or migration params
		// No need to increment bins also ...
		// **** Here I call routing code that does single- recursive layer stacking ...
		GrdSrchRouter singleHKStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, timeTag, searchParamsfn, ROUTERCODE, epbazmin, epbazmax, epicmin, epicmax, getVerbose);
	}
	
	

	
	return 0;
	
}

