/** @mainpage RecFunc.cpp
 *  
 *
 *  @author Tolulope Olugboji (tolumorayo@gmail.com)  
 *  @date Nov. 11, 2010.
 *	@version   3.1    Nov.  13, 2012 [Harmonic Stacking, ParkOlugbojiSep.key]
 *
 *  @section LICENSE 
 *       Copyright 2010 Yale University. All rights reserved. 
 *  
 *  @section DESCRIPTION
 *			Version 2.0      closed on April 14, 2011   [One option still open, retire?]
 *          Version 2.1      commenced on June  17, 2011

 *	Usage: RecFunc.cpp -LEventlist -Ffmax  -BMaxBinsize -Abazmin/bazmax/bazinc 
 *	                   -Ebazmax/bazmin/epicmin/epicmax/epicinc -Ttimewin -HheaderTag -Vverbose 
 *  
 *  Requested Update (April 26, 2011): -D option
 *                      T[his option specifies which epicentral distance range is 
 *					    acceptable for computation. Useful for discriminating records that are not
 *						teleseismic. Ongoing!]
 *  Requested Update (June 17, 2011):  -M migration module or option...
 *                                      still Thinking of How that will be implemented??
 *                                     -S statistics module, after CCP stacking.
 *                                     
 *  Requested Utility (June 17, 2011): A plotter class external to my main module, useful for 
 *                                      debugging. Make extensible & toggle capable
 *
 *  Requested Update (June 17, 2011): Write SvRecFunc.cpp : 
 *                                     S receiver function After migration & Stats module
 *  Implementing Update
 *  
 *
 *  Files Output by RecFunc:
 *     [V]     1. RFs by Azim:  Azim_Rad.xyz & Azim_Trans.xyz
 *     [V]     2. RFs by Epic:  Epic_Rad.xyz & Epic_Trans.xyz
 *     [V]     3. Spectrum:     Rad_Spec.xyz,  Tran_Spec.xyz & Noise_Spec.xyz
 *     [X]     4. Coherence:    Coher_Rad.xyz & Coher_Trans.xyz
 *     [V]     5. Event Stat:   FileNme_Event.txt
 *     [V]     6. Log File:     FileNme.log
 */


#include "error.h"
#include <cstdio>
#include <cctype>
#include <iostream>
#include <stdlib.h>
#include "RFrouter.h"

using namespace std;

void usage_and_exit(){
	cout << "\n\n\n***ScatteredWaveAnalysis(Radial-Transverse Ps Receiver Functions) Version 2.0. *** \nAuthor: Tolulope Olugboji. \nCopyright: 2013 Yale University. All rights reserved.  \nVersion: Original Nov. 11, 2010. Updated April 8, 2014. [Fix slowness units]\nUpdated April 25, 2014 [add horizontal rotation]\nUpdating April 20, 2015 [adding jacknife to all stacks with bin stat output] \n  *********************************************************** \nFunctionalities in Version 2.0: \n 1.0 Radial[R] Transverse[T] Ps receiver functions stacked by Azimuth and Epicentral Distance \n 2.0 Harmonic Decomposition of R-T Ps Receiver Functions (e.g. Olugboji et al., Gcubed 2016a,b) \n 3.0 Moving Window Migration for: \n \t 1. Direct Ps Phase \n \t 2. Reverberated PPPMS Phase \n \t 3. Reverberated PPSMS Phase \n 4.0 LQT Rotation  \n *********************************************************** \n  \n \n";
	cout << "Usage: RecFuncHarmonic  [options] \n \n";
	
	cout << "*** File Handling Options - Read, Write, Windowing and Frequency Tags ********************************\n";
	cout << "-L(inFname) \t File name of text file holding the [L]ist of SAC records. The leading [.RTZ] attached internally \n";
	cout << "-O(outFname) \t File name where the [O]utput xyz files should be stored. \n";
	cout << "-T(nSec) \t [T]ime window in seconds. Default 50 secs.\n";
	cout << "-T(-nSec) \t If negative [T]ime window used then only noise is used for deconv: OBS check! \n";
	cout << "-F(fMax) \t Maximum [F]requency for computing the Receiver Functions. Default 1 Hz. \n\n";
	
	cout << "*** Stacking Mode Options - Azim/Epic Stack or Harmonic Stack ********************************\n [-S(sType=1/nBoot)] [-S(sType=0/nBoot) -B(binsze) -A(bazmin/bazmax/bazinc) [-E(bazmin/bazmax/epicmin/epicmax/epicinc)]] \n";
	cout << "-S(sType/nBoot) \t select stack mode:\n\t\t\t sType = 0: Azimuth or Epic Stack. \n\t\t\t sType = 1: Harmonic Stack; nBoot: no. of times for bootstrap resampling\n";
	cout << "-A(bazmin/bazmax/bazinc) \t Not used if [-Sstype=1]. Parameters for restricting stack to earthquakes in particular azimuth sector.\n\t\t\t bazmin: right back azimuth \n\t\t\t bazmax: left back azimuth \n\t\t\t bazinc: increment for back azimuth binWidth  \n";
	cout << "-E(bazmin/bazmax/epicmin/epicmax/epicinc) \t  Not used if [-Sstype=1]. Parameters for restricting stack to earthquakes in particular location.\n\t\t\t bazmax: left back azimuth \n\t\t\t bazmin: right back azimuth \n\t\t\t epicmin: minimum distance in deg. \n\t\t\t epicmax: maximum distance in degrees \n";
	cout << "-B(binsze) \t minimum no. of records allowed in a stacking bin  \n\n";
	
	cout << "*** Migration Mode Options - Phase type, Target Depth, Migration algorithm ********************************\n[-M(mType=0,p/tDepth)] [-m(mType=0,1/tDepth)]  -I(velFname) \n";
	cout << "-M(mType/tDepth) \t Ps phase select migration mode:\n\t\t\t mType = 0: Moving Window Migration. \n\t\t\t mType = p: !!Print File With timing information. Skips receiver function computation! \n\t\t\t tDepth(km): conversion depth in km \n \t\t\t ***Output will print 2 files time migration and depth(_Depth...xyz) migration [Azim and Epic Only] \n";
	cout << "-m(mType/tDepth) \t Reverberated phase select migration mode:\n\t\t\t mType = 0: PpSms phase \n\t\t\t mType = 1: PpPms phase \n\t\t\t tDepth(km): conversion depth in km \n";
	cout << "-I(velFname) \t File name for nLayer velocity model. format similar to JeffPark anirecsynth.f \n \n";
	
	cout << "***Other Options - QualityCheck Discrimination, Rotation, and Verbosity ********************************\n";
	cout << "-H \t Use [H]eader tag look up. IF time tag is NOT set, skip record.  \n";
	cout << "-R \t Set 1 if LQT [R]otation is to be done.  \n";
	cout << "-Rh/horAngle \t Work in progress - horizontal rotation.  \n";
	cout << "-V \t Increase verbosity \n" << endl;
	
	exit(1);
	
	/*
	error(" Created by Tolulope Olugboji. \n Copyright 2010 Yale University. All rights reserved. Version 2.0 Migration module added, including reverberated phases. \n Updated July 1, 2013. \n Usage RecFunc: -Larraylist -Ffmax -H1 [-Bbin_size] [-Abazmin/bazmax/bazinc] [-Ebazmin/bazmax/epicmin/epicmax/epicinc] [-S[StackMode/StackModeArg]] \n [-M[MigrationMode/MigrationModeArg]] [-R[lqtRotArg]] [-Vverbose] [-m(reverbPhase/migDepthKm)]")
	 */
	
}

int main(int argc, char** argv){
	
	// Variables to hold option flags. They determine what modules to be run in RFrouter
	bool getHelp=false;
	bool getArray=false, getFreq=false, getAzim=false, getEpic=false, err=false;
	bool getBinsz=false, getTimewin=false, getOut = false, getVerbose = false;
	bool getHeaderTag=false, getLqtSwch = false, getMigrationModel = false;
	bool getMigrationMode=false, getMigReverbMode=false, getStackMode=false;
	bool getMigrationPrint=false;	// If set... Skip RF computation and print timing information only!

	
	// Variables to support options: additional option arguments. used by RFrouter ...
	int binsize=0; 
	int timeTag =1;								// Use header t1 by default: timeTag=1
    int migrationModeFlag, migrationModeArg, stackModeFlag, stackModeArg;
	char migrationPrint;
	int migDepthKm;							    // Depth for reverb phases (Expected Moho)
	int reverbPhaseFlag;							// Choose between PpSms or PpPms
	float bazmin, bazmax, bazinc, epicmin, epicmax;
	float epbazmin, epbazmax;
	float epicinc, Fcutoff, timeWin;
	float lqtRotArg = -1.0; char rotationModeChar; float horRotAng=0.0;   // rotation options
	char *eventsfn, *filename, *outfn, *velocityfn;
	
	//Retrieve arguments and ensure array file well parsed.
	for (int i=1; i<argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
				case 'L':
					eventsfn = &argv[i][2];
					getArray = true;
					break;
				case 'T':
					getTimewin = true;
					timeWin = atof(&argv[i][2]);
					break;
				case 'F':
					Fcutoff = atof(&argv[i][2]);
					getFreq = true;
					break;
				case 'B':
					// Work on passing binsize to the TraceStack routine
					binsize = atoi(&argv[i][2]);
					getBinsz = true;
					break;
				case 'A':
					if (sscanf(&argv[i][2], "%f/%f/%f", 
							   &bazmin, &bazmax, &bazinc) == 3) {
						getAzim = true;
					}else {
						error("Invalid -A option. Using default Azimuth setting");
					}
					break;
				case 'E':
					if (sscanf(&argv[i][2], "%f/%f/%f/%f/%f", 
							  &epbazmin, &epbazmax, &epicmin, &epicmax, &epicinc) == 5) {
						getEpic = true;
					}else {
						error("Invalid -E option. Using default Epic. setting");
					}
					break;
				case 'O':
					outfn = &argv[i][2];
					getOut = true;
					break;
				case 'H':
					getHeaderTag = true;
					timeTag = atoi(&argv[i][2]);
					break;
                case 'M':
                    if (sscanf(&argv[i][2], "%d/%d", &migrationModeFlag, &migrationModeArg) == 2) {
                        getMigrationMode = true;
                    }else if (sscanf(&argv[i][2], "%c/%d", &migrationPrint, &migrationModeArg) == 2){
						getMigrationPrint = true;
					}
					else{
                        error("Invalid -M option. See usage for details!");
                    }
                    break;
				case 'm':
                    if (sscanf(&argv[i][2], "%d/%d", &reverbPhaseFlag, &migDepthKm) == 2) {
                        getMigReverbMode = true;
                    }else{
                        error("Invalid -m option. Usage: [-m(reverbPhaseFlag/migDepthKm)] Code migrates to moho (or depth of reverberation) for reverberating phases: 1 for PpSms 2 for PpPms ...");
						//  Could use a more suitable flag that does not interfere with the other migration flag.. best hack so far though.	
                    }
                    break;
				case 'S':
                    if (sscanf(&argv[i][2], "%d/%d", &stackModeFlag, &stackModeArg) == 2) {
                        getStackMode = true;
                    }else{
                        error("Invalid -S option. Using default Stack mode: Azimuthal & Epicentral");
                    }
                    break;
				case 'R':
					// Working on extending this option for two cases - both horizontal and lqt rotation
					// Multiple grid rotations?
					// I also allow for backward compatibility so code doesn't break
					if (sscanf(&argv[i][2], "%c/%f", &rotationModeChar, &horRotAng) == 2) {
						cout << "New Hack: run Horizontal rotation " << horRotAng << endl;
					}else{
						cout << "Backward compatibility: run LQT rotation" << endl;
						
						getLqtSwch = true;
						lqtRotArg = atof(&argv[i][2]);
					}
					break;
				case 'V':
					getVerbose = true;
					break;
				case 'I':
					velocityfn = &argv[i][2];
					getMigrationModel = true;
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
	//End Arguments processing ..
	// This is a hack to allow for testing argument processing
	// Uncomment line below to allow code to proceed with normal execution
	// usage_and_exit();
	
	// Check options input or exit
	if (!getArray || !getFreq || !getTimewin || getHelp){
		
		if (getMigrationPrint) {
			
		}else{
			err=true;
			usage_and_exit();
		}

	}
	
	
	// Set default time header field = T0
	if ( !getHeaderTag) timeTag = 0;
	
	// If LQT rotation requested then default arg of -1 forlqtRotArg is replaced
	// This argument is used in sacread to determine if LQT rotation is required or not.
	if (getLqtSwch) {
		if (getVerbose) cout << "\n\n\n***ScatteredWaveAnalysis(RF)*** Created by Tolulope Olugboji. \n************************** Copyright 2010 Yale University. All rights reserved. \n**Option: Switch for LQT rotation activated with incidence angle bias: " << lqtRotArg << endl;
	}
    
	// ********* Do Arguments Checking for Migration Mode To Use. Ensuring that the Velocity File is Well Constructed
	if (getMigrationModel && getMigrationMode){
		if (getVerbose) cout << "**Option: Direct Phase Migration with model: " 
			<< velocityfn << endl;
		
	}else if (getMigReverbMode && getMigrationModel) {
		if (getVerbose) cout << "**Option: Reverberated Phase Migration with model: " 
			<< velocityfn << endl;
	}else if (getMigrationPrint && getMigrationModel){
		 if (getVerbose) cout << "***Option: Skip RF computation!! Print Timing Only!!!  " << endl;	
	}
	else if (  (!getMigrationModel) && (getMigrationMode || getMigReverbMode || getMigrationPrint) ){
		error(" Either Migration Mode (Reverb or Direct of Print) or Migraton Model Not Set: See usage! ");
	}
	
	// ******* Do Arguments Checking for Stacking Mode.. Otherwise Use Defaults
	if (getStackMode){
		
		// Report stack mode so user can know what is goign on.
		switch (stackModeFlag) {
			case 0:
				if (getVerbose) cout << "**Option: Azimuth-EpicentralDist Stack Selected" << endl;
				break;
			case 1:
				if (getVerbose) cout << "**Option: Harmonic Stack Selected" << endl;
				break;
			case 2:
			    if (getVerbose) cout << "**Option: Single-Prolate (P=1) ECC mode selected" << endl;
				break;
			default:
				if (getVerbose) cout << "**Option: Not recognized..." << endl;
				break;
		}
		
		
	}else if (!getStackMode && !getMigrationPrint) {
		error("Include the type of Stacking -S(sType/nBoot): \n [0. Azimuth-Epicentral] Or [1. Harmonic Stack]");
	}
	
	// Check Mode Set and Call The RF router with appropriate Parameters.
	if ( !getMigrationMode && !getMigReverbMode && !getMigrationPrint) {
		cout << "No Migration ...." << endl;
		switch (stackModeFlag) {
			case 0:
			{// Call RF router with Epicentral-Azimuth Stack - NULL for velocity model
				// added horizontal rotation angle to constructor
				int ROUTERCODE = SIMPLEAZIMEPIC;
				
				if (getAzim) {
					RFrouter noMigrationEpicAzimStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg, 0, NULL,ROUTERCODE, bazmin, bazmax, bazinc, binsize, getVerbose);
				}
				if (getEpic){
					RFrouter noMigrationEpicAzimStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg, 0, NULL, ROUTERCODE, epbazmin, epbazmax, epicmin, epicmax, epicinc, binsize, getVerbose);
					
				} else if (!getAzim && !getEpic) {
					error("Unable to determine Azimuth or Epicentral Stack set -A or -E");
				}

			}
				break;
			case 1:
			{// Call RF router with Harmonic Stack - NULL for velocity model
				// add horizontal rotation argument to constructur...
				int ROUTERCODE = SIMPLEHARMONIC;
				RFrouter stackRFs(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg,0, NULL, ROUTERCODE, getVerbose);
			}
				break;	
			case 2:
			{// Call RF router with Epicentral-Azimuth Stack - NULL for velocity model
				// force P = 1 and K = 1 ... see if this works ...
				int ROUTERCODE = ECCAZIMEPIC;
				
				if (getAzim) {
					RFrouter noMigrationEpicAzimStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg, 0, NULL,ROUTERCODE, bazmin, bazmax, bazinc, binsize, getVerbose);
				}
				if (getEpic){
					RFrouter noMigrationEpicAzimStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg, 0, NULL, ROUTERCODE, epbazmin, epbazmax, epicmin, epicmax, epicinc, binsize, getVerbose);
					
				} else if (!getAzim && !getEpic) {
					error("Routing for ECC: Unable to determine Azimuth or Epicentral Stack set -A or -E");
				}
			
			}
			    break;			
			default:
				break;
		}
	} else if (getMigrationMode && getMigrationModel){
		switch (stackModeFlag) {
			case 0:
			{	// Call RF router with Crude Stack + Migration Routine. [Only MWM Implemented]
				// Use MigrationMode to call RF router...
				cout << "[Azimuth Epicentral Stack Called with MWM migration]" << endl;
				int ROUTERCODE = MWMAZIEPIC;
				
				if (getAzim) {
					RFrouter noMigrationEpicAzimStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg, migrationModeArg, velocityfn, ROUTERCODE, bazmin, bazmax, bazinc, binsize, getVerbose);
				}
				if (getEpic){
					RFrouter noMigrationEpicAzimStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng,  timeTag, stackModeArg, migrationModeArg, velocityfn, ROUTERCODE, epbazmin, epbazmax, epicmin, epicmax, epicinc, binsize, getVerbose);
					
				} else if (!getAzim && !getEpic) {
					error("Unable to determine Azimuth or Epicentral Stack set -A or -E");
				}
			}
				break;
			case 1:
				// Call RF router with Harmonic Stack + Migration Routine [Only MWM implemented]
				cout << "[Harmonic Stack Called With MWM migration]" << endl;
			{// Call RF router with Harmonic Stack - NULL for velocity model
				int ROUTERCODE = MWMHARMONIC;
				RFrouter stackRFs(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg,migrationModeArg, velocityfn, ROUTERCODE, getVerbose);
			}
				break;
			default:
				break;
		}
	}else if (getMigReverbMode && getMigrationModel){
		// Include the Harmonic Reverberation Migration...
		switch (stackModeFlag) {
			case 0:
			{
				cout << "[Azim-Epic Called With MWM migration For MohoReverb Phase" << endl;
				switch (reverbPhaseFlag) {
					case 0:
					{
						cout << "**Option: Selected PpSms Phase Migration" << endl ;
						int ROUTERCODE = MOHOREVERBPSMS;
						
						if (getAzim) {
							RFrouter MigrationEpicAzimStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg, migDepthKm, velocityfn, ROUTERCODE, bazmin, bazmax, bazinc, binsize, getVerbose);
						}
						if (getEpic){
							RFrouter MigrationEpicAzimStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg, migDepthKm, velocityfn, ROUTERCODE, epbazmin, epbazmax, epicmin, epicmax, epicinc, binsize, getVerbose);
							
						} else if (!getAzim && !getEpic) {
							error("Unable to determine Azimuth or Epicentral Stack set -A or -E");
						}
						
						break;
					}
					case 1:
					{
						cout << "**Option: Selected PpPms Phase Migration" << endl ;
						int ROUTERCODE = MOHOREVERBPPMS;
						
						if (getAzim) {
							RFrouter MigrationEpicAzimStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg, migDepthKm, velocityfn, ROUTERCODE, bazmin, bazmax, bazinc, binsize, getVerbose);
						}
						if (getEpic){
							RFrouter MigrationEpicAzimStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg, migDepthKm, velocityfn, ROUTERCODE, epbazmin, epbazmax, epicmin, epicmax, epicinc, binsize, getVerbose);
							
						} else if (!getAzim && !getEpic) {
							error("Unable to determine Azimuth or Epicentral Stack set -A or -E");
						}
						
						break;
					}	
					default:
						break;
				}
			}
				break;
			case 1:
				// Call RF router with Harmonic Stack + Migration Routine [Only MWM implemented]
				cout << "[Harmonic Stack Called With MWM migration For MohoReverb Phase" << endl;
			{
				switch (reverbPhaseFlag) {
					case 0:
					{
						int ROUTERCODE = HARMONICREVERBPSMS;
						RFrouter stackRFs(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg,migDepthKm, velocityfn, ROUTERCODE, getVerbose);
						break;
					}
					case 1:
					{
						int ROUTERCODE = HARMONICREVERBPPMS;
						RFrouter stackRFs(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg,migDepthKm, velocityfn, ROUTERCODE, getVerbose);
						break;
					}	
					default:
						break;
				}
			}
			
				break;
			default:
				break;
		}
		
		
	}else if (getMigrationPrint && getMigrationModel){
		cout << "Skipping RF computation! Computing average timing for all Phases in the bin stack using ray parameter for each record. file in: " << outfn << "_Epic_PhaseTime.txt" << endl;
		
		if (getAzim || !getEpic) {
			error("Printing option can be used only with -E option to view moveout of phase with epicentral distance!");
		}else{
			// Pass the routing code. Implemented completely in the RFrouter class
			
			int ROUTERCODE = PRINTPHASETIME;
			RFrouter printMigrationEpicStack(eventsfn, outfn, timeWin, Fcutoff, lqtRotArg, horRotAng, timeTag, stackModeArg, migrationModeArg, velocityfn, ROUTERCODE, epbazmin, epbazmax, epicmin, epicmax, epicinc, binsize, getVerbose);
		}
		
		
		
	}
	
	
	
	return 0;
	
}


