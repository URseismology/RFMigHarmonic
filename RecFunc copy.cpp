/*
 *  RecFunc.cpp
 *  
 *
 *  Created by Tolulope Olugboji on Nov. 11, 2010.
 *			   Version 2.0      closed on April 14, 2011   [One option still open, retire?]
 *             Version 2.1      commenced on June  17, 2011
 *             Version 3.1      commenced on Nov.  13, 2012 [Harmonic Stacking, ParkOlugbojiSep.key]
 *
 *  Copyright 2010 Yale University. All rights reserved. *
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
#include "RecordList.cpp"
#include "TraceStack.h"
#include "HarmonicStack.h"
#include "MigrationParams.h"
#include "MTCDriver.h"


#define  STACKAZIM 0
#define  STACKEPIC 1

#define  CMPSHARMONIC 0
#define  EXPANSIONHARMONIC 1

int main(int argc, char** argv){
	bool getArray=false, getFreq=false, getAzim=false, getEpic=false, err=false;
	bool getBinsz=false, getTimewin=false, getOut = false, getVerbose = false;
	bool getHeaderTag=false, getMode=false, getLqtSwch = false, getMigrationModel = false;
	// Variables to hold length of noise and Data following onset of pulse
	int binsize=0; 
	int timeTag =1;									// Use header t1 by default: timeTag=1
    int modeFlag, modeArg;
	float bazmin, bazmax, bazinc, epicmin, epicmax;
	float epbazmin, epbazmax;
	float epicinc, Fcutoff, timeWin;
	float lqtRotArg = -1.0;
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
                    if (sscanf(&argv[i][2], "%d/%d", &modeFlag, &modeArg) == 2) {
                        getMode = true;
                    }else{
                        error("Invalid -M option. Using default Stack mode, Azim, Epicentral");
                    }
                    break; 
				case 'R':
					getLqtSwch = true;
					lqtRotArg = atof(&argv[i][2]);
					break;
				case 'V':
					getVerbose = true;
					break;
				case 'I':
					velocityfn = &argv[i][2];
					getMigrationModel = true;
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
	
	// Check options input or exit
	if (!getArray || !getFreq || !getTimewin ) err=true;
	if (err) error("Version 1.0 Harmonic Regression module added. Updated May 22, 2013. \n Usage RecFunc: -Larraylist -Ffmax -H1 [-Bbin_size] [-Abazmin/bazmax/bazinc] [-Eepicmin/epicmax/epicinc] [-M[modeArg]] [-R[lqtRotArg]] [-Vverbose]");
	
	// Set default time header field = T1
	if ( !getHeaderTag) timeTag = 1;
	
	// If LQT rotation requested then default arg of -1 forlqtRotArg is replaced
	// This argument is used in sacread to determine if LQT rotation is required or not.
	if (getLqtSwch) {
		if (getVerbose) cout << "Switch for LQT rotation activated with incidence angle bias: " << lqtRotArg << endl;
	}
    
	// Add Section to read isotropic model for migration: moving window etc.
	if (getMigrationModel){
		if (getVerbose) cout << "Switch to read migration model set: " << velocityfn << endl; 
		MigrationParams velModel(velocityfn);
		velModel.getTimeDelay(4, 0);
		
		if (getVerbose) cout << "No of Layers + Halfspace is = " << velModel.timeDelayLayer.size() << endl;
		return 0;
	}
    
    /**********+++++++++++++++++*****************************+++++++++++++-----------------------------------||||||||||||||||||||||||||||||||||||||||||||*/
    //** My November 13, 2012 Update. I Select Stack Logic Based On The Mode Selector.
	
	
    if (getMode && modeFlag == 1) {
        // ModeFlag == 1 is Harmonic Stacking and That's what this section does. Calls The classes relevant to this section
        if (getVerbose) cout << "Harmonic Stacking Module Using mode: " << modeFlag << " No of times to compute bootstrap: " 
        << modeArg << endl;
        
        /*Scan through records and stack through a harmonic regression: by azimuth.
         My code should have HarmonicStack.regress()?
         Last two numbers represent time windows in the final stacked sections: pre & post
         */
		
		
        // Do record analysis. Administrative details to build records.
           // Logistics of Parsing good vs. bad records and forming spectrum...
            //Build Filenames, Record Counts, Baz, Epic, ...
            RecordList allrecords(eventsfn, timeWin, Fcutoff);
            
            float tEventTrace = timeWin, tNoiseTrace;
            Int nNoiseTrace, nEventTrace, nMax;
            
            //Pick out first record for administrative details.
            std::string first;
            first.assign(allrecords.recordname[0]);
            Sacread record(first, timeTag, lqtRotArg);
            
            float dT = record.deltaT;
            
            nEventTrace = (tEventTrace) / (record.deltaT);
            
            
            /* Scan through Each SAC file and reconstitute record list based on
             header entry. If entry not set, then skip record and put name in new record list
             */
            if (getVerbose) cout << "Checking Record list for tagged headers" << 
                "Total no of records: " << allrecords.nTotrec << endl;
            for (int irec =0; irec < allrecords.nTotrec; irec++) {
                record = Sacread(allrecords.recordname[irec], timeTag, lqtRotArg);
                
                tNoiseTrace = record.timeTag - record.timeStart - 3; // bias it by 3 sec
                nNoiseTrace = (tNoiseTrace) / (record.deltaT);
                
                /* Skip record with untagged field. If Tagged, Check to see that there's enough data for time window.
                 VERY IMPORTANT!!! USE record.timeTag and deltaT to do this.*/
                if (record.timeTag > 0) {
                    if( getVerbose) cout << "Record: " << irec+1 << " has the header no:  "<< timeTag << " Tagged " <<  "PhaseName: " << record.phaseName << endl;
                    
                    bool enoughDataEvent = (record.DataCmps.ncols() - ((record.timeTag - record.timeStart) / record.deltaT) ) >= nEventTrace;
                    bool enoughDataNoise = ( (record.timeTag - record.timeStart) / record.deltaT ) >= nNoiseTrace ;
                    
                    
                    if (enoughDataNoise && enoughDataEvent ) {
                        allrecords.passNew(irec);                  // Pass Record only if there's enough data, Obviously!!
                    }
                    if ( !enoughDataNoise || !enoughDataEvent) {
                        if (getVerbose) cout << "Buggy data. Here is the cause of the Segmentation fault" << endl;
                    }
                    
                }else {
                    // Update here if 
                    if( getVerbose) cout << "Header no:  "<< timeTag << " is not tagged for record: " << irec+1 << endl;
                }
                
                
            }
            
            
            
            
            //Pick out first GOOD record for administrative details.
            first.assign(allrecords.goodRecordname[0]);
            record = Sacread(first, timeTag, lqtRotArg);
            
            dT = record.deltaT;
            
            
            //Mat3DDoub noiseTrace, postEventTrace, CoherTrace;
            //Mat3DCmplx RFTrace;
            
            tNoiseTrace = record.timeTag - record.timeStart - 3; // bias it by 3 sec
            nEventTrace = (tEventTrace) / (record.deltaT);
            nNoiseTrace = (tNoiseTrace) / (record.deltaT);
            nMax = MAX(nEventTrace, nNoiseTrace);
            
            /* Define pad length based on event number, to power of 2
             Note here I pick maximum length of either the event or noise window
             */
            
            int expo = int(log2(nMax-1)) + 1;
            int nPad = pow(2.0, expo);
            
            
            //Define 3 dimensional array for all components in Event Array
            
            Mat3DDoub noiseTrace(allrecords.nGoodrec, 3, nPad); //-> fucked up! because?
            Mat3DDoub postEventTrace(allrecords.nGoodrec, 3, nPad);
            Mat3DCmplx RFTrace(allrecords.nGoodrec, 2, nPad); 
            Mat3DDoub deltaRFTrace(allrecords.nGoodrec, 2, nPad); 
            
            // Set all the Five Headers useful for Statistics (& Selective Parsing?) Set file to Store.
            allrecords.setAzim(); allrecords.setEpic(); allrecords.setEvLong();
            allrecords.setEvLat(); allrecords.setEvMag();
            string evStatfn(outfn);
            evStatfn.append("_Event.txt");
            ofstream evStat(evStatfn.c_str()); 
            
            /* Prep 3 files for Spectrum: Noise and Event spectrum
            string specNames[] = { "Vert_Spec.xyz", "Rad_Spec.xyz", "Tran_Spec.xyz", 
                "Noise_Spec.xyz"
            }; 
            ofstream specFNmes[4];
            for (int j = 0; j < 4 ; j++) {
                specFNmes[j].open(specNames[j].c_str());
            }
            //Prep 2 Files for Coherence Plot, See frequency behaviour
            string coherRadfn[] = { "Coher_Rad.xyz", "Coher_Tran.xyz" };
            ofstream coherFNmes[2];
            for (int j = 0; j < 2 ; j++) {
                coherFNmes[j].open(coherRadfn[j].c_str());
            }
            */
            
            if (getVerbose) cout << "Total Good Records: " << allrecords.nGoodrec << endl ;
            /*Create RF traces for only the Good records and print each individual trace?
             Try adapting WriteRFdata to be used outside the Constructor... and
             to write data to file outRadial.grid & outTrans.grid
             */
            for (int irec =0; irec < allrecords.nGoodrec; irec++) {
                
                record = Sacread(allrecords.goodRecordname[irec], timeTag, lqtRotArg);
                tNoiseTrace = record.timeTag - record.timeStart - 3; // bias it by 3 sec
                nNoiseTrace = (tNoiseTrace) / (record.deltaT);
                
                /* 
                 **** URGENT UPDATE REQUIRED!!! This construction fails if record length following header is too small,
                 ****							Make check during record pass
                 Construct noiseTrace & postEventTrace
                 I do this for all the components: vertical, radial & transpose.
                 ? Do I need routine for clarity in code Prose?
                 */
                
                for (int cmp = 0; cmp < 3; cmp++) {
                    
                    for (int iter=0; iter < nPad; iter++) {
                        if (iter < nNoiseTrace) {
                            noiseTrace[irec][cmp][iter] = record.DataCmps[cmp][nNoiseTrace - iter];
                        } else {
                            noiseTrace[irec][cmp][iter] = 0.0;
                        }
                        
                    }
                    
                    
                    int posEventNxt, posEventStrt = nNoiseTrace;
                    for (int iter=0; iter < nPad; iter++) {
                        posEventNxt = posEventStrt + iter;
                        if ( iter < nEventTrace) {
                            postEventTrace[irec][cmp][iter] = record.DataCmps[cmp][posEventNxt];
                        } else {
                            postEventTrace[irec][cmp][iter] = 0.0;
                        }
                    }
                    
                }
                
                
                allrecords.Azim[irec] = record.recBaz;
                allrecords.Epic[irec] = record.recEpic;
                allrecords.EvLong[irec] = record.evLon;
                allrecords.EvLat[irec] = record.evLat;
                allrecords.EvMag[irec] = record.evMag;
                
                // Event Statistics file is updated here, Five important headers.
                evStat << record.recBaz << "\t" << record.evLon <<  "\t" << record.evLat <<  "\t" << record.evMag
                << "\t" << record.recEpic << endl;
                
                //cout << "File Name: " << allrecords.recordname[irec]
                //<< "Azimuth: " << allrecords.Azim[irec]
                //<< "Epicenter " << allrecords.Epic[irec] << endl;
                
                /*Pass Record to MTC driver, computes Spectral Estimates*/
                /*Args: Noise, Event, p, k, fmax */
                MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 2.5, 3, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
                // Code above doesn't work .. Check type
                
                if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << endl;
                MTCVals.addRF(RFTrace, irec);
                MTCVals.addCoher(deltaRFTrace, irec);
                
                // Dump Spectrum & Coherence For debug purposes
                /*
                for (int j = 0; j < 4 ; j++) {
                    MTCVals.updateSpectrum(specFNmes[j], j, irec+1);
                }
                
                for (int j = 0; j < 2 ; j++) {
                    MTCVals.updateCoher(coherFNmes[j], j, irec+1);
                }
                 */
                
            }
            
            
            /*Scan through records and stack through a harmonic regression: by azimuth.
             My code should have HarmonicStack.regress()?
             Last two numbers represent time windows in the final stacked sections: pre & post
             */
            string logfn(outfn);
            logfn.append(".log");
            ofstream logFile(logfn.c_str());   
	
        
        /*Scan through records and stack through a harmonic regression: by azimuth.
         My code should have HarmonicStack.regress()?
         Last two numbers represent time windows in the final stacked sections: pre & post
         */
        
        
        int nBoot = modeArg;
        
		int flagFile = CMPSHARMONIC;
		
		int nRecs = allrecords.nGoodrec, ncmps = 2, nFreq = nPad;
		HarmonicStack testHStack(nPad, ncmps, nRecs, allrecords.Azim, allrecords.Epic, Fcutoff, dT, 10.0,  30.0, nBoot);
        testHStack.regressNBootTimes(RFTrace, deltaRFTrace);
		testHStack.print(outfn, flagFile);
		logFile << testHStack.logReport;
        
        
        
        
    
	} 
	else if ( getMode && modeFlag == 0){
    
    /* This section was used for the crude stack initially. The mode selector can be extended to capture:
       CrudeStack, HarmonicUnmodelledStack, etc. For now this else entry just selects The CrudeStack.
     */
     
		if (getVerbose) cout << "Running  Module for Crude Azimuthal and Epicentral Stacks." << endl;
			
	//Build Filenames, Record Counts, Baz, Epic, ...
	RecordList allrecords(eventsfn, timeWin, Fcutoff);
	
	float tEventTrace = timeWin, tNoiseTrace;
	Int nNoiseTrace, nEventTrace, nMax;
	
	//Pick out first record for administrative details.
	std::string first;
	first.assign(allrecords.recordname[0]);
	Sacread record(first, timeTag, lqtRotArg);
	
	float dT = record.deltaT;
	
	nEventTrace = (tEventTrace) / (record.deltaT);
	
	
	/* Scan through Each SAC file and reconstitute record list based on
	 header entry. If entry not set, then skip record and put name in new record list
	 */
	if (getVerbose) cout << "Checking Record list for tagged headers" << 
		"Total no of records: " << allrecords.nTotrec << endl;
	for (int irec =0; irec < allrecords.nTotrec; irec++) {
		record = Sacread(allrecords.recordname[irec], timeTag, lqtRotArg);
		
		tNoiseTrace = record.timeTag - record.timeStart - 3; // bias it by 3 sec
		nNoiseTrace = (tNoiseTrace) / (record.deltaT);
		
		/* Skip record with untagged field. If Tagged, Check to see that there's enough data for time window.
		   VERY IMPORTANT!!! USE record.timeTag and deltaT to do this.*/
		if (record.timeTag > 0) {
			if( getVerbose) cout << "Record: " << irec+1 << " has the header no:  "<< timeTag << " Tagged" << endl;
			
			bool enoughDataEvent = (record.DataCmps.ncols() - ((record.timeTag - record.timeStart) / record.deltaT) ) >= nEventTrace;
			bool enoughDataNoise = ( (record.timeTag - record.timeStart) / record.deltaT ) >= nNoiseTrace ;
			
		
			if (enoughDataNoise && enoughDataEvent ) {
				allrecords.passNew(irec);                  // Pass Record only if there's enough data, Obviously!!
			}
			if ( !enoughDataNoise || !enoughDataEvent) {
				if (getVerbose) cout << "Buggy data. Here is the cause of the Segmentation fault" << endl;
			}
	
		}else {
			// Update here if 
			if( getVerbose) cout << "Header no:  "<< timeTag << " is not tagged for record: " << irec+1 << endl;
		}

		
	}
	
	
	
	
	//Pick out first GOOD record for administrative details.
	first.assign(allrecords.goodRecordname[0]);
	record = Sacread(first, timeTag, lqtRotArg);
	
	dT = record.deltaT;
	
	
	//Mat3DDoub noiseTrace, postEventTrace, CoherTrace;
	//Mat3DCmplx RFTrace;
	
	tNoiseTrace = record.timeTag - record.timeStart - 3; // bias it by 3 sec
	nEventTrace = (tEventTrace) / (record.deltaT);
	nNoiseTrace = (tNoiseTrace) / (record.deltaT);
	nMax = MAX(nEventTrace, nNoiseTrace);
	
	/* Define pad length based on event number, to power of 2
	 Note here I pick maximum length of either the event or noise window
	 */
	
	int expo = int(log2(nMax-1)) + 1;
	int nPad = pow(2.0, expo);
	

	//Define 3 dimensional array for all components in Event Array
	
	Mat3DDoub noiseTrace(allrecords.nGoodrec, 3, nPad); //-> effed up! because?
	Mat3DDoub postEventTrace(allrecords.nGoodrec, 3, nPad);
	Mat3DCmplx RFTrace(allrecords.nGoodrec, 2, nPad); 
	Mat3DDoub deltaRFTrace(allrecords.nGoodrec, 2, nPad); 
	
	// Set all the Five Headers useful for Statistics (& Selective Parsing?) Set file to Store.
	allrecords.setAzim(); allrecords.setEpic(); allrecords.setEvLong();
	allrecords.setEvLat(); allrecords.setEvMag();
	string evStatfn(outfn);
	evStatfn.append("_Event.txt");
	ofstream evStat(evStatfn.c_str()); 
	
	// Prep 3 files for Spectrum: Noise and Event spectrum
	string specNames[] = { "Vert_Spec.xyz", "Rad_Spec.xyz", "Tran_Spec.xyz", 
		"Noise_Spec.xyz"
	}; 
	ofstream specFNmes[4];
	for (int j = 0; j < 4 ; j++) {
		specFNmes[j].open(specNames[j].c_str());
	}
	//Prep 2 Files for Coherence Plot, See frequency behaviour
	string coherRadfn[] = { "Coher_Rad.xyz", "Coher_Tran.xyz" };
	ofstream coherFNmes[2];
	for (int j = 0; j < 2 ; j++) {
		coherFNmes[j].open(coherRadfn[j].c_str());
	}
	 
	
	if (getVerbose) cout << "Total Good Records: " << allrecords.nGoodrec << endl ;
	/*Create RF traces for only the Good records and print each individual trace?
	 Try adapting WriteRFdata to be used outside the Constructor... and
	 to write data to file outRadial.grid & outTrans.grid
	*/
	for (int irec =0; irec < allrecords.nGoodrec; irec++) {
		
		record = Sacread(allrecords.goodRecordname[irec], timeTag, lqtRotArg);
		tNoiseTrace = record.timeTag - record.timeStart - 3; // bias it by 3 sec
		nNoiseTrace = (tNoiseTrace) / (record.deltaT);
		
		/* 
		 **** URGENT UPDATE REQUIRED!!! This construction fails if record length following header is too small,
		 ****							Make check during record pass
		 Construct noiseTrace & postEventTrace
		 I do this for all the components: vertical, radial & transpose.
		 ? Do I need routine for clarity in code Prose?
		 */
		
		for (int cmp = 0; cmp < 3; cmp++) {
			
			for (int iter=0; iter < nPad; iter++) {
				if (iter < nNoiseTrace) {
					noiseTrace[irec][cmp][iter] = record.DataCmps[cmp][nNoiseTrace - iter];
				} else {
					noiseTrace[irec][cmp][iter] = 0.0;
				}
				
			}
			
			
			int posEventNxt, posEventStrt = nNoiseTrace;
			for (int iter=0; iter < nPad; iter++) {
				posEventNxt = posEventStrt + iter;
				if ( iter < nEventTrace) {
					postEventTrace[irec][cmp][iter] = record.DataCmps[cmp][posEventNxt];
				} else {
					postEventTrace[irec][cmp][iter] = 0.0;
				}
			}
			
		}
		
		
		allrecords.Azim[irec] = record.recBaz;
		allrecords.Epic[irec] = record.recEpic;
		allrecords.EvLong[irec] = record.evLon;
		allrecords.EvLat[irec] = record.evLat;
		allrecords.EvMag[irec] = record.evMag;
		
		// Event Statistics file is updated here, Five important headers.
		evStat << record.recBaz << "\t" << record.evLon <<  "\t" << record.evLat <<  "\t" << record.evMag
		<< "\t" << record.recEpic << endl;
		
		//cout << "File Name: " << allrecords.recordname[irec]
		//<< "Azimuth: " << allrecords.Azim[irec]
		//<< "Epicenter " << allrecords.Epic[irec] << endl;
		
		/*Pass Record to MTC driver, computes Spectral Estimates*/
		/*Args: Noise, Event, p, k, fmax */
		MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 2.5, 3, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
		// Code above doesn't work .. Check type
		
		if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << endl;
		MTCVals.addRF(RFTrace, irec);
		MTCVals.addCoher(deltaRFTrace, irec);
		
		// Dump Spectrum & Coherence For debug purposes
		for (int j = 0; j < 4 ; j++) {
			MTCVals.updateSpectrum(specFNmes[j], j, irec+1);
		}
		
		for (int j = 0; j < 2 ; j++) {
			MTCVals.updateCoher(coherFNmes[j], j, irec+1);
		}
		
	}
	
	
	/*Scan through records and stack them by azimuth and epicentre here
	 Check that getAzim & getEpic is set. TraceStack type is still unweaved.
	 Last two numbers represent time windows in the final stacked sections: pre & post
	 */
	string logfn(outfn);
	logfn.append(".log");
	ofstream logFile(logfn.c_str());   
	
	if (getAzim) {
		int flagFile = STACKAZIM;
		TraceStack azimStack(RFTrace, deltaRFTrace, allrecords.Azim, allrecords.Epic, Fcutoff, dT, 
							 10.0,  30.0, binsize);
		azimStack.stackAzim(bazmax, bazmin, bazinc, RFTrace, deltaRFTrace);
		azimStack.print(outfn, flagFile);
		logFile << azimStack.logReport;
	}
	
	
	
	if (getEpic) {
		int flagFile = STACKEPIC;
		/*Expand only when you have defined The function in TraceStack*/
		TraceStack epicStack(RFTrace, deltaRFTrace, allrecords.Azim, allrecords.Epic, 
						 Fcutoff, dT, 10.0, 30.0, binsize);
		epicStack.stackEpic(epbazmin, epbazmax, epicmin, epicmax, epicinc,RFTrace, deltaRFTrace);
		epicStack.print(outfn, flagFile);
		logFile << epicStack.logReport;
	}
	
    } 
	// end crude stack: else section! [ The general Record Analysis is same. Uses 'TraceStack' Instead of 'HarmonicStack']
	
	
	return 0;
	
}


