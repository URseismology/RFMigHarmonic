/**
 *	@file RFrouter.cpp
 *
 *
 *  @author Tolulope Olugboji
 *  @date 6/14/13.
 *  @copyright Yale University. All rights reserved.
 *
 *  @brief This class stages records, selects and determines the pre-RF methods to run.
 *
 *  The calls the required sacreader, computes RFs, stages RFs and determines stacking routine
 *	Depending on user requirements it then Routes the RF into the required stacking methods.
 *				if you want crude stack you call method.	runTraceStack.cpp
 *				if you want harmonic stack you call method.  runHarmonicStack.cpp
 *				if you want moving window migrated stack call. runMovingWindowStack.cpp
 *				if you want frequency dependend migrated stack.  migrateFrequencyStack.cpp
 *  ....   This cleans up the main module and enables reusability.
 *	It also decouples the sacread routine and allows time shifting if necessary...
 *
*/

#include "RFrouter.h"

// Define PhaseFlags for Converted Vs. Reverberated Phases.
#define PS 0
#define PPSMS 1
#define PPPMS 2

using namespace std;
		

RFrouter::RFrouter(char* eventsfn, char* outfn, float timeWin, float Fmax, float lqtArg, float horArg, int hTag, int stckArg, int migArg, char* velfn, int ROUTECODE,  bool getVbose)
:tEventTrace(timeWin), lqtRotArg(lqtArg), horRotArg(horArg), headerTag(hTag), Fcutoff(Fmax), nBoot(stckArg), targetDepthKm(migArg), velocityfn(velfn), getVerbose(getVbose){
	
	//cout << "You just called the RFrouter routine. WIP" << endl;
	
	// Pick First Record and Calculate timing details.
	// Initialize these values in class object for router...
	RecordList allRecords(eventsfn, timeWin, Fcutoff);
	
	std::string first;
	first.assign(allRecords.recordname[0]);
	Sacread record(first, headerTag, lqtRotArg);
	dT = record.deltaT;
	
	
	// Calculate Indexing of SAC Data Data
	// BUGGY Initialization of nNoiseTrace Captured.
	//   :: Initialize nNoiseTrace only for records with tagged headers. 
	nEventTrace = (tEventTrace) / (record.deltaT);
	tNoiseTrace = record.timeTag - record.timeStart - 3; // bias it by 3 sec
	nNoiseTrace = (tNoiseTrace) / (record.deltaT);
	nMax = MAX(nEventTrace, nNoiseTrace);
	
	
	/* Define pad length based on event number, to power of 2
	 Note here I pick maximum length of either the event or noise window
	 */
	
	int expo = int(log2(nMax-1)) + 1;
	nPad = pow(2.0, expo);					// Crucial Not to Re-initialize Here.
											// This is a class object.
	
	// Use Only Records that Pass Quality Check,
	bool flgChkMigrate = FALSE;
	if (ROUTECODE >= MWMAZIEPIC && ROUTECODE <= FREQMIGRATEHARMONIC) {
		flgChkMigrate = TRUE;
	}
	stageRecords(allRecords, flgChkMigrate);
	// Call Appropriate Stack Routine ...
	Router(ROUTECODE, outfn, allRecords);

	
	
};


// Constructor for Azimuth Stacking ...
RFrouter::RFrouter(char* eventsfn, char* outfn, float timeWin, float Fmax, float lqtArg, float horArg, int hTag, int stckArg, int migArg, char* velfn, int ROUTECODE, float bzmin, float bzmax, float bzinc, int bsze, bool getVbose)
:tEventTrace(timeWin), lqtRotArg(lqtArg), horRotArg(horArg), headerTag(hTag), Fcutoff(Fmax), nBoot(stckArg), targetDepthKm(migArg), velocityfn(velfn), isAzim(true), isEpic(false), bazmin(bzmin), bazmax(bzmax), bazinc(bzinc), binsize(bsze), getVerbose(getVbose){
	
	//cout << "You just called the Azimuthal Stacking Constructor. WIP" << endl;
	
	// Pick First Record and Calculate timing details.
	// Initialize these values in class object for router...
	RecordList allRecords(eventsfn, timeWin, Fcutoff);
	

	std::string first;
	first.assign(allRecords.recordname[0]);
	Sacread record(first, headerTag, lqtRotArg);
	dT = record.deltaT;
	
	
	// Calculate Indexing of SAC Data Data
	// @OLUGBOJI2016 - Reply to comment by kawakatsu, test deconvolution of noise to see if reverberations
	// are still present .... if timewindow negative then read noise.
	if (tEventTrace > 0) {
		tNoiseTrace = record.timeTag - record.timeStart - 3; // bias it by 3 sec
		nEventTrace = (tEventTrace) / (record.deltaT);
	}else{
		int startPoint = (tEventTrace + record.timeTag + 3 ); // wind backwards from start time
		tNoiseTrace = startPoint - record.timeStart; // no bias
		nEventTrace = startPoint / (record.deltaT);
	}
	
	nNoiseTrace = (tNoiseTrace) / (record.deltaT);
	nMax = MAX(nEventTrace, nNoiseTrace);
	
	
	/* Define pad length based on event number, to power of 2
	 Note here I pick maximum length of either the event or noise window
	 */
	
	int expo = int(log2(nMax-1)) + 1;
	nPad = pow(2.0, expo);				// Crucial Not to Re-initialize Here.
	
	// Use Only Records that Pass Quality Check,
	bool flgChkMigrate = FALSE;
	if (ROUTECODE >= MWMAZIEPIC && ROUTECODE <= FREQMIGRATEHARMONIC) {
		flgChkMigrate = TRUE;
	}
	stageRecords(allRecords, flgChkMigrate);
	// Call Appropriate Stack Routine
	Router(ROUTECODE, outfn, allRecords);
}


// Constructor for Epicentral Stacking ...
RFrouter::RFrouter(char* eventsfn, char* outfn, float timeWin, float Fmax, float lqtArg, float horArg, int hTag, int stckArg, int migArg, char* velfn, int ROUTECODE, float epbzmin, float epbzmax, float epmin, float epmax, float epinc, int bsze, bool getVbose)
:tEventTrace(timeWin), lqtRotArg(lqtArg), horRotArg(horArg), headerTag(hTag), Fcutoff(Fmax), nBoot(stckArg),  targetDepthKm(migArg), velocityfn(velfn), isAzim(false), isEpic(true), epbazmin(epbzmin), epbazmax(epbzmax), epicmin(epmin), epicmax(epmax), epicinc(epinc), binsize(bsze), getVerbose(getVbose){
	
	//cout << "You just called the Epicentral Stacking Constructor. WIP" << endl;
	
	// Pick First Record and Calculate timing details.
	// Initialize these values in class object for router...
	
	RecordList allRecords(eventsfn, timeWin, Fcutoff);
	std::string first;
	first.assign(allRecords.recordname[0]);
	Sacread record(first, headerTag, lqtRotArg);
	dT = record.deltaT;
	
	
	// Calculate Indexing of SAC Data Data
	tNoiseTrace = record.timeTag - record.timeStart - 3; // bias it by 3 sec
	nEventTrace = (tEventTrace) / (record.deltaT);
	nNoiseTrace = (tNoiseTrace) / (record.deltaT);
	nMax = MAX(nEventTrace, nNoiseTrace);
	
	
	/* Define pad length based on event number, to power of 2
	 Note here I pick maximum length of either the event or noise window
	 */
	
	int expo = int(log2(nMax-1)) + 1;
	nPad = pow(2.0, expo);				// Crucial Not to Re-initialize Here.
	
	
	// Use Only Records that Pass Quality Check,
	bool flgChkMigrate = FALSE;
	if (ROUTECODE >= MWMAZIEPIC && ROUTECODE <= FREQMIGRATEHARMONIC) {
		flgChkMigrate = TRUE;
	}
	stageRecords(allRecords, flgChkMigrate);
	// Call Appropriate Stack Routine.
	Router(ROUTECODE, outfn, allRecords);
}


void RFrouter::stageRecords(RecordList& parseRecords, bool flgChkMigrate){
	
	// Added an option to check that records can be migrated: flgChkMigrate
	MigrationParams velModel(velocityfn);
	

	if (getVerbose) cout << "Staging Records using tagged headers" <<
		"Total no of records: " << parseRecords.nTotrec << endl;
	
	/* Scan through Each SAC file and reconstitute record list based on
	 header entry. If entry not set, then skip record and put name in new record list
	 */
	for (int irec =0; irec < parseRecords.nTotrec; irec++) {
		Sacread record = Sacread(parseRecords.recordname[irec], headerTag, lqtRotArg);
		
		tNoiseTrace = record.timeTag - record.timeStart - 3;
		int checkNoiseTrace =  (tNoiseTrace) / (record.deltaT);
		
		/* Skip record with untagged field. If Tagged, Check to see that there's enough data for time window.
		 VERY IMPORTANT!!! USE record.timeTag and deltaT to do this.*/
		
		double timeShiftRight = 0.0;

		if (flgChkMigrate) {
			double rayParam = double(record.raySlowness);
			
			// Get time shift ... If ray is trapped then remover record ...
			// If looking for time delay for vertically impiging ray,
			// Use rayParam equals 1. since cos zero = 1.
			timeShiftRight = velModel.getTimeDelayPs(targetDepthKm, rayParam);
		}
		
		
		if (record.timeTag > 0 && timeShiftRight >= 0.0 ) {
			if( getVerbose) cout << "Record: " << irec+1 << " has the header no:  "<< headerTag << " Tagged " <<  "PhaseName: " << record.phaseName << endl;
			
			bool enoughDataEvent = (record.DataCmps.ncols() - ((record.timeTag - record.timeStart) / record.deltaT) ) >= nEventTrace;
			bool enoughDataNoise = ( (record.timeTag - record.timeStart) / record.deltaT ) >= checkNoiseTrace ;
			
			
			if (enoughDataNoise && enoughDataEvent ) {
				
				parseRecords.passNew(irec);
				// Pass Record only if there's enough data, Obviously!!
				// Less obvious... Only update nNoise if data is good.
				tNoiseTrace = record.timeTag - record.timeStart - 3; //bias it by 3 sec
				nNoiseTrace = (tNoiseTrace) / (record.deltaT);
			}
			if ( !enoughDataNoise || !enoughDataEvent) {
				if (getVerbose) cout << "Buggy data. Here is the cause of the Segmentation fault" << endl;
			}
			
		}else if (record.timeTag <= 0) {
			// Update here if
			if( getVerbose) cout << "Header no:  "<< headerTag << " is not tagged for record: " << irec+1 << endl;
		} else if (timeShiftRight < 0.0) {
			if( getVerbose) cout << "Record removed because earthquake is trapped at  "<< targetDepthKm << " for Record: " << irec+1 << endl;
		}
		
		
	}
	
	// Initialize All Header Variables.
	parseRecords.initHeaders();
	
	if (getVerbose) cout << "Total Good Records: " << parseRecords.nGoodrec << endl ;
		
	
	
}


void RFrouter::Router(int ROUTECODE, char* outfn, RecordList& allRecords){
	switch (ROUTECODE) {
		case SIMPLEAZIMEPIC:
			cout << "Routing... Azim-Epic Stack Without Migration.." << endl;
			callAzimEpicStack(outfn, allRecords);
			break;
		case SIMPLEHARMONIC:
			cout << "Routing ... Harmonic Stack Without Migration.." << endl;
			callHarmonicStack(outfn, allRecords);
			break;
		case MWMAZIEPIC:
		{
			cout << "Routing ... Azim-Epic Stack With Moving Window Migration" << endl;
			callAzimEpicMWMStack(outfn, allRecords);
			
			string fNameUpdate(outfn);
			fNameUpdate.append("timeDelay.txt");
			ofstream timedelayFile(fNameUpdate.c_str());
			timedelayFile << vertT << endl;
			
			break;
		}
		case MWMHARMONIC:
		{
			cout << "Routing ... Harmonic Stack With Moving Window Migration" << endl;
			callHarmonicMWMStack(outfn, allRecords);
			
			string fNameUpdate(outfn);
			fNameUpdate.append("timeDelay.txt");
			ofstream timedelayFile(fNameUpdate.c_str());
			timedelayFile << vertT << endl;
			
			break;
		}
		case FREQMIGRATEAZIMEPIC:
			break;
		case FREQMIGRATEHARMONIC:
			break;
		case ECCAZIMEPIC:
			cout << "Routing... Azim-Epic Stack Without Migration.." << endl;
			callAzimEpicECCStack(outfn, allRecords);
			break;
			break;
		case ECCHARMONIC:
			break;
		case MOHOREVERBPSMS:
		{
			cout << "Routing ... Azim-Epic Stack Migrating Reverberated Phase PpSms" << endl;
			int phaseFlag = PPSMS;
			callAzimEpicMOHOREVERBstack(outfn, allRecords, phaseFlag);
			break;
		}
		case MOHOREVERBPPMS:
		{
			cout << "Routing ... Azim-Epic Stack Migrating Reverberated Phase PpPms" << endl;
			int phaseFlag = PPPMS;
			callAzimEpicMOHOREVERBstack(outfn, allRecords, phaseFlag);
			break;
		}	
		case HARMONICREVERBPSMS:
		{
			cout << "Routing ... Harmonic Stack Migrating Reverberated Phase PpSms" << endl;
			int phaseFlag = PPSMS;
			callHarmonicReverbStack(outfn, allRecords, phaseFlag);
			
			string fNameUpdate(outfn);
			fNameUpdate.append("timeDelay.txt");
			ofstream timedelayFile(fNameUpdate.c_str());
			timedelayFile << vertT << endl;
			
			break;
		}
		case HARMONICREVERBPPMS:
		{
			cout << "Routing ... Harmonic Stack Migrating Reverberated Phase PpPms" << endl;
			int phaseFlag = PPPMS;
			callHarmonicReverbStack(outfn, allRecords, phaseFlag);
			
			string fNameUpdate(outfn);
			fNameUpdate.append("timeDelay.txt");
			ofstream timedelayFile(fNameUpdate.c_str());
			timedelayFile << vertT << endl;
			
			break;
		}
		case PRINTPHASETIME:
		{
			cout << "Routing ... Printing Phase Timing:" << endl;
			
			string fNameUpdate(outfn);
			fNameUpdate.append("_Epic_PhaseTime.txt");
			callPhasePrinter(const_cast<char*>(fNameUpdate.c_str()), allRecords);
			
			break;
		}
		default:
			break;
	}
	
}


// Actual Connections: Actual Modules called by the Router...

//1.
void RFrouter::callHarmonicStack(char* outfn, RecordList& parseRecords){
	
	
	Mat3DDoub noiseTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DDoub postEventTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DCmplx RFTrace(parseRecords.nGoodrec, 2, nPad);
	Mat3DDoub deltaRFTrace(parseRecords.nGoodrec, 2, nPad);
	
	
	// Build RFs Here .... Then Parse into Harmonic Stack...  
	for (int irec =0; irec < parseRecords.nGoodrec; irec++) {
		
		
		/*
		 **** URGENT UPDATE REQUIRED!!! Add single record rotation ...
		 */
		// ROTATION HACK FOR Ocean Bottom Data - April 23, 2014 . Included by OlugbojiTolulope  -
		// Japanese Data are for ocean bottom data [and Borehole Data] where the actual NE directions might be off
		
		/*
		 **** URGENT UPDATE REQUIRED!!! This construction fails if record length following header is too small,
		 ****							Make check during record pass
		 Construct noiseTrace & postEventTrace
		 I do this for all the components: vertical, radial & transpose.
		 ? Do I need routine for clarity in code Prose?
		 */
		
		// HACK - Nov. 10, 2013 . Included by OlugbojiTolulope  -
		// Japanese Data requires that predicted time changes a lot, so we need to shift
		// timing for first arrival for each record ....
		
		int doHorRotation = 0;   // Old constructor, no rotation.
		if (horRotArg > 0) doHorRotation = 1;
		
		switch (doHorRotation) {
			case 0:
			{
				Sacread record = Sacread(parseRecords.goodRecordname[irec], headerTag, lqtRotArg);
				
				tNoiseTrace = record.timeTag - record.timeStart - 3; //bias it by 3 sec
				nNoiseTrace = (tNoiseTrace) / (record.deltaT);
				
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
				
				
				parseRecords.Azim[irec] = record.recBaz;
				parseRecords.Epic[irec] = record.recEpic;
				parseRecords.EvLong[irec] = record.evLon;
				parseRecords.EvLat[irec] = record.evLat;
				parseRecords.EvMag[irec] = record.evMag;
				
				
				/*Pass Record to MTC driver, computes Spectral Estimates*/
				/*Args: Noise, Event, p, k, fmax */
				MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 2.5, 3, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
				// Code above doesn't work .. Check type
				
				if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << endl;
				MTCVals.addRF(RFTrace, irec);
				MTCVals.addCoher(deltaRFTrace, irec);
				
				break;
			}
				
			case 1:
			{
				Sacread record = Sacread(parseRecords.goodRecordname[irec], headerTag, lqtRotArg,
										 horRotArg, 1.0);
				
				tNoiseTrace = record.timeTag - record.timeStart - 3; //bias it by 3 sec
				nNoiseTrace = (tNoiseTrace) / (record.deltaT);
				
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
				
				
				parseRecords.Azim[irec] = record.recBaz;
				parseRecords.Epic[irec] = record.recEpic;
				parseRecords.EvLong[irec] = record.evLon;
				parseRecords.EvLat[irec] = record.evLat;
				parseRecords.EvMag[irec] = record.evMag;
				
				
				/*Pass Record to MTC driver, computes Spectral Estimates*/
				/*Args: Noise, Event, p, k, fmax */
				MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 2.5, 3, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
				// Code above doesn't work .. Check type
				
				if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << endl;
				MTCVals.addRF(RFTrace, irec);
				MTCVals.addCoher(deltaRFTrace, irec);
				
				break;
			}
	
		}
		
		
	}
	
	double vertTshft = 0.0;
	
	// Dump RFs into Harmonic Stack ...
	{
		
		int flagFile = CMPSHARMONIC;
		
		int nRecs = parseRecords.nGoodrec, ncmps = 2, nFreq = nPad;
		HarmonicStack testHStack(nPad, ncmps, nRecs, parseRecords.Azim, parseRecords.Epic, Fcutoff, dT, 10.0,  100.0, nBoot);
		testHStack.regressNBootTimes(RFTrace, deltaRFTrace);
		testHStack.print(outfn, flagFile, vertTshft);
		//logFile << testHStack.logReport;
	}
	

	
}

//2.
void RFrouter::callAzimEpicStack(char* outfn, RecordList& parseRecords){
	
	
	Mat3DDoub noiseTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DDoub postEventTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DCmplx RFTrace(parseRecords.nGoodrec, 2, nPad);
	Mat3DDoub deltaRFTrace(parseRecords.nGoodrec, 2, nPad);

	// Build RFs Here .... Then Parse into TraceStack...
	for (int irec =0; irec < parseRecords.nGoodrec; irec++) {
		
		
		/*
		 **** URGENT UPDATE REQUIRED!!! Add single record rotation ...
		 */
		// ROTATION HACK FOR Ocean Bottom Data - April 23, 2014 . Included by OlugbojiTolulope  -
		// Japanese Data are for ocean bottom data [and Borehole Data] where the actual NE directions might be off
		
		
		/*
		 **** URGENT UPDATE REQUIRED!!! This construction fails if record length following header is too small,
		 ****							Make check during record pass
		 Construct noiseTrace & postEventTrace
		 I do this for all the components: vertical, radial & transpose.
		 ? Do I need routine for clarity in code Prose?
		 */
		
		// HACK - Nov. 10, 2013 . Included by OlugbojiTolulope  -
		// Japanese Data requires that predicted time changes a lot, so we need to shift
		// timing for first arrival for each record ....
		
		int doHorRotation = 0;
		if (horRotArg > 0) doHorRotation = 1;
		
		switch (doHorRotation) {
			case 0:
			{
				Sacread record = Sacread(parseRecords.goodRecordname[irec], headerTag, lqtRotArg);
				
				tNoiseTrace = record.timeTag - record.timeStart - 3; //bias it by 3 sec
				nNoiseTrace = (tNoiseTrace) / (record.deltaT);
				
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
				
				
				parseRecords.Azim[irec] = record.recBaz;
				parseRecords.Epic[irec] = record.recEpic;
				parseRecords.EvLong[irec] = record.evLon;
				parseRecords.EvLat[irec] = record.evLat;
				parseRecords.EvMag[irec] = record.evMag;
				
				
				/*Pass Record to MTC driver, computes Spectral Estimates*/
				/*Args: Noise, Event, p, k, fmax */
				MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 2.5, 3, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
				// Code above doesn't work .. Check type
				
				if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << "timeTag: " << tNoiseTrace << "nTag: " << nNoiseTrace << endl;
				MTCVals.addRF(RFTrace, irec);
				MTCVals.addCoher(deltaRFTrace, irec);
				
				break;
			}
				
			case 1:
			{
				Sacread record = Sacread(parseRecords.goodRecordname[irec], headerTag, lqtRotArg,
										 horRotArg, 1.0);
				
				tNoiseTrace = record.timeTag - record.timeStart - 3; //bias it by 3 sec
				nNoiseTrace = (tNoiseTrace) / (record.deltaT);
				
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
				
				
				parseRecords.Azim[irec] = record.recBaz;
				parseRecords.Epic[irec] = record.recEpic;
				parseRecords.EvLong[irec] = record.evLon;
				parseRecords.EvLat[irec] = record.evLat;
				parseRecords.EvMag[irec] = record.evMag;
				
				
				/*Pass Record to MTC driver, computes Spectral Estimates*/
				/*Args: Noise, Event, p, k, fmax */
				MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 2.5, 3, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
				// Code above doesn't work .. Check type
				
				if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << "timeTag: " << tNoiseTrace << "nTag: " << nNoiseTrace << endl;
				MTCVals.addRF(RFTrace, irec);
				MTCVals.addCoher(deltaRFTrace, irec);
				
				break;
			}
		}
				
	
		
		
	}
		
	if (getVerbose) cout << "Done Building RFs " << endl ;
	
	// Dump RFs into TraceStack ...
	
	
		
		if (isAzim) {
			
			int flagFile = STACKAZIM;
			
			TraceStack azimStack(RFTrace, deltaRFTrace, parseRecords.Azim, parseRecords.Epic, Fcutoff, dT, 10.0,  100.0, binsize);
			azimStack.stackAzim(bazmax, bazmin, bazinc, RFTrace, deltaRFTrace);
			azimStack.print(outfn, flagFile, 0.0);
			
			//if (getVerbose) cout << "Done Azim Stacking " << endl ;
			//if (getVerbose) cout << azimStack.logReport << endl ;
			//logFile << azimStack.logReport;
		}
		
		
		
		if (isEpic) {
			
			//int flagFile = STACKEPIC;  // test Jacknife here, then  code entries to toggle jacknife option ...
			int flagFile = STACKEPICJCK;
			
			
			TraceStack epicStack(RFTrace, deltaRFTrace, parseRecords.Azim, parseRecords.Epic,
								 Fcutoff, dT, 10.0, 100.0, binsize);
			//epicStack.stackEpic(epbazmin, epbazmax, epicmin, epicmax, epicinc,RFTrace, deltaRFTrace);
			
			epicStack.jacknifeStackEpic(epbazmin, epbazmax, epicmin, epicmax, epicinc,RFTrace, deltaRFTrace);
			
			
			epicStack.print(outfn, flagFile, 0.0);
			
			//if (getVerbose) cout << "Done Epic Stacking " << endl ;
			//if (getVerbose) cout << epicStack.logReport << endl ;
			//logFile << epicStack.logReport;
		}
		
		//if (getVerbose) cout << "Done Stacking " << endl ;
		//logFile << testHStack.logReport;
	
	
	
	
}

//2b.
void RFrouter::callAzimEpicECCStack(char* outfn, RecordList& parseRecords){
	
	
	Mat3DDoub noiseTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DDoub postEventTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DCmplx RFTrace(parseRecords.nGoodrec, 2, nPad);
	Mat3DDoub deltaRFTrace(parseRecords.nGoodrec, 2, nPad);

	// Build RFs Here .... Then Parse into TraceStack...
	for (int irec =0; irec < parseRecords.nGoodrec; irec++) {
		
		
		/*
		 **** URGENT UPDATE REQUIRED!!! Add single record rotation ...
		 */
		// ROTATION HACK FOR Ocean Bottom Data - April 23, 2014 . Included by OlugbojiTolulope  -
		// Japanese Data are for ocean bottom data [and Borehole Data] where the actual NE directions might be off
		
		
		/*
		 **** URGENT UPDATE REQUIRED!!! This construction fails if record length following header is too small,
		 ****							Make check during record pass
		 Construct noiseTrace & postEventTrace
		 I do this for all the components: vertical, radial & transpose.
		 ? Do I need routine for clarity in code Prose?
		 */
		
		// HACK - Nov. 10, 2013 . Included by OlugbojiTolulope  -
		// Japanese Data requires that predicted time changes a lot, so we need to shift
		// timing for first arrival for each record ....
		
		int doHorRotation = 0;
		if (horRotArg > 0) doHorRotation = 1;
		
		switch (doHorRotation) {
			case 0:
			{
				Sacread record = Sacread(parseRecords.goodRecordname[irec], headerTag, lqtRotArg);
				
				tNoiseTrace = record.timeTag - record.timeStart - 3; //bias it by 3 sec
				nNoiseTrace = (tNoiseTrace) / (record.deltaT);
				
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
				
				
				parseRecords.Azim[irec] = record.recBaz;
				parseRecords.Epic[irec] = record.recEpic;
				parseRecords.EvLong[irec] = record.evLon;
				parseRecords.EvLat[irec] = record.evLat;
				parseRecords.EvMag[irec] = record.evMag;
				
				
				/*Pass Record to MTC driver, computes Spectral Estimates*/
				/*Args: Noise, Event, p=1, k=1, fmax */
				MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 1, 1, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
				// Code above doesn't work .. Check type
				
				if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << "timeTag: " << tNoiseTrace << "nTag: " << nNoiseTrace << endl;
				MTCVals.addRF(RFTrace, irec);
				MTCVals.addCoher(deltaRFTrace, irec);
				
				break;
			}
				
			case 1:
			{
				Sacread record = Sacread(parseRecords.goodRecordname[irec], headerTag, lqtRotArg,
										 horRotArg, 1.0);
				
				tNoiseTrace = record.timeTag - record.timeStart - 3; //bias it by 3 sec
				nNoiseTrace = (tNoiseTrace) / (record.deltaT);
				
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
				
				
				parseRecords.Azim[irec] = record.recBaz;
				parseRecords.Epic[irec] = record.recEpic;
				parseRecords.EvLong[irec] = record.evLon;
				parseRecords.EvLat[irec] = record.evLat;
				parseRecords.EvMag[irec] = record.evMag;
				
				
				/*Pass Record to MTC driver, computes Spectral Estimates*/
				/*Args: Noise, Event, p=1, k=1, fmax */
				MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 1, 1, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
				// Code above doesn't work .. Check type
				
				if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << "timeTag: " << tNoiseTrace << "nTag: " << nNoiseTrace << endl;
				MTCVals.addRF(RFTrace, irec);
				MTCVals.addCoher(deltaRFTrace, irec);
				
				break;
			}
		}
				
	
		
		
	}
		
	if (getVerbose) cout << "Done Building RFs " << endl ;
	
	// Dump RFs into TraceStack ...
	
	
		
		if (isAzim) {
			
			int flagFile = STACKAZIM;
			
			TraceStack azimStack(RFTrace, deltaRFTrace, parseRecords.Azim, parseRecords.Epic, Fcutoff, dT, 10.0,  100.0, binsize);
			azimStack.stackAzim(bazmax, bazmin, bazinc, RFTrace, deltaRFTrace);
			azimStack.print(outfn, flagFile, 0.0);
			
			//if (getVerbose) cout << "Done Azim Stacking " << endl ;
			//if (getVerbose) cout << azimStack.logReport << endl ;
			//logFile << azimStack.logReport;
		}
		
		
		
		if (isEpic) {
			
			//int flagFile = STACKEPIC;  // test Jacknife here, then  code entries to toggle jacknife option ...
			int flagFile = STACKEPICJCK;
			
			
			TraceStack epicStack(RFTrace, deltaRFTrace, parseRecords.Azim, parseRecords.Epic,
								 Fcutoff, dT, 10.0, 100.0, binsize);
			//epicStack.stackEpic(epbazmin, epbazmax, epicmin, epicmax, epicinc,RFTrace, deltaRFTrace);
			
			epicStack.jacknifeStackEpic(epbazmin, epbazmax, epicmin, epicmax, epicinc,RFTrace, deltaRFTrace);
			
			
			epicStack.print(outfn, flagFile, 0.0);
			
			//if (getVerbose) cout << "Done Epic Stacking " << endl ;
			//if (getVerbose) cout << epicStack.logReport << endl ;
			//logFile << epicStack.logReport;
		}
		
		//if (getVerbose) cout << "Done Stacking " << endl ;
		//logFile << testHStack.logReport;
	
	
	
	
}

//3.
void RFrouter::callAzimEpicMWMStack(char* outfn, RecordList& parseRecords){
	
	string eventStats; // OLUGBOJI 2016: Save event statistics for reproducible reporting.
	
	Mat3DDoub noiseTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DDoub postEventTrace(parseRecords.nGoodrec, 3, nPad);
	
	Mat3DCmplx RFTrace(parseRecords.nGoodrec, 2, nPad);
	Mat3DDoub deltaRFTrace(parseRecords.nGoodrec, 2, nPad);
	
	// Trace for save and display...
	Mat3DDoub allTracesRot(parseRecords.nGoodrec, 3, 2*nPad);
	Mat3DDoub allTracesRw(parseRecords.nGoodrec, 3, 2*nPad);
	// Save VH coherence and all spectrum for display
	Mat3DDoub coherVH(parseRecords.nGoodrec, 2, nPad);
	Mat3DDoub specNseEvnt(parseRecords.nGoodrec, 6, nPad);
	
	
	
	
	// Use target depth to find layer within the modelStack where migration should begin
	// This layer will then become the bottomStack for migration.
	MigrationParams velModel(velocityfn);
	double timeShiftRight = 0.0;
	double minTshft = 0.0;
	double maxTshft = 0.0;
	double vertTshft = velModel.getVertTimeDelay(targetDepthKm, PS);
	
	// timeshiftRight: The event window is moved(moving window migration - MWM)
	// depending on:
	// 1. The time scaling implied by the rayParameter
	// 2. The timedelay for each velocity layer and
	
	eventStats.append("# rec no \t timeTag \t phaseName \t filePath \t Date ...\n");
	// Compute timeShiftRight here...
	for (int irec =0; irec < parseRecords.nGoodrec; irec++) {
		
		
		Sacread record =
					Sacread(parseRecords.goodRecordname[irec], headerTag, lqtRotArg);
		

		
		// HACK - Nov. 10, 2013 . Included by OlugbojiTolulope  -
		// Japanese Data requires that predicted time changes a lot, so we need to shift
		// timing for first arrival for each record ....
		tNoiseTrace = record.timeTag - record.timeStart - 3; //bias it by 3 sec
		
		// @OLUGBOJI2016 - Reply to comment by kawakatsu, test deconvolution of noise to see if reverberations
		// are still present .... if timewindow negative then read noise Only.
		if (tEventTrace > 0) {
			tNoiseTrace = record.timeTag - record.timeStart - 3; // bias it by 3 sec
			nEventTrace = (tEventTrace) / (record.deltaT);
		}else{
			int startPoint = (tEventTrace + record.timeTag - 10 ); // wind backwards from start time with 10 seconds error tor start time
			tNoiseTrace = startPoint - record.timeStart ; // no bias
			nEventTrace = (startPoint -2) / (record.deltaT);
		}
		
		nNoiseTrace = (tNoiseTrace) / (record.deltaT);
		
		
		// OLUGBOJI2016 see if this works for reporting ...
		eventStats.append( to_string(irec+1) );
		eventStats.append("  ");
		eventStats.append( to_string(tNoiseTrace) );
		eventStats.append("  ");
		eventStats.append( to_string(record.phaseName) );
		eventStats.append("  ");
		// debug timing and windowing variables.
		//eventStats.append("nEv,nNse: ");
		//eventStats.append( to_string(nEventTrace) );
		//eventStats.append("  ");
		//eventStats.append( to_string(nNoiseTrace) );
		eventStats.append(record.eventStats);
		
		double rayParam = double(record.raySlowness);				// in sec/Km
		
		
		// Get time shift for each record.. store max and min. 
		// Then check with vertical time delay.
		// use this time delays to bias time record during display
		timeShiftRight = velModel.getTimeDelayPs(targetDepthKm, rayParam);
		cout << "Time Delay: " << timeShiftRight << "s" << endl;
		
		
		if (irec == 0 ) {
			minTshft = timeShiftRight;
			maxTshft = timeShiftRight;
		} else {
			if (timeShiftRight < minTshft) {
				minTshft = timeShiftRight;
			} else if (timeShiftRight > maxTshft) {
				maxTshft = timeShiftRight;
			}
		}
		
		
		// Noise Is Noise. Just Update Noise ...
		for (int cmp = 0; cmp < 3; cmp++) {
			
			for (int iter=0; iter < nPad; iter++) {
				if (iter < nNoiseTrace) {
					noiseTrace[irec][cmp][iter] = record.DataCmps[cmp][nNoiseTrace - iter];
					
					// just put full traces [Noise] for display, no need to pad.
					allTracesRot[irec][cmp][iter] = record.DataCmps[cmp][nNoiseTrace - iter];
					allTracesRw[irec][cmp][iter] = record.DataCmpsRw[cmp][nNoiseTrace - iter];
					// end data print construct ...
					
					
				} else {
					noiseTrace[irec][cmp][iter] = 0.0;
				}
				
			}
			
			// Data is Migrated Using timeShift ... Only for R, T. Skip shift for Z.
			int nShift = (timeShiftRight) / (record.deltaT);
			if (cmp == 0) {
				nShift = 0;
			}
			int posEventNxt, posEventStrt = nNoiseTrace + nShift;
			
			cout << "i: " << irec << "cmp: " << cmp << "nse: " << nNoiseTrace  << "end: " << posEventStrt+nEventTrace << "nPad: " << nPad << endl;
			
			for (int iter=0; iter < nPad; iter++) {
				posEventNxt = posEventStrt + iter;
				if ( iter < nEventTrace) {
					postEventTrace[irec][cmp][iter] = record.DataCmps[cmp][posEventNxt];
					
					// just put full traces [event] for display, no need to pad.
					allTracesRot[irec][cmp][iter+nNoiseTrace] = record.DataCmps[cmp][posEventNxt];
					allTracesRw[irec][cmp][iter+nNoiseTrace] = record.DataCmpsRw[cmp][posEventNxt];
					// end data print construct ...
					
				} else {
					postEventTrace[irec][cmp][iter] = 0.0;
				}
			}
			
		}
		
		parseRecords.setHeaders(record.recBaz, record.recEpic, record.evLon,
								record.evLat, record.evMag, irec);
		
		/*Pass Record to MTC driver, computes Spectral Estimates*/
		/*Args: Noise, Event, p, k, fmax */
		MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 2.5, 3, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
		
		if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << endl;
		MTCVals.addRF(RFTrace, irec);
		MTCVals.addCoher(deltaRFTrace, irec);
		
		MTCVals.addCorrVH(coherVH, irec);
		MTCVals.addSpecAll(specNseEvnt, irec);
		
		
		
	}
	
	minT = minTshft; maxT = maxTshft; vertT = vertTshft;
	
	if (getVerbose) cout << "No of Layers + Halfspace is = " << velModel.getNoLayers() + 1 << endl;
	
	if (getVerbose) cout << "Done Building RFs " << endl ;
	
	// Dump RFs into TraceStack ...
	
	string fNameUpdate(outfn);
	fNameUpdate.append("Migrate");
	
	outfn = const_cast<char*>( fNameUpdate.c_str() );
	
	if (isAzim) {
		
		int flagFile = STACKAZIM;
		
		TraceStack azimStack(RFTrace, deltaRFTrace, parseRecords.Azim, parseRecords.Epic, Fcutoff, dT, 10.0,  100.0, binsize);
		azimStack.stackAzim(bazmax, bazmin, bazinc, RFTrace, deltaRFTrace);
		azimStack.print(outfn, flagFile, 0.0);
		
		// Dump Vertical Time delay in File for use by plotting utility..
		
		
		//if (getVerbose) cout << "Done Azim Stacking " << endl ;
		//if (getVerbose) cout << azimStack.logReport << endl ;
		//logFile << azimStack.logReport;
	}
	
	
	
	if (isEpic) {
		
		//int flagFile = STACKEPIC;
		int flagFile = STACKEPICJCK;
		
		TraceStack epicStack(RFTrace, deltaRFTrace, parseRecords.Azim, parseRecords.Epic,
							 Fcutoff, dT, 10.0, 100.0, binsize);
		//epicStack.stackEpic(epbazmin, epbazmax, epicmin, epicmax, epicinc,RFTrace, deltaRFTrace);
		
		epicStack.jacknifeStackEpic(epbazmin, epbazmax, epicmin, epicmax, epicinc,RFTrace, deltaRFTrace);
		
		
		MatDoub parseTime, parseDepth;
		parseTime = epicStack.getTimeAscii(); // return time from trace stack.
		
		cout << "Testing Migrating Time to Depth. All Times, Size: " <<
		parseTime.nrows() << "   " << parseTime.ncols() <<
		"First and Last: " << parseTime[0][0] << "  " << parseTime[0][parseTime.ncols()-3] << endl;
		
		// Write internal function to migrate time to depth ...
		parseDepth = migrateTime2Depth(parseTime, vertT);
		//parseDepth = migrateTime2Depth(parseTime, 0.0);
		epicStack.setDepthAscii(parseDepth);
		epicStack.print(outfn, flagFile, 0.0);
		
		//if (getVerbose) cout << "Done Epic Stacking " << endl ;
		//if (getVerbose) cout << epicStack.logReport << endl ;
		//logFile << epicStack.logReport;
	}
	
	//if (getVerbose) cout << "Done Stacking " << endl ;
	//logFile << testHStack.logReport;
	
	
	// ************* When all is done, print the statistics ...
	cout << "Saving eqstats, waveforms and coherence" << endl;
	printEQStats(outfn, eventStats);
	printDataTraces(outfn, parseRecords.nGoodrec, allTracesRw, allTracesRot);
	printCoherSpec(outfn, parseRecords.nGoodrec,  coherVH, specNseEvnt);
	cout << "end save" << endl;


}

//4.
void RFrouter::callAzimEpicMOHOREVERBstack(char* outfn, RecordList& parseRecords, int phaseFlag){
	
	Mat3DDoub noiseTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DDoub postEventTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DCmplx RFTrace(parseRecords.nGoodrec, 2, nPad);
	Mat3DDoub deltaRFTrace(parseRecords.nGoodrec, 2, nPad);
	
	
	// Use target depth to find layer within the modelStack where migration should begin
	// This layer will then become the bottomStack for migration.
	MigrationParams velModel(velocityfn);
	double timeShiftRight = 0.0;
	double minTshft = 0.0;
	double maxTshft = 0.0;
	double vertTshft = 0.0;
	
	switch (phaseFlag) {
		case PPSMS:
			vertTshft = velModel.getVertTimeDelay(targetDepthKm, PPSMS);
			break;
		case PPPMS:
			vertTshft = velModel.getVertTimeDelay(targetDepthKm, PPPMS);
			break;
		default:
			break;
	}
	
	// timeshiftRight: The event window is moved(moving window migration - MWM)
	// depending on:
	// 1. The time scaling implied by the rayParameter
	// 2. The timedelay for each velocity layer and
	
	
	// Compute timeShiftRight here...
	for (int irec =0; irec < parseRecords.nGoodrec; irec++) {
		
		
		Sacread record =
		Sacread(parseRecords.goodRecordname[irec], headerTag, lqtRotArg);
		
		// HACK - Nov. 10, 2013 . Included by OlugbojiTolulope  -
		// Japanese Data requires that predicted time changes a lot, so we need to shift
		// timing for first arrival for each record ....
		tNoiseTrace = record.timeTag - record.timeStart - 3; //bias it by 3 sec
		nNoiseTrace = (tNoiseTrace) / (record.deltaT);
		
		double rayParam = double(record.raySlowness);				// in sec/Km
		
		
		// Get time shift ... If ray is trapped then remover record ...
		// If looking for time delay for vertically impiging ray,
		// Use rayParam equals 1. since cos zero = 1.
		switch (phaseFlag) {
			case PPSMS:
				timeShiftRight = velModel.getTimeDelayPpSms(targetDepthKm, rayParam);
				cout << "Time Delay PpSms: " << timeShiftRight << "s" << endl;
				break;
			case PPPMS:
				timeShiftRight = velModel.getTimeDelayPpPms(targetDepthKm, rayParam);
				cout << "Time Delay PpPms: " << timeShiftRight << "s" << endl;
				break;
			default:
				break;
		}
		
		// Update minimum and maximum time shifts.
		if (irec == 0 ) {
			minTshft = timeShiftRight;
			maxTshft = timeShiftRight;
		} else {
			if (timeShiftRight < minTshft) {
				minTshft = timeShiftRight;
			} else if (timeShiftRight > maxTshft) {
				maxTshft = timeShiftRight;
			}
		}
		
		// Noise Is Noise. Just Update Noise ...
		for (int cmp = 0; cmp < 3; cmp++) {
			
			for (int iter=0; iter < nPad; iter++) {
				if (iter < nNoiseTrace) {
					noiseTrace[irec][cmp][iter] = record.DataCmps[cmp][nNoiseTrace - iter];
				} else {
					noiseTrace[irec][cmp][iter] = 0.0;
				}
				
			}
			
			// Data is Migrated Using timeShift ... Only for R, T. Skip shift for Z.
			int nShift = (timeShiftRight) / (record.deltaT);
			if (cmp == 0) {
				nShift = 0;
			}
			int posEventNxt, posEventStrt = nNoiseTrace + nShift;
			
			for (int iter=0; iter < nPad; iter++) {
				posEventNxt = posEventStrt + iter;
				if ( iter < nEventTrace) {
					postEventTrace[irec][cmp][iter] = record.DataCmps[cmp][posEventNxt];
				} else {
					postEventTrace[irec][cmp][iter] = 0.0;
				}
			}
			
		}
		
		parseRecords.setHeaders(record.recBaz, record.recEpic, record.evLon,
								record.evLat, record.evMag, irec);
		
		/*Pass Record to MTC driver, computes Spectral Estimates*/
		/*Args: Noise, Event, p, k, fmax */
		MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 2.5, 3, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
		
		if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << endl;
		MTCVals.addRF(RFTrace, irec);
		MTCVals.addCoher(deltaRFTrace, irec);
		
		
		
		
		
		
	}
	
	minT = minTshft; maxT = maxTshft; vertT = vertTshft;
	
	if (getVerbose) cout << "No of Layers + Halfspace is = " << velModel.getNoLayers() + 1 << endl;
	
	if (getVerbose) cout << "Done Building RFs " << endl ;
	
	// Dump RFs into TraceStack ...
	
	string fNameUpdate(outfn);
	fNameUpdate.append("MigrateReverb");
	
	outfn = const_cast<char*>( fNameUpdate.c_str() );
	
	if (isAzim) {
		
		int flagFile = STACKAZIM;
		
		TraceStack azimStack(RFTrace, deltaRFTrace, parseRecords.Azim, parseRecords.Epic, Fcutoff, dT, 10.0,  100.0, binsize);
		azimStack.stackAzim(bazmax, bazmin, bazinc, RFTrace, deltaRFTrace);
		azimStack.print(outfn, flagFile, 0.0);
		
		//if (getVerbose) cout << "Done Azim Stacking " << endl ;
		//if (getVerbose) cout << azimStack.logReport << endl ;
		//logFile << azimStack.logReport;
	}
	
	
	
	if (isEpic) {
		
		int flagFile = STACKEPIC;
		
		
		TraceStack epicStack(RFTrace, deltaRFTrace, parseRecords.Azim, parseRecords.Epic,
							 Fcutoff, dT, 10.0, 100.0, binsize);
		epicStack.stackEpic(epbazmin, epbazmax, epicmin, epicmax, epicinc,RFTrace, deltaRFTrace);
		epicStack.print(outfn, flagFile, 0.0);

		
		//if (getVerbose) cout << "Done Epic Stacking " << endl ;
		//if (getVerbose) cout << epicStack.logReport << endl ;
		//logFile << epicStack.logReport;
	}
	
	//if (getVerbose) cout << "Done Stacking " << endl ;
	//logFile << testHStack.logReport;
	
	
}

//5. Migrated Harmonic Stacks...  Add Depth Migration [Also Amplitude and Azimuth -]
void RFrouter::callHarmonicMWMStack(char* outfn, RecordList& parseRecords){
	
	Mat3DDoub noiseTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DDoub postEventTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DCmplx RFTrace(parseRecords.nGoodrec, 2, nPad);
	Mat3DDoub deltaRFTrace(parseRecords.nGoodrec, 2, nPad);
	
	
	// Use target depth to find layer within the modelStack where migration should begin
	// This layer will then become the bottomStack for migration.
	MigrationParams velModel(velocityfn);
	double timeShiftRight = 0.0;
	double minTshft = 0.0;
	double maxTshft = 0.0;
	double vertTshft = velModel.getVertTimeDelay(targetDepthKm, PS);
	
	// timeshiftRight: The event window is moved(moving window migration - MWM)
	// depending on:
	// 1. The time scaling implied by the rayParameter
	// 2. The timedelay for each velocity layer and
	
	
	// Compute timeShiftRight here...
	for (int irec =0; irec < parseRecords.nGoodrec; irec++) {
		
		
		Sacread record =
		Sacread(parseRecords.goodRecordname[irec], headerTag, lqtRotArg);
		
		// HACK - Nov. 10, 2013 . Included by OlugbojiTolulope  -
		// Japanese Data requires that predicted time changes a lot, so we need to shift
		// timing for first arrival for each record ....
		tNoiseTrace = record.timeTag - record.timeStart - 3; //bias it by 3 sec - NOT????
		
		// @OLUGBOJI2016 - Reply to comment by kawakatsu, test deconvolution of noise to see if reverberations
		// are still present .... if timewindow negative then read noise Only.
		if (tEventTrace > 0) {
			tNoiseTrace = record.timeTag - record.timeStart - 3; // bias it by 3 sec
			nEventTrace = (tEventTrace) / (record.deltaT);
		}else{
			int startPoint = (tEventTrace + record.timeTag - 10 ); // wind backwards from start time with 10 seconds error tor start time
			tNoiseTrace = startPoint - record.timeStart ; // no bias
			nEventTrace = startPoint / (record.deltaT);
		}
		
		nNoiseTrace = (tNoiseTrace) / (record.deltaT);
		
		double rayParam = double(record.raySlowness);				// in sec/Km
		
		
		// Get time shift ... If ray is trapped then remover record ...
		// If looking for time delay for vertically impiging ray,
		// Use rayParam equals 1. since cos zero = 1.
		timeShiftRight = velModel.getTimeDelayPs(targetDepthKm, rayParam);
		cout << "Time Delay: " << timeShiftRight << "s"  << "tNoise: "<< tNoiseTrace << endl;
		
		if (irec == 0 ) {
			minTshft = timeShiftRight;
			maxTshft = timeShiftRight;
		} else {
			if (timeShiftRight < minTshft) {
				minTshft = timeShiftRight;
			} else if (timeShiftRight > maxTshft) {
				maxTshft = timeShiftRight;
			}
		}
		
		// Noise Is Noise. Just Update Noise ...
		for (int cmp = 0; cmp < 3; cmp++) {
			
			for (int iter=0; iter < nPad; iter++) {
				if (iter < nNoiseTrace) {
					noiseTrace[irec][cmp][iter] = record.DataCmps[cmp][nNoiseTrace - iter];
				} else {
					noiseTrace[irec][cmp][iter] = 0.0;
				}
				
			}
			
			// Data is Migrated Using timeShift ... Only for R, T. Skip shift for Z.
			int nShift = (timeShiftRight) / (record.deltaT);
			cout << "start: " << record.timeStart << endl;
			if (cmp == 0) {
				nShift = 0;
			}
			int posEventNxt, posEventStrt = nNoiseTrace + nShift;
			
			for (int iter=0; iter < nPad; iter++) {
				posEventNxt = posEventStrt + iter;
				if ( iter < nEventTrace) {
					postEventTrace[irec][cmp][iter] = record.DataCmps[cmp][posEventNxt];
				} else {
					postEventTrace[irec][cmp][iter] = 0.0;
				}
			}
			
		}
		
		parseRecords.setHeaders(record.recBaz, record.recEpic, record.evLon,
								record.evLat, record.evMag, irec);
		
		/*Pass Record to MTC driver, computes Spectral Estimates*/
		/*Args: Noise, Event, p, k, fmax */
		MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 2.5, 3, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
		
		if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << endl;
		MTCVals.addRF(RFTrace, irec);
		MTCVals.addCoher(deltaRFTrace, irec);
		
		
		
		
		
		
	}
	
	minT = minTshft; maxT = maxTshft; vertT = vertTshft;
	
	if (getVerbose) cout << "No of Layers + Halfspace is = " << velModel.getNoLayers() + 1 << endl;
	
	if (getVerbose) cout << "Done Building RFs " << endl ;
	
	// Dump RFs into TraceStack ...
	
	string fNameUpdate(outfn);
	fNameUpdate.append("MigrateHarmonic");
	
	outfn = const_cast<char*>( fNameUpdate.c_str() );
	
	// Dump RFs into Harmonic Stack ...
	{
		
		int flagFile = CMPSHARMONIC;
		
		int nRecs = parseRecords.nGoodrec, ncmps = 2, nFreq = nPad;
		HarmonicStack testHStack(nPad, ncmps, nRecs, parseRecords.Azim, parseRecords.Epic, Fcutoff, dT, 10.0,  100.0, nBoot);
		testHStack.regressNBootTimes(RFTrace, deltaRFTrace);
		
		
		//  *************************************** Compute Depth Migration Here ...
		MatDoub parseTime, parseDepth;
		parseTime = testHStack.getTimeAscii(); // return time from trace stack.
		
		
		// Write internal function to migrate time to depth ...
		parseDepth = migrateTime2Depth(parseTime, vertT);
		testHStack.setDepthAscii(parseDepth);
		
		// **************************************** End Depth Migration Here ....
		
		testHStack.print(outfn, flagFile, 0.0);
		//logFile << testHStack.logReport;
		
		// after regression plot for debug purposes the spectrum of the constant term
		// @OLUGBOJI 2016 -- remove this printer if not needed ...
		int indxConstTerm = 0;
		testHStack.saveSpecMeanBAZcmps(indxConstTerm, outfn);
		// @END OLUGBOJI2016 ...

	}

	
}


//6. Harmonic Stacks Migrated With Reverberated Phases.....
void RFrouter::callHarmonicReverbStack(char* outfn, RecordList& parseRecords, int phaseFlag){
	
	Mat3DDoub noiseTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DDoub postEventTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DCmplx RFTrace(parseRecords.nGoodrec, 2, nPad);
	Mat3DDoub deltaRFTrace(parseRecords.nGoodrec, 2, nPad);
	
	
	// Use target depth to find layer within the modelStack where migration should begin
	// This layer will then become the bottomStack for migration.
	MigrationParams velModel(velocityfn);
	double timeShiftRight = 0.0;
	double minTshft = 0.0;
	double maxTshft = 0.0;
	double vertTshft = 0.0;
	
	switch (phaseFlag) {
		case PPSMS:
			vertTshft = velModel.getVertTimeDelay(targetDepthKm, PPSMS);
			break;
		case PPPMS:
			vertTshft = velModel.getVertTimeDelay(targetDepthKm, PPPMS);
			break;
		default:
			break;
	}
	
	// timeshiftRight: The event window is moved(moving window migration - MWM)
	// depending on:
	// 1. The time scaling implied by the rayParameter
	// 2. The timedelay for each velocity layer and
	
	
	// Compute timeShiftRight here...
	for (int irec =0; irec < parseRecords.nGoodrec; irec++) {
		
		
		Sacread record =
		Sacread(parseRecords.goodRecordname[irec], headerTag, lqtRotArg);
		
		// HACK - Nov. 10, 2013 . Included by OlugbojiTolulope  -
		// Japanese Data requires that predicted time changes a lot, so we need to shift
		// timing for first arrival for each record ....
		tNoiseTrace = record.timeTag - record.timeStart - 3; //bias it by 3 sec
		nNoiseTrace = (tNoiseTrace) / (record.deltaT);
		
		double rayParam = double(record.raySlowness);				// in sec/Km
		
		
		// Get time shift ... If ray is trapped then remover record ...
		// If looking for time delay for vertically impiging ray,
		// Use rayParam equals 1. since cos zero = 1.
		switch (phaseFlag) {
			case PPSMS:
				timeShiftRight = velModel.getTimeDelayPpSms(targetDepthKm, rayParam);
				cout << "Time Delay PpSms: " << timeShiftRight << "s" << endl;
				break;
			case PPPMS:
				timeShiftRight = velModel.getTimeDelayPpPms(targetDepthKm, rayParam);
				cout << "Time Delay PpPms: " << timeShiftRight << "s" << endl;
				break;
			default:
				break;
		}
		
		// Update minimum and maximum time shifts.
		if (irec == 0 ) {
			minTshft = timeShiftRight;
			maxTshft = timeShiftRight;
		} else {
			if (timeShiftRight < minTshft) {
				minTshft = timeShiftRight;
			} else if (timeShiftRight > maxTshft) {
				maxTshft = timeShiftRight;
			}
		}
		
		// Noise Is Noise. Just Update Noise ...
		for (int cmp = 0; cmp < 3; cmp++) {
			
			for (int iter=0; iter < nPad; iter++) {
				if (iter < nNoiseTrace) {
					noiseTrace[irec][cmp][iter] = record.DataCmps[cmp][nNoiseTrace - iter];
				} else {
					noiseTrace[irec][cmp][iter] = 0.0;
				}
				
			}
			
			// Data is Migrated Using timeShift ... Only for R, T. Skip shift for Z.
			int nShift = (timeShiftRight) / (record.deltaT);
			if (cmp == 0) {
				nShift = 0;
			}
			int posEventNxt, posEventStrt = nNoiseTrace + nShift;
			
			for (int iter=0; iter < nPad; iter++) {
				posEventNxt = posEventStrt + iter;
				if ( iter < nEventTrace) {
					postEventTrace[irec][cmp][iter] = record.DataCmps[cmp][posEventNxt];
				} else {
					postEventTrace[irec][cmp][iter] = 0.0;
				}
			}
			
		}
		
		parseRecords.setHeaders(record.recBaz, record.recEpic, record.evLon,
								record.evLat, record.evMag, irec);
		
		/*Pass Record to MTC driver, computes Spectral Estimates*/
		/*Args: Noise, Event, p, k, fmax */
		MTCDriver MTCVals(irec, noiseTrace, postEventTrace, 2.5, 3, Fcutoff, record.deltaT,nNoiseTrace,nEventTrace);
		
		if (getVerbose) cout << "Calculting RF for record:  " << irec+1 << endl;
		MTCVals.addRF(RFTrace, irec);
		MTCVals.addCoher(deltaRFTrace, irec);
		
		
		
		
		
		
	}
	
	minT = minTshft; maxT = maxTshft; vertT = vertTshft;
	
	if (getVerbose) cout << "No of Layers + Halfspace is = " << velModel.getNoLayers() + 1 << endl;
	
	if (getVerbose) cout << "Done Building RFs. Regressing Nboot: " << nBoot << " times"<< endl ;
	
	// Dump RFs into TraceStack ...
	
	string fNameUpdate(outfn);
	fNameUpdate.append("MigrateReverbHarmonic");
	
	outfn = const_cast<char*>( fNameUpdate.c_str() );
	
	// Dump RFs into Harmonic Stack ...
	{
		
		int flagFile = CMPSHARMONIC;
		
		int nRecs = parseRecords.nGoodrec, ncmps = 2, nFreq = nPad;
		HarmonicStack testHStack(nPad, ncmps, nRecs, parseRecords.Azim, parseRecords.Epic, Fcutoff, dT, 10.0,  100.0, nBoot);
		testHStack.regressNBootTimes(RFTrace, deltaRFTrace);
		testHStack.print(outfn, flagFile, 0.0);
		//logFile << testHStack.logReport;
	}
	
	
}

//7. Printer for Epicentral Distance vs. Pred Phase time dumped to text file...
void RFrouter::callPhasePrinter(char* outfn, RecordList& parseRecords){
	cout << "Target Depth: " << targetDepthKm << endl;
	
	int nRecs = parseRecords.nGoodrec;
	vector<double> Epic, Azim, RayP;  // Record metaData for computing timing!
	Epic.assign(nRecs, 0.0); Azim.assign(nRecs, 0.0); RayP.assign(nRecs, 0.0);
	
	// Load Records Once and obtain local copy of metaData
	for (int irec =0; irec < nRecs; irec++) {
		
		Sacread record =
		Sacread(parseRecords.goodRecordname[irec], headerTag, lqtRotArg);
		
		double rayParam = double(record.raySlowness);				// in sec/Km
		
		Azim[irec] = record.recBaz;
		Epic[irec] = record.recEpic;
		RayP[irec] = rayParam;
		
		
		
	}
	
	
	// Logic to sort records into epicentral bins and aggregate timing ..
	/* Parameters to track if RFs are in wedge and stack 'em *****************/
	float dist2Bin;  // Distance between current distance and bin edge
	float episcan;
	float epi2bin;
	int minBinSze = binsize;
	
	
	episcan = epicmax;
	if (epicinc > 0.0) epicinc = -1.0 * epicinc;
	
	//Determine the stack size using bin parameters
	int stackSze = 0;
	while (episcan >= epicmin) {
		stackSze++;
		episcan += epicinc;
	}
	// Initialize Timing store - Bin Epic Distance, timing PS, timing PMS, timing SMS
	MatDoub timeEpic(stackSze,4,0.0);
	
	// reset the scan variable
	episcan = epicmax;
	int iterStck = 0, wdgeCnt = 0;
	bool isInWedge, isInBin, parseOnce = false;
	/************************************************************************/
	
	
	
	// Use target depth to find layer within the modelStack where migration should begin
	// This layer will then become the bottomStack for migration.
	MigrationParams velModel(velocityfn);
	double timeShiftRight = 0.0;
	double minTshft = 0.0;
	double maxTshft = 0.0;
	double vertTshft = velModel.getVertTimeDelay(targetDepthKm, PS);
	
	
	while (episcan >= epicmin) {
		int binCnt = 0; // No of traces in bin .. Min no for bin hit is in minBinsze
		
		
		// Run through all Records And Pick drop value in Epicentral bin ...
		for (int irec = 0; irec < Azim.size() ; irec++) {
			bool isInWedge = (Azim[irec]>= epbazmin && Azim[irec] <= epbazmax);
			dist2Bin = abs( Epic[irec] - episcan );
			bool isInBin = dist2Bin <= abs(epicinc);
			
			if( isInWedge && !parseOnce){
				wdgeCnt++;
			} 
			
			// Verify presence in wedge and in bin and then stack! dependent on status of flasgs. Review or debug
			if (isInWedge && isInBin) {
				binCnt++;
				
				timeEpic[iterStck][0] = episcan;
				timeEpic[iterStck][1] += velModel.getTimeDelayPs(targetDepthKm, RayP[irec]); 
				timeEpic[iterStck][2] += velModel.getTimeDelayPpPms(targetDepthKm, RayP[irec]);
				timeEpic[iterStck][3] += velModel.getTimeDelayPpSms(targetDepthKm, RayP[irec]);
				// timing PS..  how about others???
			
				
			}
		}
		/***********************  END Bin Scan ******************************/
		parseOnce = true;
		
		// Average out the timing for the bin. if no data in bin, then set zero.
		if (binCnt >= minBinSze) {
			timeEpic[iterStck][0] = episcan;
			timeEpic[iterStck][1] = timeEpic[iterStck][1] / binCnt; 
			timeEpic[iterStck][2] = timeEpic[iterStck][2] / binCnt;
			timeEpic[iterStck][3] = timeEpic[iterStck][3] / binCnt;

		}else {
			timeEpic[iterStck][0] = episcan;
			timeEpic[iterStck][1] = 0;
			timeEpic[iterStck][2] = 0;
			timeEpic[iterStck][3] = 0;

		}

		episcan += epicinc;
		iterStck++;
	}
	
	// Dump data in text file ...
	ofstream timedelayFile(outfn);
	
	for (int iterStck = 0; iterStck < stackSze; iterStck++) {
		timedelayFile << timeEpic[iterStck][0] << "   "
					<< timeEpic[iterStck][1] << "   "
					<< timeEpic[iterStck][2] << "   "
					<< timeEpic[iterStck][3] << "   " << endl;
	}

	
	

	

	
}

//8. Run reverse time to depth time migration. Implement in migration module
MatDoub RFrouter::migrateTime2Depth(MatDoub parseTime, double vertTshft){
	int nrows = parseTime.nrows(); int ncols = parseTime.ncols();
	MatDoub DepthMig; DepthMig.assign(nrows, ncols, 0.0);
	double ithMigTime = 0.0;
	MigrationParams velModel(velocityfn);
	
	// run through first row vector for TIME.
	for (int ithCol = 0; ithCol < ncols; ithCol++) {
		
		// Migrate only positive time values ...
		ithMigTime = parseTime[0][ithCol] + vertTshft;
		if ( ithMigTime >= 0.0 ) {
			
			cout << "T:  " << ithMigTime << "s  ";
			DepthMig[0][ithCol] = velModel.getVertDepthPs(ithMigTime, PS) / 1000.0;
			
			// Update all Rows ...
			for (int ithRow = 1; ithRow < nrows; ithRow++) {
				DepthMig[ithRow][ithCol] = DepthMig[0][ithCol];
			}
	
			cout << "D:  " << DepthMig[0][ithCol] / 1000.0 << "km  " << endl;
			
		}else{
			cout << "*********** Less than zero: " << ithCol << endl;
		}
	}
	cout << endl;
	return DepthMig;
	
}


//P1 Print Earthquake statistics here
void RFrouter::printEQStats(char* outfn, string evntStats){
	
	string fNameUpdate(outfn);
	fNameUpdate.append("_EQStats.txt");
	
	// Dump stats formatted as string into text file ...
	ofstream eqstatsOut(fNameUpdate);
	
	eqstatsOut << evntStats << endl;
	
}

//P2 Print Earthquake statistics here
void RFrouter::printDataTraces(char* outfn, int nRec, Mat3DDoub &traces, Mat3DDoub &traces2){

	// print raw
	string fNameUpdate;
	for (int iCmp = 0; iCmp < 3; iCmp++) {
		switch (iCmp) {
			case 0:
				// Dump stats formatted as string into text file ...
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_WvfrmRw_Z.txt");
				break;
			case 1:
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_WvfrmRw_R.txt");
				break;
			case 2:
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_WvfrmRw_T.txt");
				break;
			default:
				break;
		}
		
		ofstream wvfrmOut(fNameUpdate);
		
		//cout << "Debug timing" << nNoiseTrace+nEventTrace << endl;
		
		for (int irec = 0; irec < nRec; irec++) {
			float timeNxt = 0;
			
			for (int iNext = 0; iNext < nNoiseTrace+nEventTrace; iNext++) {
				
				wvfrmOut  << timeNxt << " \t "
				<< irec << "\t" << traces[irec][iCmp][iNext] << endl;
				
				timeNxt += dT;
			}
			wvfrmOut << ">" << endl;
		}
		
	}
	
	// print rotated ...
	for (int iCmp = 0; iCmp < 3; iCmp++) {
		switch (iCmp) {
			case 0:
				// Dump stats formatted as string into text file ...
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_WvfrmRot_Z.txt");
				break;
			case 1:
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_WvfrmRot_R.txt");
				break;
			case 2:
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_WvfrmRot_T.txt");
				break;
			default:
				break;
		}
		
		ofstream wvfrmOut(fNameUpdate);
		
		for (int irec = 0; irec < nRec; irec++) {
			float timeNxt = 0;
			
			for (int iNext = 0; iNext < nNoiseTrace+nEventTrace; iNext++) {
				
				wvfrmOut  << timeNxt << " \t " << irec << "\t"
				<< traces2[irec][iCmp][iNext] << endl;
				
				timeNxt += dT;
			}
			wvfrmOut << ">" << endl;
		}
		
	}
	

	
}

//P3 Print Earthquake statistics here
void RFrouter::printCoherSpec(char* outfn, int nRec, Mat3DDoub &traces, Mat3DDoub &traces2){
	
	int FreqNyq = nPad/2 + 1;
 	double FreqRlg = 1.0 / (dT * nPad);
	
	//int nFreqMax = (Fcutoff / FreqRlg) + 1;
	int nFreqMax = (Fcutoff *4.0 / FreqRlg) + 1; // Plot up to Nyquist
	int FreqSkip = 1;						// No skipping

	string fNameUpdate;
	
	// print  coherence
	for (int iCmp = 0; iCmp < 2; iCmp++) {
		switch (iCmp) {
			case 0:
				// Dump stats formatted as string into text file ...
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_Coher_ZR.txt");
				break;
			case 1:
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_Coher_ZT.txt");
				break;
			default:
				break;
		}
		
		ofstream wvfrmOut(fNameUpdate);
		
		for (int irec = 0; irec < nRec; irec++) {
			
			for (int iterF =0 ; iterF < nFreqMax; iterF+=FreqSkip) {
				
				wvfrmOut << iterF*FreqRlg << "\t"  << irec  << " \t "
				<< traces[irec][iCmp][iterF] << endl;
				
			}
			wvfrmOut << ">" << endl;
		}
		
	}
	
	// print all spectrum
	for (int iCmp = 0; iCmp < 6; iCmp++) {
		switch (iCmp) {
			case 0:
				// Dump stats formatted as string into text file ...
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_EvSpec_Z.txt");
				break;
			case 1:
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_EvSpec_R.txt");
				break;
			case 2:
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_EvSpec_T.txt");
				break;
			case 3:
				// Dump stats formatted as string into text file ...
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_NoiseSpec_Z.txt");
				break;
			case 4:
				fNameUpdate.assign(outfn);
				fNameUpdate.append("_NoiseSpec_R.txt");
				break;
			case 5:
				fNameUpdate.assign(outfn);				
				fNameUpdate.append("_NoiseSpec_T.txt");
				break;
			default:
				break;
		}
		
		ofstream wvfrmOut(fNameUpdate);
		
		for (int irec = 0; irec < nRec; irec++) {
			
			for (int iterF =0 ; iterF < nFreqMax; iterF+=FreqSkip) {
				
				wvfrmOut << iterF*FreqRlg << "\t" <<
				irec << "\t" << traces2[irec][iCmp][iterF] << endl;
				
			}
			wvfrmOut << ">" << endl;
		}
		
	}
	
}
