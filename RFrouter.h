/**
 *	@file RFrouter.h
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
 *  @bug no known bugs
*/         

#ifndef __RFHarmonicStacking__RFrouter__
#define __RFHarmonicStacking__RFrouter__

#include <iostream>
#include <assert.h>
#include "RecordList.h"
#include "sacread.h"
#include "TraceStack.h"
#include "HarmonicStack.h"
#include "MigrationParams.h"
#include "MTCDriver.h"


// Harmonic DisplayTimes
#define  CMPSHARMONIC 0
#define  EXPANSIONHARMONIC 1

// Router Specific Codes ....

/** Simple Azimuthal- Epicentral Stacks without migration */
#define SIMPLEAZIMEPIC 1
/** Simple Harmonic Stacks without migration */
#define SIMPLEHARMONIC 2
/** Azimuthal- Epicentral Stacks with moving window migration */
#define MWMAZIEPIC 3
/** Harmonic Stacks with moving window migration */
#define MWMHARMONIC 4
/** Azimuthal-Epicentral stacks with frequency migration */
#define FREQMIGRATEAZIMEPIC 5
/** Harmonic Stacks with frequency migration */
#define FREQMIGRATEHARMONIC 6
/** Simple Azimuthal- Epicentral Stacks with event cross correlation */
#define ECCAZIMEPIC 7
/** Simple Harmonic Stacks with event cross correlation */
#define ECCHARMONIC 8

/** Migrate Reverberated PpSms */
#define MOHOREVERBPSMS 9	
/** Migrate Reverberated PpPms. */
#define MOHOREVERBPPMS 10
/** Migrate Reverberated PpSms + Harmonic Stack */
#define HARMONICREVERBPSMS 11  
/** Migrate Reverberated PpPms + Harmonic Stack */
#define HARMONICREVERBPPMS 12  

/** Skip RF computation and print timing for all phases at target depth */
#define PRINTPHASETIME 13

class RFrouter{
public:
	// Simplest Constructor: 11 Parameters: Minimum No. Of Parameters Needed.
	// Add 1 more parameter - horizontal rotation ... after lqt rotation
	RFrouter(char* eventsfn, char* outfn,float timeWin, float Fcutoff, float lqtArg, float horArg, int hTag, int stckArg, int migArg, char* velfn, int routCode, bool getVbose);
	
	// Constructor for Azim-Stacks. --- WIProgress...
	// Added arguments to collect stacking parameters
	// Azimuth stack - Decimation parameters: 
	// Epicentral stack - Sector wedge and decimation parameters
	
	// Azimuth Stacks: 11 + 4 = 15 Parameters [bazmin, bazmax, bazinc, binsize]
	// Add 1 more parameter - horizontal rotation ... after lqt rotation
	RFrouter(char* eventsfn, char* outfn,float timeWin, float Fcutoff, float lqtArg, float horArg, int hTag, int stckArg, int migArg, char* velfn, int routCode, float bazmin, float bazmax, float bazinc, int bsze, bool getVbose);
	
	// Epicentral Stacks: 11 + 6 = 17 Parameters [ epbazmin, epbazmax, epicmin, epicmax, epicinc, binsize]
	// Add 1 more parameter - horizontal rotation ... after lqt rotation
	RFrouter(char* eventsfn, char* outfn,float timeWin, float Fcutoff, float lqtArg, float horArg, int hTag, int stckArg, int migArg, char* velfn, int routCode, float epbazmin, float epbazmax, float epicmin, float epicmax, float epicinc, int bsze, bool getVbose);
	
	
	
	
private:
	
	// Data Dependent Variables Tells How To Manipulate
	// SAC Data. And How Much Data to Pick... 
	float dT, tEventTrace, tNoiseTrace;
	int nNoiseTrace, nEventTrace, nMax, nPad;
	
	
	// Determines How to Stage Records. Header Tag
	// Specifies which record to keep which to discard
	int headerTag;
	
	// Determine to Rotate or Not Rotate SAC
	float lqtRotArg;
	
	// Determine to DO horizontal rotations
	float horRotArg;
	
	
	// Miscellaneous
	bool getVerbose;
	int nBoot;
	float Fcutoff;
	
	
	// Determine if stack type: Azim or Epic.
	bool isAzim, isEpic;
	
	// Parameter for Azimuth and Epicentral sweeps
	float bazmin, bazmax, bazinc;
	float epbazmin, epbazmax, epicmin, epicmax, epicinc;
	int binsize;
	
	// Migration: FileName 4 velocity & targetDepth in Km
	// NULL if not in migration mode.
	char* velocityfn;
	int targetDepthKm;
	double displayTimeShift, minT, maxT, vertT;
	
	
	
	
	// For Single Event Cross Correlation .. Parameters Here
	// 1 pi Taper. Bandwidth-Frequency Product Here ...
	
	void Router(int rtcode, char* outfn, RecordList& rlist);
	
	void stageRecords(RecordList& parseRecords, bool flgChkMigrate);
	
	// Listed in order of implementation...
	//1.
	void callHarmonicStack(char* outfn, RecordList& parseRecords);
	
	//2.
	void callAzimEpicStack(char* outfn, RecordList& parseRecords);
	
	//3. 
	void callAzimEpicMWMStack(char* outfn, RecordList& parseRecords);
	
	//4. Code to Route for Reverberated Phase Migration
	void callAzimEpicMOHOREVERBstack(char* outfn, RecordList& parseRecords, int phaseFlag);
	
	//5. Migrated Harmonic Stacks...
	void callHarmonicMWMStack(char* outfn, RecordList& parseRecords);
	
	//6. Code to Route for Reverberated Phase Migration
	void callHarmonicReverbStack(char* outfn, RecordList& parseRecords, int phaseFlag);
	
	
	//7. Printer for Epicentral Distance vs. Pred Phase time dumped to text file...
	void callPhasePrinter(char* outfn, RecordList& parseRecords);
	
	//8. Run reverse time to depth time migration. Implement in migration module
	MatDoub migrateTime2Depth(MatDoub parseTime, double vertTshft);
	
	//9. Helper module todo iterated search for true horizontal ldirections ...
	/* Rotate horizontals by single angle - clockwise or anticlockwise, then spit out harmonics  */
	void singleRotateComputeHarmonic(float rotateAngle);
	
	//10. 
	void callAzimEpicECCStack(char* outfn, RecordList& parseRecords);
	
	/* Do multiple rotation and spit out a grid estimate of the rms constant unmodelled amplitudes */
	void scanRotateComputeHarmonic(float nStartAngle, float nStopAngle, float nStepAngle);
	
	
	// New print functions
	
	// P1. Print event statistics
	void printEQStats(char* fName, string evntStats);
	
	// P2. Print traces before and after rotation
	void printDataTraces(char*fName, int nRec, Mat3DDoub &traces, Mat3DDoub &traces2);
	
	// P3. Print spectrum before and after rotation
	void printCoherSpec(char*fName, int nRec, Mat3DDoub &traces, Mat3DDoub &traces2);

	// P4. Print coherence before and after rotation
	
	// string concatenation class..
	template <class T>
	inline std::string to_string (const T& t)
	{
		std::stringstream ss;
		ss << t;
		return ss.str();
	}
	
};



#endif /* defined(__RFHarmonicStacking__RFrouter__) */
