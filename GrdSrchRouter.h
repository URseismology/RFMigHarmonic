/**  
 *  @file GrdSrchRouter.h
 *  
 *  @author Tolulope  Olugboji 
 *	@date  7/30/13.
 *  @copyright 2013 Yale University. All rights reserved.
 *
 *	@brief Performs Sequential-Multi H-K stacks.
 *
 *  This class reads the grdsearch paramter file, and determines the bounds
 *	for single (or multi-layer) properties - then uses this to call:
 *				1. MigrationParams - Constructor for single layer migration
 *				2. TraceStack - Member function to return zero time amplitude
 *
 *	Before code does this. It uses paramter file to set up grid... then loops
 *  grid...
 *
 *  Most other concepts are borrowed from RFrouter -
 *			e.g. 1.RecordStaging
 *				 2. ROUTECODES - Type of HK stacking etc.
 *  .. Add a module in trace stack: 
			[TraceStack.grdStck()] that returns a zero time amplitude for S(h,k).
 *  @warning Implementation in progress
 *  @bug Work in progress ... log bug reports here
 */



#ifndef RFHarmonicStacking_GrdSrchRouter_h
#define RFHarmonicStacking_GrdSrchRouter_h


#include <iostream>
#include <assert.h>
#include "RecordList.h"
#include "sacread.h"
#include "TraceStack.h"   
#include "MigrationParams.h"
#include "MTCDriver.h"


/** 1. HK Stack Single ... */
#define HKSTCKSINGLE 1	
/** HK Stack Multiple - in devpt phase*/
#define HKSTCKMULTIPLE 2
/** HK Stack Mutliple with Spline Interpolation - in devpt phase */
#define HKSTCKMULTIPLEINTERPOLATE 3

/** Recursive HK in Freq Domain - Joint algo with JJpark */
#define HKSTCKRECURSEFREQ 4


class GrdSrchRouter{
public:
	
	/** Main Constructor 
	 * Epicentral Summary Stacks: 13 Parameters [ 8 + 4(define sector) + 1(routecode)]
	 */
	GrdSrchRouter(char* eventsfn, char* outfn,float timeWin, float Fcutoff, float lqtArg, int hTag, char* searchParamsFn, int routCode, float epbazmin, float epbazmax, float epicmin, float epicmax, bool getVbose);
	
	/** Support Constructor
	 * Add interpolate arg: "nMainConstructor args + 1" = 14 Parameters
	 */
	GrdSrchRouter(char* eventsfn, char* outfn,float timeWin, float Fcutoff, float lqtArg, int hTag, char* searchParamsFn, int routCode, float epbazmin, float epbazmax, float epicmin, float epicmax, int decim, bool getVbose);
	
	
	/**
	 @brief select from a list of routines to conduct grid stack.
	 @param outfn file name for output files
	 @param RecordList object for storing list of SAC records
	 */
	void Router(int rtcode, char* outfn, RecordList& rlist);
	
	/**
	 @brief load or unload SAC files based on QC flag check
	 @param parseRecords list of SAC files 
	 @param flgChkMigrate set to skip waves trapped in layer
	 */
	void stageRecords(RecordList& parseRecords, bool flgChkMigrate);
	
	/**
	 @brief initialize Grid parameters using searchParams file
	 @param searchParams name of file to load search parameters
	 */
	void loadGrdParams(char* searchParams);
	
	/**
	 @brief initialize Grid values after loading Parameter files
	 */
	void initializeGrd();
	
	/**
	 @brief initialize ithGrid or ithSolved Layer Parameter files  ...
	 */
	void setLayerParams(Layer& sLayer, double H, double Vp, double Vs);
	
	/**
	 @brief update preceeding layers with solved HK values for multi-layer stacking
	 */
	void updateLayerAbove(vecLayers& multiLayers, Layer ithSolved);
	
	/**
	 @brief print Layer parameter in file ...
	 @param mutl
	 */
	void printLayerAbove(vecLayers& LayerAbove);
	
	/**
	 @brief initialize the 2D grids
	 @param Xmin minimum Xvalue
	 @param Xmax maximum Xvalue
	 @param Ymin minimum Yvalue
	 @param Ymax maximum Yvalue
	 @param nDivs number of points in X and Y dimension
	 @return outGrd 2D grid for stack values - PS
	 @return xGrd 2D grid for Xvalue
	 @return yGrd 2D grid for Yvalue
	 */
	void make2DimGrd(double Xmin, double Xmax, double Ymin, double Ymax, int nDivs, MatDoub& xGrd, MatDoub& yGrd);
	
	/**
	 @brief compute 2D grid stack.
	 @param outfn file name for grid output
	 @param parseRecords list of SAC files that pass QC check
	 @param phaseFlag MACRO code to distinguish PS, PPPMS, or PPSMS phase
	 */
	void doSectorGrdStack(char* outfn, RecordList& parseRecords, int phaseFlag);
	
		
	/**
	 @brief save singlePhase grd to output file in xyz GMT format
	 in comparison, saveGrdSumStack uses weighting parameters to construct final
	 summary grid and then saves to the file ...
	 @param outfn file name for grid file
	 @param outGrd grid file to be saved
	 */
	void saveGrid(char* outfn, MatDoub& outGrd);
	
	/**
	 @brief look up sumary stack for predicted values - used ONLY by saveGrdSumStack
	 @return - stores hPred( hVal (indxMax) ), kPred, and tPred in member variables
	 @bugs no bugs detected ...
	 */
	void lkUpStck4PrdVals();
	
	
	/**
	 @brief saves the full Grd Stack with weighting parameters (w1 + w2 + w3 = 1)
	 
	
	 uses lkUpStck4PrdVals to find predicted values using indx of maximum of grdStck
	 
	 @param w1 weighting for phase PS
	 @param w2 weighting for phase PPSMS
	 @param w3 weighting for phase PPPMS
	 @bugs no bugs deteced
	 */
	void saveGrdSumStack(char* outfn, double w1, double w2, double w3);
	
	/**
	 @brief saves the predicted values for hPred, kPred in ParkLevinStyle outfn_"HKFwd.txt"
	 
	 saves the timing of the phase in outfn_"Timing" .txt
	 
	 @param outfn file name in the -O tag...
	 @warning working on updating this for multi-layer models.
	 */
	void saveModelPred(char* outfn);
	
	
	/**
	 @brief print matrix in readable format
	 @param inMat Matrix
	 @param matName name of matrix
	 @return indirect return of grid values
	 */
	void printMatrix(MatDoub& inMat, string matName);
	
	/**
	 @brief return maximum value in 2D array...
	 @param inMat Matrix
	 @return maximum value ...
	 */
	double valMax(MatDoub& inMat);
	
	/**
	 @brief return absolute maximum value in 2D array...
	 @param inMat Matrix
	 @return maximum value ...
	 */
	double valAbsMax(MatDoub& inMat);
	
	
	/**
	 @brief normalize matrix by max value.
	 @param inMat Matrix
	 @return maximum value ...
	 */
	void normMatMax(MatDoub& inMat);
	
	/**
	 @brief normalize matrix by min value.
	 @param inMat Matrix
	 @return maximum value ...
	 */
	void normMatAbsMax(MatDoub& inMat);
	

	/**
	 @brief perform bounds check on search params.
	 @param Hmin minimum thickness for layer
	 @param Hmax maximum thickness for layer
	 @param Kmin minimum vp/vs ratio for layer
	 @param Kmax maximum vp/vs ratio for layer
	 @param Vs s velocity for layer
	 @param nPoints no. of grid points for layer
	 @param wtPs weighting for Ps phase in layer
	 @param wtPms weighting for PpPms phase in layer
	 @param wtSms weighting for PpSms phase in layer
	 */
	void boundsCheck();
	
	
	
private:
	
	/** Data Dependent Variables Tells How To Manipulate
	 * SAC Data. And How Much Data to Pick... 
	 */
	float dT, tEventTrace, tNoiseTrace;
	int nNoiseTrace, nEventTrace, nMax, nPad;
	
	
	// Determines How to Stage Records. Header Tag
	// Specifies which record to keep which to discard
	int headerTag;
	
	// Determine to Rotate or Not Rotate SAC
	float lqtRotArg;
	
	
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
	char* searchParamsfn; int noLayers;
	int targetDepthKm;
	double displayTimeShift, minT, maxT, vertT;
	
	
	// Grid DataValues initialized by loadGrdParams
	MatDoub stackGrdPS, stackGrdPPSMS, stackGrdPPPMS;		// Vertical time for each Grid
	MatDoub stackGrdPSresmple, stackGrdPPSMSresmple, stackGrdPPPMSresmple;
	
	MatDoub vertTGrdPS, vertTGrdPPSMS, vertTGrdPPPMS;		// Vertical time for each Grid
	MatDoub vertTGrdPSresmple, vertTGrdPPSMSresmple, vertTGrdPPPMSresmple;		// Vertical time for each Grid
	
	MatDoub stackGrdFULL, hVals, kVals;		// GrdStack Full & H,K dimension
	VecDoub newHVals, newKVals;
	MatDoub stackGrdFULLresmple;
	vector<double> hPred, kPred, tPredPS, tPredPPSMS, tPredPPPMS;	 // Use Stack to find predicted H,K, and Phase timing... all 3 main phases.
	
	
	// Grid point values for the ith layer
	double Hmin, Hmax, Kmin, Kmax, Vs, rho;
	
	// Phase weights for the ith Layer - must sum to 1.
	double wtPs, wtPms, wtSms;
	
	int nPoints, newNpoints;
	int nPointsStretch; // provided by user - used for higher decimation before interpolation ...
	bool runInterpolate;
	
	// Grid point values for each layer ...
	vector<double> layersHmin, layersHmax, layersKmin, layersKmax, layersVs;
	// Phase weights for each layer
	vector<double> layerWtPs, layerWtPms, layerWtSms;
	vector<int> layersNpoints ;
	
	Layer ithGridParams, ithSolvedParams;
	vecLayers solvedLayerAbove;
	int noLayerAbove;
	
		

	
	
};


#endif
