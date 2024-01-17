/**
 *  @file GrdSrchRouter.cpp
 *
 *  @author Tolulope  Olugboji
 *	@date  7/30/13.
 *  @copyright 2013 Yale University. All rights reserved.
 *
 *	@brief Performs Sequential-Multi H-K stacks.
 *
 *  This class reads the grdsearch paramter file, and determines the bounds
 *	for single (or multi-layer) properties - then uses the to call:
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

#include "GrdSrchRouter.h"
#include "interp_2d.h"


// Define PhaseFlags for Converted Vs. Reverberated Phases.
#define PS 0
#define PPSMS 1
#define PPPMS 2

using namespace std;


//HK Stacks are done in sector bins
//     Future work tests for dipping interfaces ...

GrdSrchRouter::GrdSrchRouter(char* eventsfn, char* outfn, float timeWin, float Fmax, float lqtArg, int hTag,  char* searchfn, int ROUTECODE, float epbzmin, float epbzmax, float epmin, float epmax,  bool getVbose)
:tEventTrace(timeWin), lqtRotArg(lqtArg), headerTag(hTag), Fcutoff(Fmax), searchParamsfn(searchfn), isAzim(false), isEpic(true), epbazmin(epbzmin), epbazmax(epbzmax), epicmin(epmin), epicmax(epmax), getVerbose(getVbose){
	
	//cout << "You just called the Epicentral Stacking Constructor. WIP" << endl;
	
	// Pick First Record and Calculate timing details.
	// Initialize these values in class object for router...
	runInterpolate = false;
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
	
	// Read  GRDsearch paramter file and determine what stack code to use
	loadGrdParams(searchParamsfn); 
	bool flgChkMigrate = false;
	stageRecords(allRecords, flgChkMigrate);
	
	// Call Appropriate GridSearchRouter
	Router(ROUTECODE, outfn, allRecords);
	
	// run nPoint by nPoint times and pick out elastic parameters iTimes
	/*
	for (int xDim = 0; xDim < nPoints; xDim++) {
		for (int yDim = 0; yDim < nPoints; yDim++) {
			
			double ithLayerZ  = hVals[xDim][yDim];
			double ithLayerVs = Vs;
			double ithLayerVp = Vs * kVals[xDim][yDim]; 
			
		}
	}
	 */
	
}

// Helper custructor for interpolation . Just one extra argument - decim
GrdSrchRouter::GrdSrchRouter(char* eventsfn, char* outfn, float timeWin, float Fmax, float lqtArg, int hTag,  char* searchfn, int ROUTECODE, float epbzmin, float epbzmax, float epmin, float epmax, int decim, bool getVbose)
:tEventTrace(timeWin), lqtRotArg(lqtArg), headerTag(hTag), Fcutoff(Fmax), searchParamsfn(searchfn), isAzim(false), isEpic(true), epbazmin(epbzmin), epbazmax(epbzmax), epicmin(epmin), epicmax(epmax), nPointsStretch(decim), getVerbose(getVbose){
	
	cout << "Constructor: Do similar data preparation before interpolation! nStretch:  " << nPointsStretch << endl;
	runInterpolate = true;
	
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
	
	// Read  GRDsearch paramter file and determine what stack code to use
	loadGrdParams(searchParamsfn);
	bool flgChkMigrate = false;
	stageRecords(allRecords, flgChkMigrate);
	
	// Call Appropriate GridSearchRouter
	 Router(ROUTECODE, outfn, allRecords);
	
}

void GrdSrchRouter::loadGrdParams(char* searchParamsfn){
	
	// Mimick Migration Constructor with minor diferences 
	// Here Read Parameters in the format: seqential layer bounds.
	
	std:string line, title;
	
	if (searchParamsfn != NULL) {
		ifstream stream(searchParamsfn, std::ios_base::in);
		
		if (!stream) {
			cerr << "in GrdSrchRouter::loadGrdParams: Cannot open ;" << searchParamsfn << "\n";
			exit(1);
		}
		
		getline(stream, title);
		stream >> noLayers;
		
		cout << "Title:  " << title << "No of Layers: "  << noLayers << endl;
		
		
		int iLayer = 1;
		while (!stream.eof()) {
			
			
			//stream >> iTheta >> iPhi;
			stream >> Hmin >> Hmax >> Kmin >> Kmax >> Vs >>  nPoints >> wtPs >> wtPms >> wtSms;
			

			cout << "Layer " << iLayer << endl;
			cout << " Thickness Bounds " << Hmin  << " " << Hmax << endl;
			//Update Layers ...
			layersHmin.push_back(Hmin); layersHmax.push_back(Hmax);
			
			cout << "Bounds on Vp/Vs Ratio " << Kmin  << " " << Kmax << endl;
			// update layers
			layersKmin.push_back(Kmin); layersKmax.push_back(Kmax);
			
			cout << "Vs value " << Vs << endl;
			// update layers
			layersVs.push_back(Vs);
			
			cout << "no. of GrdPoints: nPoints x nPoints: " << nPoints*nPoints << endl;
			// update layers
			layersNpoints.push_back(nPoints);
			
			cout  << "Phase weights used for composite stacking - sums to 1.";
			cout << "Weight for Ps phase: " << wtPs*100  << "%"
			<< " weight for PpPms phase: " << wtPms*100  << "%"
			<< " weight for PpSms phase: " << wtSms*100  << "%" << endl;
			layerWtPs.push_back(wtPs); layerWtPms.push_back(wtPms);
			layerWtSms.push_back(wtSms);
			
			// For each layer, do bounds check. Exit if bounds check fails!
			boundsCheck();
			
			if (iLayer >= noLayers) {
				break;
			}
			iLayer++;
		}
		stream.close();
		bool isRayGood = true;
	}else {
		cout << "No Grd Search Parameter File Found. -S option " << endl;
		cout << "Format: \n Line FOr File Description \n No layers \n thicknessMin thicknessMax kMin(VpVs ratio) kMax  Vs nPoints(grd size = nPnts x nPnts)" << endl;
		noLayers = -1;
	}
	
	
	
	
	// stackGrd stores the stack values ...
	
	
}

void GrdSrchRouter::stageRecords(RecordList& parseRecords, bool flgChkMigrate){
	
	// Added an option to check that records can be migrated: flgChkMigrate
	MigrationParams velModel(NULL);
	
	
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

void GrdSrchRouter::Router(int ROUTECODE, char* outfn, RecordList& allRecords){
	switch (ROUTECODE) {
		case HKSTCKSINGLE:
		{
			cout << "Using single layer H-K Stacking ..." << endl;
			
			doSectorGrdStack(outfn, allRecords, PS);
			doSectorGrdStack(outfn, allRecords, PPSMS);
			doSectorGrdStack(outfn, allRecords, PPPMS);
			
			saveGrdSumStack(outfn, 0.5, 0.4, 0.1);
			saveModelPred(outfn);
			
			break;
		}
		case HKSTCKMULTIPLE:
		{
			cout << "testing the multiple H-K Stacking incrementally..." << endl;
			if (noLayers == 1) {
				cout << "Param file single layer" << endl;
				
				setLayerParams(ithGridParams, 1.0, 2.0, 3.0);
				cout << ithGridParams.H << "  " << ithGridParams.Vp << "  " << ithGridParams.Vs << endl;
			}else{
				cout << "Param file multi layer " << endl;
				

				
				for (int iterLayer = 0; iterLayer < noLayers; iterLayer++) {
					
					// Initialize current and previous layer
					noLayerAbove = iterLayer;
					nPoints = layersNpoints[iterLayer];
					Vs = layersVs[iterLayer];
					initializeGrd();
					
					// reset H-Kk search grids for the currentLayer ..
					make2DimGrd(layersHmin[iterLayer], layersHmax[iterLayer], 
								layersKmin[iterLayer], layersKmax[iterLayer], layersNpoints[iterLayer],  hVals, kVals);
					
					//printMatrix(hVals, string("Thickness hVals"));
					//printMatrix(kVals, string("VpVs ratio kVals"));
					//printMatrix(stackGrdPPSMS, string("S(H,K)  stackGrdPS"));
					
					
					if (iterLayer > 0 ) {
						
						//Vs = layersVs[iterLayer-1];
						
						setLayerParams(ithSolvedParams,hPred[iterLayer-1],kPred[iterLayer-1]*layersVs[iterLayer-1], layersVs[iterLayer-1]);
						//setLayerParams(ithSolvedParams, iterLayer,iterLayer+1, iterLayer+2);
						updateLayerAbove(solvedLayerAbove, ithSolvedParams);
					}
					
					
					printLayerAbove(solvedLayerAbove);
					
					// Update layer name before save ...
					std::ostringstream ss;
					ss << (iterLayer + 1);
					
					string fNameUpdate(outfn); 
					fNameUpdate.append( "_Layer_" + ss.str() );
					
					char * outfnLayer = const_cast<char*>( fNameUpdate.c_str() );
					cout << outfnLayer << endl;
					
					
					
					doSectorGrdStack(outfn, allRecords, PS);
					doSectorGrdStack(outfn, allRecords, PPSMS);
					doSectorGrdStack(outfn, allRecords, PPPMS);
					
					double wt1 = layerWtPs[iterLayer] ;
					double wt2 = layerWtSms[iterLayer] ;
					double wt3 = layerWtPms[iterLayer];
					saveGrdSumStack(outfnLayer, wt1, wt2, wt3);
					printMatrix(vertTGrdPS, string("Timing for PS, After LookUp"));
				}
				 
			}
			
			// Save multi-layer model ...
			saveModelPred(outfn);
			
			break;
		}
		case HKSTCKMULTIPLEINTERPOLATE:
		{
			cout << "Adds interpolation functionality to stack matrices" << endl;
			if (noLayers == 1) {
				// Dummy code for testing dataStruc. Parameters.
				cout << "Param file single layer" << endl;
				
				setLayerParams(ithGridParams, 1.0, 2.0, 3.0);
				cout << ithGridParams.H << "  " << ithGridParams.Vp << "  " << ithGridParams.Vs << endl;
			}else{
				cout << "Param file multi layer with interpolation " << endl;
				
				for (int iterLayer = 0; iterLayer < noLayers; iterLayer++) {
					
					// Initialize current and previous layer
					noLayerAbove = iterLayer;
					nPoints = layersNpoints[iterLayer];
					Vs = layersVs[iterLayer];
					initializeGrd();
					
					// reset H-Kk search grids for the currentLayer ..
					make2DimGrd(layersHmin[iterLayer], layersHmax[iterLayer],
								layersKmin[iterLayer], layersKmax[iterLayer], layersNpoints[iterLayer],  hVals, kVals);
					
					//printMatrix(hVals, string("Thickness hVals"));
					//printMatrix(kVals, string("VpVs ratio kVals"));
					//printMatrix(stackGrdPPSMS, string("S(H,K)  stackGrdPS"));
					
					
					if (iterLayer > 0 ) {
						
						//Vs = layersVs[iterLayer-1];
						
						setLayerParams(ithSolvedParams,hPred[iterLayer-1],kPred[iterLayer-1]*layersVs[iterLayer-1], layersVs[iterLayer-1]);
						//setLayerParams(ithSolvedParams, iterLayer,iterLayer+1, iterLayer+2);
						updateLayerAbove(solvedLayerAbove, ithSolvedParams);
					}
					
					
					printLayerAbove(solvedLayerAbove);
					
					// Update layer name before save ...
					std::ostringstream ss;
					ss << (iterLayer + 1);
					
					string fNameUpdate(outfn);
					fNameUpdate.append( "_Layer_" + ss.str() );
					
					char * outfnLayer = const_cast<char*>( fNameUpdate.c_str() );
					cout << outfnLayer << endl;
					
					
					
					doSectorGrdStack(outfn, allRecords, PS);
					doSectorGrdStack(outfn, allRecords, PPSMS);
					doSectorGrdStack(outfn, allRecords, PPPMS);
					
					// **** Since you save here. Interpolate before saving files..
					double wt1 = layerWtPs[iterLayer] ;
					double wt2 = layerWtPms[iterLayer];
					double wt3 = layerWtSms[iterLayer] ;
					
					saveGrdSumStack(outfnLayer, wt1, wt2, wt3);
					printMatrix(vertTGrdPS, string("Timing for PS, After LookUp"));
				}
				
			}
			
			// Save multi-layer model ...
			saveModelPred(outfn);
			
			break;
		}
		default:
			break;
	}
}

void GrdSrchRouter::doSectorGrdStack(char* outfn, RecordList& parseRecords, int phaseFlag){
	
	
	Mat3DDoub noiseTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DDoub postEventTrace(parseRecords.nGoodrec, 3, nPad);
	Mat3DCmplx RFTrace(parseRecords.nGoodrec, 2, nPad);
	Mat3DDoub deltaRFTrace(parseRecords.nGoodrec, 2, nPad);
	
	
	double timeShiftRight = 0.0;
	double minTshft = 0.0;
	double maxTshft = 0.0;
	double vertTshft = 0.0;
	
	// **  2 loops for Grid Iteration. 1 Loop For Record and Migration 
	// Another single loop for  Layered Structure..
	
	// ? Use the grid panel to load the target depth and velocity migration parameters.
	// run nPoint by nPoint times and pick out elastic parameters iTimes -
	// in test, set nPoint = 1
	
	// Status update for grid search ...	
	int totPoints = nPoints * nPoints;
	int iPoints = 1;
	
	for (int xDim = 0; xDim < nPoints; xDim++) {
		for (int yDim = 0; yDim < nPoints; yDim++) {
			
			double ithLayerZ  = hVals[xDim][yDim];
			double ithLayerVs = layersVs[noLayerAbove];
			double VpVsRatio = kVals[xDim][yDim];
			double ithLayerVp = layersVs[noLayerAbove] * VpVsRatio;
			
			if (getVerbose) {
				cout << "Iterating through Grid at: " << ithLayerZ << "km, Vs " << ithLayerVs << " km/s, Vp" << ithLayerVp << "km/s" << endl;
			}
			
			
			setLayerParams(ithGridParams, ithLayerZ, ithLayerVp, ithLayerVs);
			MigrationParams ithVelModel(ithGridParams, solvedLayerAbove, noLayerAbove );
			
			//MigrationParams ithVelModel(ithLayerZ, ithLayerVp, ithLayerVs, 0.0);
			
			
			// Pick out vertical time shift for display and debug purposes
			// the timing information will be stored  and saved for each stackPhase 
			switch (phaseFlag) {
				case PS:
					vertTshft = ithVelModel.getVertTimeDelay(ithLayerZ/1000, PS);
					break;
				case PPSMS:
					vertTshft = ithVelModel.getVertTimeDelay(ithLayerZ/1000, PPSMS);
					break;
				case PPPMS:
					vertTshft = ithVelModel.getVertTimeDelay(ithLayerZ/1000, PPPMS);
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
				
				double rayParam = double(record.raySlowness);				// in sec/Km
				
				
				// Get time shift ... If ray is trapped then remover record ...
				// If looking for time delay for vertically impiging ray,
				// Use rayParam equals 1. since cos zero = 1.
				switch (phaseFlag) {
					case PS:
						timeShiftRight = ithVelModel.getTimeDelayPs(ithLayerZ/1000, rayParam);
						//cout << "Time Delay PS: " << timeShiftRight << "s" << endl;
						break;
					case PPSMS:
						timeShiftRight = ithVelModel.getTimeDelayPpSms(ithLayerZ/1000, rayParam);
						//cout << "Time Delay PpSms: " << timeShiftRight << "s" << endl;
						break;
					case PPPMS:
						timeShiftRight = ithVelModel.getTimeDelayPpPms(ithLayerZ/1000, rayParam);
						//cout << "Time Delay PpPms: " << timeShiftRight << "s" << endl;
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
			
			/* Done Calculating RFs. Stack and Store */
			{
				minT = minTshft; maxT = maxTshft; vertT = vertTshft;
				
				if (getVerbose) cout << "No of Layers + Halfspace is = " << ithVelModel.getNoLayers() + 1 << endl;
				
				//cout << "Done Building RFs, for: " << parseRecords.nGoodrec << "recs . Grid search at " << iPoints << " of " << totPoints  << endl ;
				
				// Write Code in TraceStack For Summary stacks....
				int preTime = 2, postTime = 2;
				
				int flagFile = STACKGRD;				// Used in printing - rewrite for GRDSTACK.
				
				int binsize = 1;
				TraceStack epicStack(RFTrace, deltaRFTrace, parseRecords.Azim, parseRecords.Epic,
									 Fcutoff, dT, preTime, postTime, binsize);
				
				double epicinc = epicmax - epicmin + 1.0;			// Force single stack ...
				double HKStackval = epicStack.grdStack(epbazmin, epbazmax, epicmin, epicmax, epicinc,RFTrace, deltaRFTrace);
				
				//save stack and timing information depending on phase value 
				switch (phaseFlag) {
					case PS:
						stackGrdPS[xDim][yDim] = HKStackval;
						vertTGrdPS[xDim][yDim] = vertTshft;
						break;
					case PPSMS:
						stackGrdPPSMS[xDim][yDim] = HKStackval;
						vertTGrdPPSMS[xDim][yDim] = vertTshft;
						break;
					case PPPMS:
						stackGrdPPPMS[xDim][yDim] = HKStackval;
						vertTGrdPPPMS[xDim][yDim] = vertTshft;
						break;
					default:
						break;
				}
				
			}
			
			iPoints = iPoints + 1;
		}
		
	}
	
	
	switch (phaseFlag) {
		case PS:
		{
			// Update layer name before save ...
			std::ostringstream ss;
			ss << (noLayerAbove + 1);
			string fNameUpdate(outfn); 
			char * outfnLayer = const_cast<char*>( fNameUpdate.c_str() );
			
			fNameUpdate.append("_Layer_" + ss.str() + "_HKGrdStackPS.txt");
			outfn = const_cast<char*>( fNameUpdate.c_str() );
			saveGrid(outfn, stackGrdPS);
			break;
		}
		case PPSMS:
		{
			// Update layer name before save ...
			std::ostringstream ss;
			ss << (noLayerAbove + 1);
			string fNameUpdate(outfn); 
			char * outfnLayer = const_cast<char*>( fNameUpdate.c_str() );
			
			fNameUpdate.append("_Layer_" + ss.str() + "_HKGrdStackPPSMS.txt");
			outfn = const_cast<char*>( fNameUpdate.c_str() );
			saveGrid(outfn, stackGrdPPSMS);
			break;
		}
		case PPPMS:
		{
			// Update layer name before save ...
			std::ostringstream ss;
			ss << (noLayerAbove + 1);
			string fNameUpdate(outfn); 
			char * outfnLayer = const_cast<char*>( fNameUpdate.c_str() );
			
			fNameUpdate.append("_Layer_" + ss.str() + "_HKGrdStackPPPMS.txt");
			outfn = const_cast<char*>( fNameUpdate.c_str() );
			saveGrid(outfn, stackGrdPPPMS);
			break;
		}
		default:
			break;
	}
	
	
	
	
	
	
	//cout << epicStack.logReport << endl;
	// test the max function ...
	
}

void GrdSrchRouter::initializeGrd(){
	
	// Computed Matrices ...
	stackGrdPS = MatDoub(nPoints, nPoints, 0.0);
	stackGrdPPSMS = MatDoub(nPoints, nPoints, 0.0);
	stackGrdPPPMS = MatDoub(nPoints, nPoints, 0.0);
	
	stackGrdFULL = MatDoub(nPoints,nPoints, 0.0);
	
	vertTGrdPS = MatDoub(nPoints, nPoints, 0.0);
	vertTGrdPPSMS = MatDoub(nPoints, nPoints, 0.0);
	vertTGrdPPPMS = MatDoub(nPoints, nPoints, 0.0);
	
	// Resmampled grid using interpolation ...
	if (runInterpolate == true) {
		newNpoints = nPoints*nPointsStretch;
		
		stackGrdPSresmple = MatDoub(newNpoints, newNpoints, 0.0);
		stackGrdPPSMSresmple = MatDoub(newNpoints, newNpoints, 0.0);
		stackGrdPPPMSresmple = MatDoub(newNpoints, newNpoints, 0.0);
		
		stackGrdFULLresmple = MatDoub(newNpoints,newNpoints, 0.0);
		
		vertTGrdPSresmple = MatDoub(newNpoints, newNpoints, 0.0);
		vertTGrdPPSMSresmple = MatDoub(newNpoints, newNpoints, 0.0);
		vertTGrdPPPMSresmple = MatDoub(newNpoints, newNpoints, 0.0);
		
		// New vectors for interpolation
		 newHVals =  VecDoub(newNpoints);
		 newKVals = VecDoub(newNpoints);
	}

}

void GrdSrchRouter::updateLayerAbove(vecLayers& multiLayers, Layer prevSolvedLayer){
	
	multiLayers.push_back(prevSolvedLayer);
	
}

void GrdSrchRouter::printLayerAbove(vecLayers& LayerAbove){
	if (LayerAbove.empty()) {
		cout << "Currently Empty! No Layer Above" << endl;
	}else{
		
		cout << "Solution to solved layers are:" << endl;
		cout << "Thickness  Vp     Vs      "<< endl;
		for (int iLayer = 0; iLayer < LayerAbove.size(); iLayer++) {
		cout << LayerAbove[iLayer].H << "  " << LayerAbove[iLayer].Vp << "  " << LayerAbove[iLayer].Vs << endl;
		}
		
	}
		
}

void GrdSrchRouter::setLayerParams(Layer& sLayer, double h, double vp, double vs){
	sLayer.H = h;
	sLayer.Vp = vp;
	sLayer.Vs = vs;
}

void GrdSrchRouter::make2DimGrd(double Xmin, double Xmax, double Ymin, double Ymax, int nPoints, MatDoub& xGrd, MatDoub& yGrd){
	
	// Declare S(H,K) Panel to Store Maximum Values? Replaced in initializeGrd
	//MatDoub XYGrid(nPoints, nPoints, 0.0);
	MatDoub YVals(nPoints, nPoints, 0.0);
	MatDoub XVals(nPoints, nPoints, 0.0);
	
	double Ystep = (Ymax - Ymin) / (nPoints - 1 );
    double	Xstep = (Xmax - Xmin) / (nPoints - 1 );
	
	//Initialize YVals Panels For Iteration ..
	for (int xDim = 0; xDim < nPoints; xDim++ ) {
		for (int yDim = 0; yDim < nPoints; yDim++) {
			
			XVals[xDim][yDim] = Xmin + Xstep*yDim;
			
		}
	
	}
	
	cout << "~~~~~~~~~~~~~~~~~~~" << endl;
	//Initialize XVals Panels For Iteration ..
	for (int yDim = 0; yDim < nPoints; yDim++ ) {
		for (int xDim = 0; xDim < nPoints; xDim++) {
			
			YVals[xDim][yDim] = Ymin + Ystep*xDim;
		
		}
		
	}
	
	//printMatrix(XVals);
	//printMatrix(YVals);
	
	
	//outGrd = XYGrid;
	xGrd = XVals;
	yGrd = YVals;
	
}

void GrdSrchRouter::printMatrix(MatDoub& inMat, string matName){

	
	cout << "\n\n printMatrix for Debug.." << matName << endl;
	for (int xDim = 0; xDim <inMat.nrows(); xDim++ ) {
		cout << "[ " ;
		for (int yDim = 0; yDim < inMat.ncols(); yDim++) {
			
			cout << inMat[xDim][yDim] << "  ";
		}
		cout << "] " << endl;
	}
	
}

void GrdSrchRouter::lkUpStck4PrdVals(){

	double findH = 0, findK = 0;
	int findXindx = 0, findYindx=0;
	
	// Discriminate between resampled stack data or true stack data.
	
	if (runInterpolate == TRUE) {
		
		double maxStckVal = stackGrdFULLresmple[0][0];
		// local store of predicted layer values: H and K values.
		
		for (int xDim = 0; xDim < newNpoints; xDim++) {
			for (int yDim = 0; yDim < newNpoints; yDim++) {
				
				if (stackGrdFULLresmple[xDim][yDim] >= maxStckVal) {
					maxStckVal = stackGrdFULLresmple[xDim][yDim];
					findXindx = xDim; findYindx = yDim;
				}
			}
		}
		
		// Store predicted values from inversion
		hPred.push_back(newHVals[findYindx]);
		kPred.push_back(newKVals[findXindx]);
		
		// predicted timing for the phases ... This files are also stored?
		tPredPS.push_back(vertTGrdPSresmple[findXindx][findYindx]);
		tPredPPSMS.push_back(vertTGrdPPSMSresmple[findXindx][findYindx]);
		tPredPPPMS.push_back(vertTGrdPPPMSresmple[findXindx][findYindx]);
	}else {
		
		double maxStckVal = stackGrdFULL[0][0];
		// local store of predicted layer values: H and K values.
		
		for (int xDim = 0; xDim < nPoints; xDim++) {
			for (int yDim = 1; yDim < nPoints; yDim++) {
				
				if (stackGrdFULL[xDim][yDim] >= maxStckVal) {
					maxStckVal = stackGrdFULL[xDim][yDim];
					findXindx = xDim; findYindx = yDim;
				}
			}
		}
		
		// Store predicted values from inversion
		hPred.push_back(hVals[findXindx][findYindx]);
		kPred.push_back(kVals[findXindx][findYindx]);
		
		// predicted timing for the phases ...
		tPredPS.push_back(vertTGrdPS[findXindx][findYindx]);
		tPredPPSMS.push_back(vertTGrdPPSMS[findXindx][findYindx]);
		tPredPPPMS.push_back(vertTGrdPPPMS[findXindx][findYindx]);
		
	}
	


	
}


// Extending this routine with interpolation
void GrdSrchRouter::saveGrdSumStack(char* outfn, double w1, double w2, double w3){
	
	
	double stPs, stSMS, stPMS;
	
	// normalize all matrices
	normMatMax(stackGrdPS);
	normMatMax(stackGrdPPPMS);
	normMatAbsMax(stackGrdPPSMS);
			   
	
	// Stack already conducted ... Interpolate before summing !
	if (runInterpolate == true) {
		
		// Stack Then Interpolate. Interpolate then Stack?? Any difference?
		// Define new dimension and initialize with vectorDoub - used by interp2d.h

		
		double hMin = hVals[0][0];
		double hMax = hVals[0][hVals.nrows()-1];
		
		double kMin = kVals[0][0];
		double kMax = kVals[kVals.ncols()-1][0];
		
		double hStep = (hMax - hMin) / (newNpoints - 1);
		double kStep = (kMax - kMin) / (newNpoints - 1);
		
		cout << "hMin:  " << hMin << "hMax:  " << hMax << endl;
		cout << "kMin:  " << kMin << "kMax:  " << kMax << endl;
		
	
		for (int iterNpoints = 0; iterNpoints < newNpoints; iterNpoints++) {
			newHVals[iterNpoints] = hMin + hStep*iterNpoints;
			newKVals[iterNpoints] = kMin + kStep*iterNpoints;
		}
		
		// display new values ...
		for (int iterNpoints = 0; iterNpoints < newNpoints; iterNpoints++) {
			cout << newHVals[iterNpoints] << "  ";
			//cout << newKVals[iterNpoints] << "  ";
		}
		cout << endl;
		
		// display new values ...
		for (int iterNpoints = 0; iterNpoints < newNpoints; iterNpoints++) {
			//cout << newHVals[iterNpoints] << "  ";
			cout << newKVals[iterNpoints] << "  ";
		}
		cout << endl;
		
		cout << "New grid dimensions: " << newNpoints << endl;
		
		// Check normalization of matrices ...
		cout << "PS:  " << valMax(stackGrdPS) << "   PPSMS:   " << valAbsMax(stackGrdPPSMS) << "  PPPMS:   " << valMax(stackGrdPPPMS) << endl;
		
		// Conduct bilinear or cubic spline interpolation here ...
		// I need a new data structure or data structures with the right dimension
		// and then spit out the data. Don't forget LAB migration
		VecDoub hVec(nPoints); VecDoub kVec(nPoints); int vecSze = hVals.nrows();
		for (int iterVec = 0; iterVec < vecSze; iterVec++) {
			hVec[iterVec] = hVals[0][iterVec];
			kVec[iterVec] = kVals[iterVec][0];
		}
		
		// RF AMPLITUDE stacks!
		Spline2D_interp PSfunc(hVec, kVec, stackGrdPS);
		Spline2D_interp SMSfunc(hVec, kVec, stackGrdPPSMS);
		Spline2D_interp PMSfunc(hVec, kVec, stackGrdPPPMS);
		
		//Bilin_interp PSfunc(hVec, kVec, stackGrdPS);
		//Bilin_interp SMSfunc(hVec, kVec, stackGrdPPSMS);
		//Bilin_interp PMSfunc(hVec, kVec, stackGrdPPPMS);
		
		// RF TIMING for the stacks!
		Spline2D_interp TPSfunc(hVec, kVec, vertTGrdPS);
		Spline2D_interp TSMSfunc(hVec, kVec, vertTGrdPPSMS);
		Spline2D_interp TPMSfunc(hVec, kVec, vertTGrdPPPMS);
		
		//Bilin_interp TPSfunc(hVec, kVec, vertTGrdPS);
		//Bilin_interp TSMSfunc(hVec, kVec, vertTGrdPPSMS);
		//Bilin_interp TPMSfunc(hVec, kVec, vertTGrdPPPMS);
		
		
		int grdSze = stackGrdPSresmple.nrows();
		
		cout << "grdSze: " << grdSze << endl;
		
		for (int xDim = 0; xDim < grdSze; xDim++) {
			for (int yDim = 0; yDim < grdSze; yDim++) {
				
				// resample RF AMPLITUDE stacks!
				stackGrdPSresmple[xDim][yDim] = PSfunc.interp(newHVals[xDim], newKVals[yDim]);
				stackGrdPPSMSresmple[xDim][yDim] = SMSfunc.interp(newHVals[xDim], newKVals[yDim]);
				stackGrdPPPMSresmple[xDim][yDim] = PMSfunc.interp(newHVals[xDim], newKVals[yDim]);
				
				// resample RF TIMING for the stacks!
				vertTGrdPSresmple[xDim][yDim] = TPSfunc.interp(newHVals[xDim], newKVals[yDim]);
				vertTGrdPPSMSresmple[xDim][yDim] = TSMSfunc.interp(newHVals[xDim], newKVals[yDim]);
				vertTGrdPPPMSresmple[xDim][yDim] = TPMSfunc.interp(newHVals[xDim], newKVals[yDim]);
				
			}
		}
		
		
		
		// Stack After interpolation
		for (int xDim = 0; xDim < grdSze; xDim++) {
			for (int yDim = 0; yDim < grdSze; yDim++) {
				
				
				stPs = stackGrdPSresmple[xDim][yDim];
				stSMS = stackGrdPPSMSresmple[xDim][yDim];
				stPMS =  stackGrdPPPMSresmple[xDim][yDim];
				
				stackGrdFULLresmple[xDim][yDim] = (w1 * stPs) + (w2* stPMS) + (-1.0*w3*stSMS );
				
			}
		}
		 
		
		// Look up predictions here. As in Max Values.
		lkUpStck4PrdVals();
		
		// Save full grd stack ....
		string fNameUpdate(outfn); string fNameUpdate2(outfn);
		fNameUpdate.append("_Int_HKGrdStackFULL.txt");
		outfn = const_cast<char*>( fNameUpdate.c_str() );
		saveGrid(outfn, stackGrdFULLresmple);
		
		// Save interpolated stacks... All 3 of them
		for (int iterPhase = 0; iterPhase < 3; iterPhase++) {
			string outfnLayer(fNameUpdate2);
			
			switch (iterPhase) {
				case PS:
				{
					// Update layer name before save ...
					std::ostringstream ss;
					ss << (noLayerAbove + 1);
					
					outfnLayer.append("_Int_HKGrdStackPS.txt");
					outfn = const_cast<char*>( outfnLayer.c_str() );
					saveGrid(outfn, stackGrdPSresmple);
					break;
				}
				case PPSMS:
				{
					// Update layer name before save ...
					std::ostringstream ss;
					ss << (noLayerAbove + 1);
					
					outfnLayer.append("_Int_HKGrdStackPPSMS.txt");
					outfn = const_cast<char*>( outfnLayer.c_str() );
					saveGrid(outfn, stackGrdPPSMSresmple);
					break;
				}
				case PPPMS:
				{
					// Update layer name before save ...
					std::ostringstream ss;
					ss << (noLayerAbove + 1);
					string fNameUpdate(outfn);

					outfnLayer.append("_Int_HKGrdStackPPPMS.txt");
					outfn = const_cast<char*>( outfnLayer.c_str() );
					saveGrid(outfn, stackGrdPPPMSresmple);
					break;
				}
				default:
					break;
			}
		}
		
		
	}else{
		
		// Stack without interpolating ...
		for (int xDim = 0; xDim < nPoints; xDim++) {
			for (int yDim = 0; yDim < nPoints; yDim++) {
				
				stPs = stackGrdPS[xDim][yDim];
				stSMS = stackGrdPPSMS[xDim][yDim];
				stPMS =  stackGrdPPPMS[xDim][yDim];
				
				stackGrdFULL[xDim][yDim] = (w1 * stPs) + (w2* stPMS) + ( 1.0*w3*stSMS );
				
			}
		}
		
		// Look up predictions here. As in Max Values.
		lkUpStck4PrdVals();
		
		// Save full grd stack ....
		string fNameUpdate(outfn);
		fNameUpdate.append("HKGrdStackFULL.txt");
		outfn = const_cast<char*>( fNameUpdate.c_str() );
		saveGrid(outfn, stackGrdFULL);
	}
	

	
	
}

void GrdSrchRouter::saveModelPred(char* outfn){
	
	double 	density = 2700;		//Kg per mcubed, review this with lookup table.
	double  Vp  = 0.0;			// Calculated using Vp-Vs ratio...
	
	// Initialize model name for ParkLevin Style ...
	string fNameUpdate(outfn);
	fNameUpdate.append("_IsoHKmodels.txt");
	ofstream modelOut( fNameUpdate.c_str() );
	

	
	// First Line is the description of model ...
	modelOut << "Model produced from RFVelStck [HK Stacking] (c)Olugboji" << endl;
	modelOut << noLayers << endl;				//noLayers - updated for MutliHK
	
	for (int iLayer = 0; iLayer < noLayers;  iLayer++) {
		
		
		modelOut << "0    0" << endl;		    //Isotropic - dip & tilt
		
		Vp = layersVs[iLayer] * kPred[iLayer];
		modelOut << hPred[iLayer] << "  " << Vp << " 0   0 "  << layersVs[iLayer] << "  0  " << density << endl;
		
		
	}
	// Mantle HalfSpace Underneath ...
	
	modelOut << "0    0" << endl;		    //Isotropic - dip & tilt
	modelOut << 80000 << "  " << 8160 << "  " <<  0 <<  "  " << 0  << "  "  <<  4750 << "  " <<  0  << "  "  << 30000 << endl;
	
	
	// Dump timing information as well ...
	fNameUpdate.append("_PhaseTiming.txt");
	ofstream timingOut( fNameUpdate.c_str() );
	
	// Dump timing information ... Update for multilayers...
	timingOut << "Timing info from RFVelStck [HK Stacking] (c)Olugboji" << endl;

	for (int iLayer = 0; iLayer < noLayers; iLayer++) {
		
		timingOut << "Layer: " << iLayer <<  endl;
		timingOut << "Predicted PS Timing (ray param = 0): " << tPredPS[iLayer] << endl;
		timingOut << "Predicted PPSMS Timing (ray param = 0): " << tPredPPSMS[iLayer] << endl;
		timingOut << "Predicted PPPMS Timing (ray param = 0): " << tPredPPPMS[iLayer] << endl;
		timingOut << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ " <<  endl;
		
	}

	
}

void GrdSrchRouter::saveGrid(char* outfn, MatDoub& stackGrd){
	
	//Prep file name for file save
	string fileHkGrd(outfn);
	ofstream grdOut(fileHkGrd.c_str());
	
	int nSze = stackGrd.nrows();
	
	// Do a dimensionality test here. If interpolation implied, use different save
	
	if (nSze > nPoints) {
		
		for (int xDim = 0; xDim < nSze; xDim++) {
			for (int yDim = 0; yDim < nSze; yDim++) {
				grdOut << newHVals[yDim]/1000 << "  " << newKVals[xDim] << "  " <<stackGrd[xDim][yDim] << endl;
			}
			//grdOut << ">" << endl;
		}
	}else{
		
		for (int xDim = 0; xDim < nSze; xDim++) {
			for (int yDim = 0; yDim < nSze; yDim++) {
				grdOut << hVals[xDim][yDim]/1000 << "  " << kVals[xDim][yDim] << "  " <<stackGrd[xDim][yDim] << endl;
			}
			//grdOut << ">" << endl;
		}
	}

	
}

double GrdSrchRouter::valMax(MatDoub& inGrd){
	
	// Store administrative details of matrix in here.
	int nSze = inGrd.nrows();
	double tempMax = inGrd[0][0];
	
	
	for (int xDim = 0; xDim < nSze; xDim++) {
		for (int yDim = 0; yDim < nSze; yDim++) {
			
			if (tempMax < inGrd[xDim][yDim]) {
				tempMax = inGrd[xDim][yDim];
			}

			
		}
	}
	
	return tempMax;
	
}

double GrdSrchRouter::valAbsMax(MatDoub& inGrd){
	
	// Store administrative details of matrix in here.
	int nSze = inGrd.nrows();
	double tempMax = abs(inGrd[0][0]);
	
	
	for (int xDim = 0; xDim < nSze; xDim++) {
		for (int yDim = 0; yDim < nSze; yDim++) {
			
			if (tempMax < abs(inGrd[xDim][yDim]) ) {
				tempMax = abs(inGrd[xDim][yDim]);
			}
			
			
		}
	}
	
	return tempMax;
	
}

void GrdSrchRouter::normMatMax(MatDoub& inGrd){
	
	// Store administrative details of matrix in here.
	int nSze = inGrd.nrows();
	double maxVal = valMax(inGrd);
	double tempVal = 0.0;
	
	
	for (int xDim = 0; xDim < nSze; xDim++) {
		for (int yDim = 0; yDim < nSze; yDim++) {
			
			tempVal = inGrd[xDim][yDim];
			inGrd[xDim][yDim] = tempVal / maxVal;
		}
	}
	
	
}

void GrdSrchRouter::normMatAbsMax(MatDoub& inGrd){
	
	// Store administrative details of matrix in here.
	int nSze = inGrd.nrows();
	double maxVal = valAbsMax(inGrd);
	double tempVal = 0.0;
	
	
	for (int xDim = 0; xDim < nSze; xDim++) {
		for (int yDim = 0; yDim < nSze; yDim++) {
			
			tempVal = inGrd[xDim][yDim];
			inGrd[xDim][yDim] = tempVal / maxVal;
		}
	}

	
}

void GrdSrchRouter::boundsCheck(){
	bool validateParams = false;
	bool validateThick = false;
	bool validateVpVsRatio = false;
	bool validateVs = false;
	bool validateWeights = false;
	string failMsg("");
	
	//1.  bounds check on thickness
	if(layersHmax.size() > 1) {
		int indxLast = layersHmax.size() - 2;
		cout << " Indx last " << indxLast << " "  << layersHmax[indxLast] <<  "  " << Hmin << endl;
		//
		if ( (layersHmax[indxLast] <= Hmin) && (Hmax > Hmin) ){
			validateThick = true;
		}else{
			failMsg.append("previous layer overlaps current layer or layer bounds are not open!");
		}
	}else{
		
		if ( Hmax > Hmin ){
			validateThick = true;
		}else{
			failMsg.append("bounds on layer thickness not open!");
		}
		
	}
	
	//2. bounds check on vp/vs ratio
	if (Kmax > Kmin) {
		validateVpVsRatio = true;
	}else{
		failMsg.append("bounds on Vp/Vs ratio not open");
	}
	
	
	//3. bounds check on vs
	if(layersVs.size() > 1) {
		int indxLast = layersVs.size() - 2;
		
		//
		if ( (layersVs[indxLast] <= Vs) && (Vs > 0.0) ){
			validateVs = true;
		}else{
			failMsg.append("Expected: crustal velocity should increase with depth? and Vs greater than zero ");
		}
	}else{
		
		if ( Vs > 0 ){
			validateVs = true;
		}else{
			failMsg.append("Vs should be a positive number");
		}
		
	}
	
	//4. bounds check on weights for the phase stack
	if ( (wtPs + wtSms + wtPms) == 1.0 ) {
		validateWeights = true;
	}else{
		failMsg.append("Phase weights must sum to 1");
	}
	
	
	validateParams = validateThick && validateVpVsRatio && validateWeights && validateVs;
	
	if (!validateParams) {
		cout << "\n\n\n *****@@@@@@!!!!!***@@@@@**!Warning! Grid Search Parameters Not Well Formed. \n ";
		cout << failMsg << endl;
		cout << "File format below: \n";
		cout << "File Title \nnLayers \n[hmin, hmax, kmin, kmax, nPts, wt1, wt2, wt3] \n.... " << endl;
		exit(1);
	}
	
}


