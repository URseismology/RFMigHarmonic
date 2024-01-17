/**
 *  @file TraceStack.h
 *
 *
 *  @author Tolulope Olugboji
 *	@date 12/26/10.
 *  @copyright Yale University. All rights reserved.
 *
 *  @brief This class is a specialist class that stacks the Reciever functions using the coherence estimates.
 */

#include "TraceStack.h"
#include "MTCDriver.h"
#include <string>

TraceStack::TraceStack(Mat3DCmplx &RF, Mat3DDoub &delta, 
					   VecDoub &A, VecDoub &E, float f, float delT, float pre, 
					   float post, int bSze)
: Azim(A), Epic(E), binmax(355.0), binmin(0.0),
 bininc(-5.0), zero(0.0,0.0), FreqMax(f), DeltaT(delT), preTime(pre), postTime(post), outLen(0),
binSze(bSze), flagDepthMig(false), setPrintStats(false)
{
	
	/* Initialize your 3Dimensional arrays. Extension: 
	 Write data structure driver in nr3.h
	 Here I allocate the same memory structure ..
	 Solution Point: Tentative problem: I can't initialize NRMat3D look into
	 this. For now put this initialization in stackAzim routine*/
	Mat3DCmplx RFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	Mat3DDoub sigmaRFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	
	
	/* Deep Copy Here*/
	for (int first = 0; first < RF.dim1(); first++) {
		for (int second = 0; second < RF.dim2(); second++) {
			for (int third = 0; third < RF.dim3(); third++) {
				RFTrace[first][second][third] = RF[first][second][third];
				sigmaRFTrace[first][second][third] = pow(delta[first][second][third],2.0);
			}
		}
	}
	// Done initializing the 3Dimensional arrays.
	
	int recLen = RFTrace.dim3();
	FreqNyq = recLen/2 + 1;
	FreqRlg = 1.0 / (DeltaT * recLen);
	nFreqMax = (FreqMax / FreqRlg) + 1;
	Pih = M_PI/2.0;
	iZero = 0;
	nPad = recLen;
	padNyq = nPad/2 + 1;
	//outStream << iZero << endl;
	
	logReport.append("#No of records stacked: ");
	int nRecords = Azim.size();
	logReport.append( to_string( nRecords ) );
	logReport.append( "\n" );
}

/* Fills stackRFTrace, stacAzim & calls the inverseFFT, okay consult
   MTCDriver to see how to do this properly
*/
void TraceStack::stackAzim(float max, float min, float inc, 
						   Mat3DCmplx &RF, Mat3DDoub &delta) {
	
	logReport.append("Azim Stack Called. Summarry by Bin \n ");
		
	/* Put sections here into my TraceStack class*/
	
	/* Initialize your 3Dimensional arrays. Extension: 
	 Write data structure driver in nr3.h
	 Here I allocate the same memory structure ..
	 Solution Point: Tentative problem: I can't initialize NRMat3D look into
	 this. For now put this initialization in stackAzim routine*/
	Mat3DCmplx RFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	Mat3DDoub sigmaRFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	
	
	/* Deep Copy Here, Dump file here so that you can plot [Debug Purpose]*/
	 
	for (int first = 0; first < RF.dim1(); first++) {
		int cnt = 0;
		for (int third = 0; third < RF.dim3(); third++) {
			for (int second = 0; second < RF.dim2(); second++) {
				RFTrace[first][second][third] = RF[first][second][third];
				sigmaRFTrace[first][second][third] = pow(delta[first][second][third],2.0);
				
			}
		}
	}
	// Done initializing the 3Dimensional arrays.
	float bazscan, binscan;
	float Azim2bin;			// Distance between Azim & Current Bin Boundary
	int minBinSize = binSze;		// Min. No. of traces in bin during Stack
	
	binmin = min;
	binmax = max;
	bininc = inc;
	
	bazscan = binmax;
	if (bininc > 0.0) bininc = -1.0 * bininc;
	binscan = binmax;
	
	/* Determine stack size using bin parameters ...
	*/
	int stackSze = 0;
	while (bazscan >= binmin) {
		stackSze++;
		bazscan += bininc;
	}
	bazscan = binmax;
	
	/* Initialize 3D Stack Matrix, note that some of the bins might not contain any traces,
	   So I use another vector 'azimStack' to check if a bin within the stack contains traces:
	 To do this I initialize and fill the azimStack matrix: 
	 column1: Left Azim boundary in Bin 
	 column2: status of 'bin sum', -1.0 if a bin is not filled with traces; 1.0 if it is.
	 */
	Mat3DCmplx stackRFTrace(stackSze, RFTrace.dim2(), nPad); //[dim1: bin][dim2: rad,trans][dim3:freq]
	for (int first = 0; first < stackSze; first++) {
		for (int third = 0; third < nPad; third++) {
			for (int second = 0; second < 2; second++) {
				stackRFTrace[first][second][third] = zero;
			}
		}
	}
	
	
	azimStack.assign(stackSze, 2, -1.0);
	
	
	
	int iterStck = 0;
	
	while (bazscan >= binmin) {
		
		// Variable binCnt holds number of traces within a bin
		int binCnt = 0; // No of Traces in bin .. Min. no for a BinHit determined by minBinSze
		// Dimensions of RF in dim2 & nPad important for inverse fourier
		
		MatComplex StackRFHolder(RFTrace.dim2(), nPad, zero);  
		MatDoub StackSigmaHolder(RFTrace.dim2(), nPad, iZero);
		
		// Run through entire Azim. Scan here
		for (int irec = 0; irec < Azim.size(); irec++) {
			Azim2bin = abs(Azim[irec] - bazscan);
			
			// Verify presence in bin & stack all RFs in this bin.
			if ( Azim2bin < abs(bininc) || Azim2bin > 360.0 - abs(bininc) ) {
				binCnt++;
				
				//Stack Reciever Functions using their coherence estimates
				for (int itercmp = 0; itercmp < 2; itercmp++) {
					for (int iternf = 0; iternf < nFreqMax; iternf++) { 
						StackRFHolder[itercmp][iternf] += 
								RFTrace[irec][itercmp][iternf] / sigmaRFTrace[irec][itercmp][iternf];
						StackSigmaHolder[itercmp][iternf] += 1.0 / sigmaRFTrace[irec][itercmp][iternf];
					}
				}
				
			}
		}
		
		// Push in new Stack Trace, depending on Scan ...
		if (binCnt > minBinSize) {
			for (int itercmp = 0; itercmp < 2; itercmp++) {
				for (int iternf = 0; iternf < nFreqMax; iternf++) {
					
					Hann = pow( cos( Pih * (double(iternf)/double(nFreqMax)) ), 2.0);
					stackRFTrace[iterStck][itercmp][iternf] = 
					Hann * StackRFHolder[itercmp][iternf] / StackSigmaHolder[itercmp][iternf];
					
				}
			}
			azimStack[iterStck][1] = 1.0;      // Particular Azimuth bin hit. Yay!
		}else {
			// Should I not surpress the entry in stak here, if I can't get any
			// Reciever functions within this azimuth bin?
			azimStack[iterStck][1] = -1.0;      // No Traces in this Azimuth bin. Aaaww..
		}

		
		azimStack[iterStck][0] = bazscan;
		/***********  LOG SUMMARY ************/
		logReport.append( "Azimuth: " );
		logReport.append( to_string(bazscan) );
		logReport.append( " , No of RFs in bin:" );
		logReport.append( to_string( binCnt ) );
		logReport.append( " \n");
		/************ END SUMMARY ************/
		
		bazscan += bininc;
		iterStck++;
	}
	
	/* Iterate through Stack and use MTCDriver to invert the stack. 
	   I put this outside Stack construction because decisions need to be 
	   made about whether to print an azimuth bin depending on its HIT status	  
	 */
	int npre = preTime / DeltaT, npost = postTime / DeltaT;
	outLen = npre + npost;                 // Variable holds output legth
	radialRFAscii.assign(stackSze, outLen, 0.0);
	transRFAscii.assign(stackSze, outLen, 0.0);
	timeAscii.assign(stackSze, outLen, 0.0);
	//depthAscii.assign(stackSze, outLen, 0.0);
	
	for (int iterStck = 0; iterStck < stackSze; iterStck++) {
		
		for (int itercmp = 0; itercmp < 2; itercmp++) {
			
			VecComplex RFfreq(stackRFTrace.dim3(), stackRFTrace[iterStck][itercmp]);
			VecDoub_IO RFtimebyAzim;
			MTCDriver invertStack;
			
			invertStack.icmplxfft(RFfreq, RFtimebyAzim, padNyq);
			
			// Dump RF here.. big bug. Yish!!
			/*for (int i = 0; i < RFtimebyAzim.size(); i++) {
				debugFile << abs(RFtimebyAzim[i]) << "  " << endl;
			} */
			//Build Reciever Function in time domain? 
			//Causative & Non Causative Build
			if (itercmp == 0) {
				buildAsciiOut(preTime, postTime, RFtimebyAzim, radialRFAscii, timeAscii,
							  iterStck);
			}else {
				buildAsciiOut(preTime, postTime, RFtimebyAzim, transRFAscii, timeAscii,
							  iterStck);
			}
			
		}
	}
	
	
}

/* This member method helps build the ascii output before using the out function 
   All RFs in this section are in the time domain, what I need to do is get more
   printfns to print RFs in the frequency domain for test cases ...
 */
void TraceStack::buildAsciiOut(double pretime, double posttime, VecDoub_IO &RFin, 
							  MatDoub &RFout, MatDoub &timeOut, int ithAzim)
{
	int recLen = RFin.size();
	int npre = pretime / DeltaT, npost = posttime / DeltaT;
 
	double scaleRF = 2.0 * FreqNyq / (FreqMax/FreqRlg);
	
	for (int i = 1; i <= (npre+npost) ; i++) {
		timeOut[ithAzim][i-1] =  -npre*DeltaT + (i-1)*DeltaT ;
	}
	
	for (int i =1 ; i <= npost; i++) {
		RFout[ithAzim][npre+i-1] = scaleRF * RFin[i-1];
	}
	
	for (int i = 1; i <= npre; i++) {
		RFout[ithAzim][npre-i] = scaleRF * RFin[recLen - i-1];
	}
	
		
}

/* This is  ahelper function that retrurns maximum radial amplitude and timing. 
 I will tweak it in future to just return the maximum at zero time - used by grdStack */
void TraceStack::findMAX(MatDoub &RFout, MatDoub &timeOut, double& maxOut, double& maxInd)
{
	int npre = preTime / DeltaT, npost = postTime / DeltaT;
	outLen = npre + npost;
	
	int stckSize = RFout.nrows();     // Make sure RF is a single vector.
	int RFsze = RFout.ncols();        // RFvector - Size .
	
	if (stckSize > 1) {
		cout << "Can't do summary stack, RF has more than 1 rows" << endl;
	}
	
	//cout << "*********** RF testing MAX function ..." <<  "  " << endl;
	for (int iterRF = 0; iterRF < (outLen); iterRF++) {
		
		if ( timeOut[0][iterRF] == 0.0 || abs(timeOut[0][iterRF] - 0.0) <= 0.5*DeltaT  ) {
			//cout << RFout[0][iterRF] << "  " << timeOut[0][iterRF] << endl;
			maxOut = RFout[0][iterRF]; maxInd = timeOut[0][iterRF];
			break;
		}
		
		
	}
	
}

/* Fills stackRFTrace, stacAzim & calls the inverseFFT, okay consult
 Instead of Azimuth on vertical axis, we have Epic Distance ...
 View this as a wedge cut out from the a section close to the station
 */
void TraceStack::stackEpic(float epbazmin, float epbazmax, float min, float max, float inc,
						   Mat3DCmplx &RF, Mat3DDoub &delta) {
	
	logReport.append("Epicentral Stack Called. \n Epicentral dist:  ");
	logReport.append("Baz min:  ");logReport.append(to_string(epbazmin)); logReport.append("\n  ");
	logReport.append("Baz max:  ");logReport.append(to_string(epbazmax)); logReport.append("\n  ");
	logReport.append("Epic min: "); logReport.append(to_string(min)); logReport.append("\n  ");
	logReport.append("Epic max:  "); logReport.append(to_string(max)); logReport.append("\n  ");
	
	/* Put sections here into my TraceStack class*/
	
	/* Initialize your 3Dimensional arrays. Extension: 
	 Write data structure driver in nr3.h
	 Here I allocate the same memory structure ..
	 Solution Point: Tentative problem: I can't initialize NRMat3D look into
	 this. For now put this initialization in stackAzim routine*/
	Mat3DCmplx RFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	Mat3DDoub sigmaRFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	
	
	/* Deep Copy Here */
	for (int first = 0; first < RF.dim1(); first++) {
		for (int third = 0; third < RF.dim3(); third++) {
			for (int second = 0; second < RF.dim2(); second++) {
				RFTrace[first][second][third] = RF[first][second][third];
				sigmaRFTrace[first][second][third] = pow(delta[first][second][third],2.0);
			}
		}
	}
	
	/* Parameters to track if RFs are in wedge and stack 'em */
	float dist2Bin;  // Distance between current distance and bin edge
	float episcan;
	float epi2bin;
	int minBinSze = binSze;
	
	binmin = min;
	binmax = max;
	bininc = inc;
	
	episcan = binmax;
	if (bininc > 0.0) bininc = -1.0 * bininc;
	
	//Determine the stack size using bin parameters
	int stackSze = 0;
	while (episcan >= binmin) {
		stackSze++;
		episcan += bininc;
	}
	// reset the scan variable 
	episcan = binmax;
	
	/* Initialize 3D Stack Matrix, note that some of the bins might not contain any traces,
	 So I use another vector 'epicStack' to check if a bin within the stack contains traces:
	 To do this I initialize and fill the azimStack matrix:
	 Column1: Left Epic boundary in Bin
	 Column2: status of 'bin sum', -1.0 if a bin is not filled with traces; 1.0 f it is
	 */
	Mat3DCmplx stackRFTrace(stackSze, RFTrace.dim2(), nPad); //[dim1: bin][dim2: rad,trans][dim3:freq]
	for (int first = 0; first < stackSze; first++) {
		for (int third = 0; third < nPad; third++) {
			for (int second = 0; second < 2; second++) {
				stackRFTrace[first][second][third] = zero;
			}
		}
	}	  // Initialize the Array here ... revisit nr3.h
	epicStack.assign(stackSze,2,-1.0); // Read comment above on why this is initialized to -1.0
	
	int iterStck = 0;
	int wdgeCnt = 0; // No of Traces that fall into sector wedge
	bool parseOnce = false;
	
	// Scan variable scans down, that's why I have it initialized to outer boundary of the wedge
	// This scan is currently inefficient. Revisit. Think of ways to make more efficient.
	while (episcan >= binmin) { 
		int binCnt = 0; // No of Traces in bin .. Min. no for a BinHit is determined by minBinSze
		
		/* The two variables below hold the stack within this bin. 
		Obviously they'll be empty if there is no RF in this Epicentral scan */
		MatComplex StackRFHolder(RFTrace.dim2(), nPad, zero);  
		MatDoub StackSigmaHolder(RFTrace.dim2(), nPad, iZero);
		
		// Run through all Records ...
		for (int irec = 0; irec < Azim.size() ; irec++) {
			bool isInWedge = (Azim[irec]>= epbazmin && Azim[irec] <= epbazmax);
			dist2Bin = abs( Epic[irec] - episcan );
			bool isInBin = dist2Bin <= abs(bininc);
			
			if( isInWedge && !parseOnce){
				wdgeCnt++;
			} 
				
			// Verify presence in wedge and in bin and then stack! dependent on status of flasgs. Review or debug
			if (isInWedge && isInBin) {
				binCnt++;
				// Stack RFs using their coherence estimates 
				for (int itercmp = 0; itercmp < 2; itercmp++) {
					for (int iternf = 0; iternf < nFreqMax; iternf++) { 
						StackRFHolder[itercmp][iternf] += 
						RFTrace[irec][itercmp][iternf] / sigmaRFTrace[irec][itercmp][iternf];
						StackSigmaHolder[itercmp][iternf] += 1.0 / sigmaRFTrace[irec][itercmp][iternf];
					}
				}
				
			}
		}
		/***********************  END Bin Scan ******************************/
		parseOnce = true;
		// Push in new Stack Trace, depending on Scan ...
		if (binCnt >= minBinSze) {
			for (int itercmp = 0; itercmp < 2; itercmp++) {
				for (int iternf = 0; iternf < nFreqMax; iternf++) {
					
					Hann = pow( cos( Pih * double(iternf/nFreqMax) ), 2.0);
					stackRFTrace[iterStck][itercmp][iternf] = 
					Hann * StackRFHolder[itercmp][iternf] / StackSigmaHolder[itercmp][iternf];
					
				}
			}
			epicStack[iterStck][1] = 1.0;      // Particular Azimuth bin hit. Yay!
		}else {
			// Should I not surpress the entry in stak here, if I can't get any
			// Reciever functions within this azimuth bin?
			epicStack[iterStck][1] = -1.0;      // No Traces in this Azimuth bin. Aaaww..
		}
		
		/***********  LOG SUMMARY ************/
		logReport.append( "Bin Size: " ); logReport.append( to_string(minBinSze) );
		logReport.append( "Epicentre: " );
		logReport.append( to_string(episcan) );
		logReport.append( " , No of RFs in bin:" );
		logReport.append( to_string( binCnt ) );
		logReport.append( " \n");
		/************ END SUMMARY ************/
		
		epicStack[iterStck][0] = episcan;
		episcan += bininc;
		iterStck++;

	}
	/***********************  END Epic Stack ******************************/
	
	/* Iterate through Stack and use MTCDriver to invert the stack. 
	 I put this outside Stack construction because decisions need to be 
	 made about whether to print an azimuth bin depending on its HIT status	  
	 */
	int npre = preTime / DeltaT, npost = postTime / DeltaT;
	outLen = npre + npost;                 // Variable holds output legth
	radialRFAscii.assign(stackSze, outLen, 0.0);
	transRFAscii.assign(stackSze, outLen, 0.0);
	timeAscii.assign(stackSze, outLen, 0.0);
	
	for (int iterStck = 0; iterStck < stackSze; iterStck++) {
		
		for (int itercmp = 0; itercmp < 2; itercmp++) {
			
			VecComplex RFfreq(stackRFTrace.dim3(), stackRFTrace[iterStck][itercmp]);
			VecDoub_IO RFtimebyAzim;
			MTCDriver invertStack;
			
			invertStack.icmplxfft(RFfreq, RFtimebyAzim, padNyq);
			
			// Dump RF here.. big bug. Yish!!
			/*for (int i = 0; i < RFtimebyAzim.size(); i++) {
			 debugFile << abs(RFtimebyAzim[i]) << "  " << endl;
			 } */
			//Build Reciever Function in time domain? 
			//Causative & Non Causative Build
			if (itercmp == 0) {
				buildAsciiOut(preTime, postTime, RFtimebyAzim, radialRFAscii, timeAscii,
							  iterStck);
			}else {
				buildAsciiOut(preTime, postTime, RFtimebyAzim, transRFAscii, timeAscii,
							  iterStck);
			}
			
		}
	}
	
	/***************************** Done inverting and building *******************************/
	logReport.append("Total no of RFs in wedge:  ");
	logReport.append( to_string( wdgeCnt ) );
	logReport.append( "\n" );
	
	
}


/* @OLUGBOJI2016: 
 Do jacknife stack by re-estimating RF stacks nRB (number of records in bin bucket) times
 for every nRB times we drop 1 record in the bin sample using only nRB - 1 records
 This results in a jacknife sample from which 3 new traces can be reconstructed:
     Trace 1:  average
	 Trace 2 & 3: Trace1 +- stdev RF
 Unlike stackEpic above, I need a temporary store of record indices in bin bucket
 Fills stackRFTrace, stacAzim & calls the inverseFFT, okay consult
 Instead of Azimuth on vertical axis, we have Epic Distance ...
 View this as a wedge cut out from the a section close to the station
 */
void TraceStack::jacknifeStackEpic(float epbazmin, float epbazmax, float min, float max, float inc,
						   Mat3DCmplx &RF, Mat3DDoub &delta) {
	

	/* Put sections here into my TraceStack class*/
	
	/* Initialize your 3Dimensional arrays. Extension:
	 Write data structure driver in nr3.h
	 Here I allocate the same memory structure ..
	 Solution Point: Tentative problem: I can't initialize NRMat3D look into
	 this. For now put this initialization in stackAzim routine*/
	Mat3DCmplx RFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	Mat3DDoub sigmaRFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	
	
	/* Deep Copy Here */
	for (int first = 0; first < RF.dim1(); first++) {
		for (int third = 0; third < RF.dim3(); third++) {
			for (int second = 0; second < RF.dim2(); second++) {
				RFTrace[first][second][third] = RF[first][second][third];
				sigmaRFTrace[first][second][third] = pow(delta[first][second][third],2.0);
			}
		}
	}
	
	/* Parameters to track if RFs are in wedge and stack 'em */
	float dist2Bin;  // Distance between current distance and bin edge
	float episcan;
	float epi2bin;
	int minBinSze = binSze;
	
	binmin = min;
	binmax = max;
	bininc = inc;
	
	episcan = binmax;
	if (bininc > 0.0) bininc = -1.0 * bininc;
	
	//Determine the stack size using bin parameters
	int stackSze = 0;
	while (episcan >= binmin) {
		stackSze++;
		episcan += bininc;
	}
	
	// *********************************************************** stats report header
	logReport.append("#Epicentral Stack Called. \n#Wedge Azimuth,");
	logReport.append("Min:");logReport.append(to_string(epbazmin));
	logReport.append("  Max:");logReport.append(to_string(epbazmax)); logReport.append("\n");
	
	logReport.append("#Epicentral distance bounds(deg),Min:"); logReport.append(to_string(min));logReport.append("  Max:"); logReport.append(to_string(max)); logReport.append( "  Bin Width:" ); logReport.append( to_string(bininc) ); logReport.append("\n");
	
	logReport.append( "#Bin centre (deg) , No. of RFs in bin \n " );
	
	// ************************************************************** end stats header
	
	// use information above to initialize bin statistic variables used in output
	setPrintStats = true;
	BinCntStat.assign(stackSze, 0);
	BinCntrStat.assign(stackSze, 0);
	
	// reset the scan variable
	episcan = binmax;
	
	/***********************  Initialize time domain RF traces (sample, mean and stdev from jacknife)  ******************************/
	/*
	 Iterate through Stack and use MTCDriver to invert the stack.
	 I put this outside Stack construction because decisions need to be
	 made about whether to print an azimuth bin depending on its HIT status
	 */
	int npre = preTime / DeltaT, npost = postTime / DeltaT;
	outLen = npre + npost;                 // Variable holds output legth
	
	// record index in particular bin. should be same size as bin count
	recIndxInBin.assign(Azim.size(), 0);
	
	// add extra variance biases from the jacknife estimation ...
	radialRFAscii.assign(stackSze, outLen, 0.0);
	avgRadialRFAscii.assign(stackSze, outLen, 0.0);
	devRadialRFAscii.assign(stackSze, outLen, 0.0);
	
	transRFAscii.assign(stackSze, outLen, 0.0);
	avgTransRFAscii.assign(stackSze, outLen, 0.0);
	devTransRFAscii.assign(stackSze, outLen, 0.0);
	
	timeAscii.assign(stackSze, outLen, 0.0);
	
	/***********************  Initialize 3D Stack Matrix, note that some of the bins might not contain any traces  ************************/
	/*
	 So I use another vector 'epicStack' to check if a bin within the stack contains traces:
	 To do this I initialize and fill the azimStack matrix:
	 Column1: Left Epic boundary in Bin
	 Column2: status of 'bin sum', -1.0 if a bin is not filled with traces; 1.0 f it is
	 */
	Mat3DCmplx stackRFTrace(stackSze, RFTrace.dim2(), nPad); //[dim1: bin][dim2: rad,trans][dim3:freq]
	for (int first = 0; first < stackSze; first++) {
		for (int third = 0; third < nPad; third++) {
			for (int second = 0; second < 2; second++) {
				stackRFTrace[first][second][third] = zero;
			}
		}
	}
	
	/***********************  Initialize bookeeping arrray here ... revisit nr3.h ************************/
	// Initialize the Array here ... revisit nr3.h
	epicStack.assign(stackSze,2,-1.0); // Read comment above on why this is initialized to -1.0
	
	int iterStck = 0;
	int wdgeCnt = 0; // No of Traces that fall into sector wedge
	bool parseOnce = false;
	
	// Scan variable scans down, that's why I have it initialized to outer boundary of the wedge
	// This scan is currently inefficient. Revisit. Think of ways to make more efficient.
	while (episcan >= binmin) {
		int binCnt = 0; // No of Traces in bin .. Min. no for a BinHit is determined by minBinSze
		

		
		// Run through all Records ...
		for (int irec = 0; irec < Azim.size() ; irec++) {
			bool isInWedge = (Azim[irec]>= epbazmin && Azim[irec] <= epbazmax);
			dist2Bin = abs( Epic[irec] - episcan );
			bool isInBin = dist2Bin <= abs(bininc);
			
			if( isInWedge && !parseOnce){
				wdgeCnt++;
			}
			
			// Verify presence in wedge and in bin and then stack! dependent on status of flasgs. Review or debug
			// @ OLUGBOJI2016: Since stack is done with jacknife. Need to build size and record index in bin bucket first
			if (isInWedge && isInBin) {
				
				// save record index for use by jacknife loop
				recIndxInBin[binCnt] = irec;
				binCnt++;
				
				// Stack RFs using their coherence estimates
				/*
				for (int itercmp = 0; itercmp < 2; itercmp++) {
					for (int iternf = 0; iternf < nFreqMax; iternf++) {
						StackRFHolder[itercmp][iternf] +=
						RFTrace[irec][itercmp][iternf] / sigmaRFTrace[irec][itercmp][iternf];
						StackSigmaHolder[itercmp][iternf] += 1.0 / sigmaRFTrace[irec][itercmp][iternf];
					}
				}
				 */
				
			}
		}
		/***********************  END Bin Scan ******************************/
		
		
		parseOnce = true;   // use this to see if to re-run jacknife binCnt times?
		
		// update statistics for future use ...
		BinCntStat[iterStck] = binCnt;
		BinCntrStat[iterStck] = episcan;
		
		// Push in new Stack Trace, depending on Scan ...
		
		
		if (binCnt >= minBinSze) {
			cout << "Demo" << endl;
			
			/***********************  Run Jacknife for the remaning binCnt-1 traces here ... ******************************/
			iJckRadialRFAscii.assign(binCnt, outLen, 0.0);
			iJckTransRFAscii.assign(binCnt, outLen, 0.0);
			int nxtJck = 0;
			
			for (int nxtJck = 0; nxtJck < binCnt; nxtJck++) {
				/* The two variables below hold the stack within this bin.
				 Obviously they'll be empty if there is no RF in this Epicentral scan */
				MatComplex StackRFHolder(RFTrace.dim2(), nPad, zero);
				MatDoub StackSigmaHolder(RFTrace.dim2(), nPad, iZero);
				
				// sum sample .
				for (int iJck = 0; iJck < binCnt; iJck++) {
					if (iJck != nxtJck) {
						// use this records in the sum
						int irec = recIndxInBin[iJck];
						cout << "rec " << irec << endl;
						
						for (int itercmp = 0; itercmp < 2; itercmp++) {
							for (int iternf = 0; iternf < nFreqMax; iternf++) {
								StackRFHolder[itercmp][iternf] +=
								RFTrace[irec][itercmp][iternf] / sigmaRFTrace[irec][itercmp][iternf];
								StackSigmaHolder[itercmp][iternf] += 1.0 / sigmaRFTrace[irec][itercmp][iternf];
							}
						}
					} else {
						///ignore the record in the sum
					}
				}
				
				// Stack RFs using their coherence estimates: used all records from first pass above - 1st Jack
				for (int itercmp = 0; itercmp < 2; itercmp++) {
					for (int iternf = 0; iternf < nFreqMax; iternf++) {
						
						Hann = pow( cos( Pih * double(iternf/nFreqMax) ), 2.0);
						stackRFTrace[iterStck][itercmp][iternf] =
						Hann * StackRFHolder[itercmp][iternf] / StackSigmaHolder[itercmp][iternf];
						
					}
				}
				
				
				// Insert here the time domain code and then update the RF traces in time domain ...
				// Time domain after 1st Jack
				for (int itercmp = 0; itercmp < 2; itercmp++) {
					
					VecComplex RFfreq(stackRFTrace.dim3(), stackRFTrace[iterStck][itercmp]);
					VecDoub_IO RFtimebyAzim;
					MTCDriver invertStack;
					
					invertStack.icmplxfft(RFfreq, RFtimebyAzim, padNyq);
					
					// Dump RF here.. big bug. Yish!!
					/*for (int i = 0; i < RFtimebyAzim.size(); i++) {
					 debugFile << abs(RFtimebyAzim[i]) << "  " << endl;
					 } */
					//Build Reciever Function in time domain?
					//Causative & Non Causative Build
					if (itercmp == 0) {
						buildAsciiOut(preTime, postTime, RFtimebyAzim, radialRFAscii, timeAscii,
									  iterStck);
					}else {
						buildAsciiOut(preTime, postTime, RFtimebyAzim, transRFAscii, timeAscii,
									  iterStck);
					}
					
				}
				
				// Update the average in frequency or time domain (In time domain!) ??
				for (int itime = 0; itime < outLen; itime++) {
					
					iJckRadialRFAscii[nxtJck][itime] = radialRFAscii[iterStck][itime];
					iJckTransRFAscii[nxtJck][itime] = transRFAscii[iterStck][itime];
					
					avgRadialRFAscii[iterStck][itime] = avgRadialRFAscii[iterStck][itime] + radialRFAscii[iterStck][itime];
					avgTransRFAscii[iterStck][itime] = avgTransRFAscii[iterStck][itime] + transRFAscii[iterStck][itime];
				}
    
			}
			

			
			cout << "Jack all complete " << binCnt << endl;
			//cout << "Jack done once: " << endl;
			// Now rebuild stack before recomputing - re-initialize data stack matrices
			
			
			/*
			while (nxtJck < binCnt-1){
				cout << "Jack " << nxtJck << endl;
				MatComplex StackRFHolder(RFTrace.dim2(), nPad, zero);
				MatDoub StackSigmaHolder(RFTrace.dim2(), nPad, iZero);
				
				//cout << "Next Jack" << endl;
				
				// Re-stack RFs using their coherence estimates, dropping one record every time
				for (int binIndx = 0; binIndx < (binCnt - nxtJck - 1) ; binIndx++) {
					cout  << "binindx " << binIndx << endl;
					int irec = recIndxInBin[binIndx];
					cout << "rec " << irec << endl;
					
					for (int itercmp = 0; itercmp < 2; itercmp++) {
						for (int iternf = 0; iternf < nFreqMax; iternf++) {
							StackRFHolder[itercmp][iternf] +=
							RFTrace[irec][itercmp][iternf] / sigmaRFTrace[irec][itercmp][iternf];
							StackSigmaHolder[itercmp][iternf] += 1.0 / sigmaRFTrace[irec][itercmp][iternf];
						}
					}
					
				}
				
				for (int itercmp = 0; itercmp < 2; itercmp++) {
					for (int iternf = 0; iternf < nFreqMax; iternf++) {
						
						Hann = pow( cos( Pih * double(iternf/nFreqMax) ), 2.0);
						stackRFTrace[iterStck][itercmp][iternf] =
						Hann * StackRFHolder[itercmp][iternf] / StackSigmaHolder[itercmp][iternf];
						
					}
				}
				
				
				// Insert here the time domain code and then update the RF traces in time domain ...
				// Time domain after each binCnt -1 Jack
				for (int itercmp = 0; itercmp < 2; itercmp++) {
					
					VecComplex RFfreq(stackRFTrace.dim3(), stackRFTrace[iterStck][itercmp]);
					VecDoub_IO RFtimebyAzim;
					MTCDriver invertStack;
					
					invertStack.icmplxfft(RFfreq, RFtimebyAzim, padNyq);
					
					// Dump RF here.. big bug. Yish!!
					//for (int i = 0; i < RFtimebyAzim.size(); i++) {
					// debugFile << abs(RFtimebyAzim[i]) << "  " << endl;
					// }
					//Build Reciever Function in time domain?
					//Causative & Non Causative Build
					if (itercmp == 0) {
						buildAsciiOut(preTime, postTime, RFtimebyAzim, radialRFAscii, timeAscii,
									  iterStck);
					}else {
						buildAsciiOut(preTime, postTime, RFtimebyAzim, transRFAscii, timeAscii,
									  iterStck);
					}
					
				}
				
				cout << "before update" << endl;
				// Update the average in frequency or time domain (In time domain!) next jack Estimate
				nxtJck = nxtJck + 1;
				for (int itime = 0; itime < outLen; itime++) {
					iJckRadialRFAscii[nxtJck][itime] = radialRFAscii[iterStck][itime];
					iJckTransRFAscii[nxtJck][itime] = transRFAscii[iterStck][itime];
					
					avgRadialRFAscii[iterStck][itime] = avgRadialRFAscii[iterStck][itime] + radialRFAscii[iterStck][itime];
					avgTransRFAscii[iterStck][itime] = avgTransRFAscii[iterStck][itime] + transRFAscii[iterStck][itime];
				}
				
			}
			*/
			/***********************  End Jacknife Sum here ... ******************************/
			
			cout << "aggregating" << endl;
			
			// scale by binCnt before variance computation ....
			for (int itime = 0; itime < outLen; itime++) {
				avgRadialRFAscii[iterStck][itime] = avgRadialRFAscii[iterStck][itime] /binCnt;
				avgTransRFAscii[iterStck][itime] = avgTransRFAscii[iterStck][itime] /binCnt;
			}
			
			// build jacknife variance after sum above ... use ith jack stored.
			for (int nxtJck = 0; nxtJck < binCnt; nxtJck++) {
				for (int itime = 0; itime < outLen; itime++) {
					
					devRadialRFAscii[iterStck][itime] +=
					         pow((iJckRadialRFAscii[nxtJck][itime] - avgRadialRFAscii[iterStck][itime]), 2.0);
					
					devTransRFAscii[iterStck][itime] +=
					        pow( (iJckTransRFAscii[nxtJck][itime] - avgTransRFAscii[iterStck][itime]), 2.0);
				}
			}
			
			// scale variance by jacknife normalization parameters ...
			for (int itime = 0; itime < outLen; itime++) {
				devRadialRFAscii[iterStck][itime] =  sqrt( (binCnt - 1) * devRadialRFAscii[iterStck][itime]/binCnt );
				devTransRFAscii[iterStck][itime] = sqrt( (binCnt - 1) * devTransRFAscii[iterStck][itime]/binCnt );
			}
			
			/***********************  End Jacknife Variance here ... ******************************/
			
			epicStack[iterStck][1] = 1.0;      // Particular Azimuth bin hit. Yay!
		}else {
			// Should I not surpress the entry in stak here, if I can't get any
			// Reciever functions within this azimuth bin?
			epicStack[iterStck][1] = -1.0;      // No Traces in this Azimuth bin. Aaaww..
		}
		
		
		/***********  LOG SUMMARY ************/
		//logReport.append( "Epicentre: " );
		logReport.append( to_string(episcan) );
		logReport.append( " , " );
		logReport.append( to_string( binCnt ) );
		logReport.append( " \n");
		cout << "binCnt " << binCnt << endl;
		/************ END SUMMARY ************/
		
		epicStack[iterStck][0] = episcan;
		episcan += bininc;
		iterStck++;
		
	}
	/***********************  END Epic Stack ******************************/
	
	
	/***************************** Done inverting and building *******************************/
	/***************************** Done this nRBin times to estimate the jacknife ?? *******************************/
	logReport.append("#Total no of RFs in wedge:  ");
	logReport.append( to_string( wdgeCnt ) );
	logReport.append( "\n" );
	
	
	
	
	
}


/* Hear I do a summary stack. Collapse only the Radial RFs into A single Number in 
   [zero +- delta]time
 */
double TraceStack::grdStack(float epbazmin, float epbazmax, float min, float max, float inc,
						   Mat3DCmplx &RF, Mat3DDoub &delta) {
	
	logReport.append(" Grid Stack Called With Sector Parameters ... : \n ");
	logReport.append("Baz min:  ");logReport.append(to_string(epbazmin)); logReport.append("\n  ");
	logReport.append("Baz max:  ");logReport.append(to_string(epbazmax)); logReport.append("\n  ");
	logReport.append("Epic min: "); logReport.append(to_string(min)); logReport.append("\n  ");
	logReport.append("Epic max:  "); logReport.append(to_string(max)); logReport.append("\n  ");
	
	/* Put sections here into my TraceStack class*/
	
	/* Initialize your 3Dimensional arrays. Extension: 
	 Write data structure driver in nr3.h
	 Here I allocate the same memory structure ..
	 Solution Point: Tentative problem: I can't initialize NRMat3D look into
	 this. For now put this initialization in stackAzim routine*/
	Mat3DCmplx RFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	Mat3DDoub sigmaRFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	
	
	/* Deep Copy Here */
	for (int first = 0; first < RF.dim1(); first++) {
		for (int third = 0; third < RF.dim3(); third++) {
			for (int second = 0; second < RF.dim2(); second++) {
				RFTrace[first][second][third] = RF[first][second][third];
				sigmaRFTrace[first][second][third] = pow(delta[first][second][third],2.0);
			}
		}
	}
	
	/* Parameters to track if RFs are in wedge and stack 'em */
	float dist2Bin;  // Distance between current distance and bin edge
	float episcan;
	float epi2bin;
	int minBinSze = binSze;
	
	binmin = min;
	binmax = max;
	bininc = inc;
	
	episcan = binmax;
	if (bininc > 0.0) bininc = -1.0 * bininc;
	
	//Determine the stack size using bin parameters
	int stackSze = 0;
	while (episcan >= binmin) {
		stackSze++;
		episcan += bininc;
	}
	// reset the scan variable 
	episcan = binmax;
	
	/* Initialize 3D Stack Matrix, note that some of the bins might not contain any traces,
	 So I use another vector 'epicStack' to check if a bin within the stack contains traces:
	 To do this I initialize and fill the azimStack matrix:
	 Column1: Left Epic boundary in Bin
	 Column2: status of 'bin sum', -1.0 if a bin is not filled with traces; 1.0 f it is
	 */
	Mat3DCmplx stackRFTrace(stackSze, RFTrace.dim2(), nPad); //[dim1: bin][dim2: rad,trans][dim3:freq]
	for (int first = 0; first < stackSze; first++) {
		for (int third = 0; third < nPad; third++) {
			for (int second = 0; second < 2; second++) {
				stackRFTrace[first][second][third] = zero;
			}
		}
	}	  // Initialize the Array here ... revisit nr3.h
	epicStack.assign(stackSze,2,-1.0); // Read comment above on why this is initialized to -1.0
	
	int iterStck = 0;
	int wdgeCnt = 0; // No of Traces that fall into sector wedge
	bool parseOnce = false;
	
	// Scan variable scans down, that's why I have it initialized to outer boundary of the wedge
	// This scan is currently inefficient. Revisit. Think of ways to make more efficient.
	while (episcan >= binmin) { 
		int binCnt = 0; // No of Traces in bin .. Min. no for a BinHit is determined by minBinSze
		
		/* The two variables below hold the stack within this bin. 
		 Obviously they'll be empty if there is no RF in this Epicentral scan */
		MatComplex StackRFHolder(RFTrace.dim2(), nPad, zero);  
		MatDoub StackSigmaHolder(RFTrace.dim2(), nPad, iZero);
		
		// Run through all Records ...
		for (int irec = 0; irec < Azim.size() ; irec++) {
			bool isInWedge = (Azim[irec]>= epbazmin && Azim[irec] <= epbazmax);
			dist2Bin = abs( Epic[irec] - episcan );
			bool isInBin = dist2Bin <= abs(bininc);
			
			if( isInWedge && !parseOnce){
				wdgeCnt++;
			} 
			
			// Verify presence in wedge and in bin and then stack! dependent on status of flasgs. Review or debug
			if (isInWedge && isInBin) {
				binCnt++;
				// Stack RFs using their coherence estimates 
				for (int itercmp = 0; itercmp < 2; itercmp++) {
					for (int iternf = 0; iternf < nFreqMax; iternf++) { 
						StackRFHolder[itercmp][iternf] += 
						RFTrace[irec][itercmp][iternf] / sigmaRFTrace[irec][itercmp][iternf];
						StackSigmaHolder[itercmp][iternf] += 1.0 / sigmaRFTrace[irec][itercmp][iternf];
					}
				}
				
			}
		}
		/***********************  END Bin Scan ******************************/
		parseOnce = true;
		// Push in new Stack Trace, depending on Scan ...
		if (binCnt >= minBinSze) {
			for (int itercmp = 0; itercmp < 2; itercmp++) {
				for (int iternf = 0; iternf < nFreqMax; iternf++) {
					
					Hann = pow( cos( Pih * double(iternf/nFreqMax) ), 2.0);
					stackRFTrace[iterStck][itercmp][iternf] = 
					Hann * StackRFHolder[itercmp][iternf] / StackSigmaHolder[itercmp][iternf];
					
				}
			}
			epicStack[iterStck][1] = 1.0;      // Particular Azimuth bin hit. Yay!
		}else {
			// Should I not surpress the entry in stak here, if I can't get any
			// Reciever functions within this azimuth bin?
			epicStack[iterStck][1] = -1.0;      // No Traces in this Azimuth bin. Aaaww..
		}
		
		/***********  LOG SUMMARY ************/
		logReport.append( "StackSize: " ); logReport.append( to_string(stackSze) );
		logReport.append( "Epicentre: " );
		logReport.append( to_string(episcan) );
		logReport.append( " , No of RFs in bin:" );
		logReport.append( to_string( binCnt ) );
		logReport.append( " \n");
		/************ END SUMMARY ************/
		
		epicStack[iterStck][0] = episcan;
		episcan += bininc;
		iterStck++;
		
	}
	/***********************  END Epic Stack ******************************/
	
	/* Iterate through Stack and use MTCDriver to invert the stack. 
	 I put this outside Stack construction because decisions need to be 
	 made about whether to print an azimuth bin depending on its HIT status	 
	 !!!!!!!!!!!!!!!!!!!!!! UUUUSEEEE ONLY Radial Stacks, and pick only zero time!!
	 */
	int npre = preTime / DeltaT, npost = postTime / DeltaT;
	outLen = npre + npost;                 // Variable holds output legth
	radialRFAscii.assign(stackSze, outLen, 0.0);
	transRFAscii.assign(stackSze, outLen, 0.0);
	timeAscii.assign(stackSze, outLen, 0.0);
	
	
	for (int iterStck = 0; iterStck < stackSze; iterStck++) {
		
		
		
		//!!! Hardwired for only radial stack since isotropic signal needed
		int maxIterCmp = 1;
		for (int itercmp = 0; itercmp < maxIterCmp; itercmp++) {
			
			VecComplex RFfreq(stackRFTrace.dim3(), stackRFTrace[iterStck][itercmp]);
			VecDoub_IO RFtimebyAzim;
			MTCDriver invertStack;
			
			invertStack.icmplxfft(RFfreq, RFtimebyAzim, padNyq);
			
			// Dump RF here.. big bug. Yish!!
			/*for (int i = 0; i < RFtimebyAzim.size(); i++) {
			 debugFile << abs(RFtimebyAzim[i]) << "  " << endl;
			 } */
			//Build Reciever Function in time domain? 
			//Causative & Non Causative Build
			preTime = 2; postTime = 2;
			if (itercmp == 0) {
				buildAsciiOut(preTime, postTime, RFtimebyAzim, radialRFAscii, timeAscii,
							  iterStck);
			}else {
				buildAsciiOut(preTime, postTime, RFtimebyAzim, transRFAscii, timeAscii,
							  iterStck);
			}
			
		}
	}
	
	/***** Find Maximum Radial and Timing?
	 */
	double maxRad, maxTime;
	
	findMAX(radialRFAscii, timeAscii, maxRad, maxTime);
	 
	
	/***************************** Done inverting and building *******************************/
	logReport.append("Total no of RFs in wedge:  ");
	logReport.append( to_string( wdgeCnt ) );
	logReport.append( "\n" );
	return maxRad;
	
	
}


/* Utility method for printing the ascii traces into a file .. 
   New Adaptation: I need to convert the frequency decimation to time decimation
 */
void TraceStack::print(char* outfn, int fileFlag, int timeShift){

	bool outdata = true; // Check here for file, or outside the fucctiion?
						// mute if already printed in case check - jacknife ...
	
	
	string fileRadial(outfn);
	string  fileTrans(outfn);
	
	// Open new file if depth migration routine activated ...
	string fileRadialMig(outfn);
	string fileTransMig(outfn);
	
	//4 files for jacknife uncertainties
	string fileRadialDevMin(outfn);
	string fileRadialDevMax(outfn);
	string fileTransDevMin(outfn);
	string fileTransDevMax(outfn);
	
	// File to save bin statistics, use in conjuction with uncertainty to measure data quality. use only if setPrintStats is true
	string filePrintStats(outfn);
	
	switch (fileFlag) {
		case STACKAZIM:
			fileRadial.append(".Azim_Rad.xyz");
			fileTrans.append(".Azim_Trans.xyz");
			
			fileRadialMig.append("_Depth.Azim_Rad.xyz");
			fileTransMig.append("_Depth.Azim_Trans.xyz");
			
			break;
		case STACKEPIC:
			fileRadial.append(".Epic_Rad.xyz");
			fileTrans.append(".Epic_Trans.xyz");
			
			fileRadialMig.append("_Depth.Epic_Rad.xyz");
			fileTransMig.append("_Depth.Epic_Trans.xyz");
			
			break;
		case STACKGRD:
			fileRadial.append(".GrdStckRad.xyz");
			fileTrans.append(".GrdStckTrans.xyz");
		
			break;
		case STACKEPICJCK:
			// 4 extra files. this is why jacknife is a pain in the rear...
			fileRadial.append(".Epic_Rad_Jack.xyz");
			fileTrans.append(".Epic_Trans_Jack.xyz");
			
			//ofstream radialout(fileRadial.c_str());
			//ofstream transout(fileTrans.c_str());

			fileRadialDevMin.append(".Epic_Rad_Jack.devMin.xyz");
			fileRadialDevMax.append(".Epic_Rad_Jack.devMax.xyz");
			
			//ofstream radialOutMin(fileRadialDevMin.c_str());
			//ofstream radialOutMax(fileRadialDevMax.c_str());
			
			fileTransDevMin.append(".Epic_Trans_Jack.devMin.xyz");
			fileTransDevMax.append(".Epic_Trans_Jack.devMax.xyz");
			
			//ofstream transOutMin(fileTransDevMin.c_str());
			//ofstream transOutMax(fileTransDevMax.c_str());
		
			break;
			
		default:
			break;
	}
	
	ofstream radialout(fileRadial.c_str());
	ofstream transout(fileTrans.c_str());
	
	ofstream radialMigOut(fileRadialMig.c_str());
	ofstream transMigOut(fileTransMig.c_str());
	
	ofstream radialOutMin(fileRadialDevMin.c_str());
	ofstream radialOutMax(fileRadialDevMax.c_str());
	
	ofstream transOutMin(fileTransDevMin.c_str());
	ofstream transOutMax(fileTransDevMax.c_str());
	
	cout << "Opened file: " << outfn << endl;
	
	// Write files using templates from Epicentral Dump below ...
	int stackSze = radialRFAscii.nrows();   // Hopefully, this holds row size
	for (int iterStck = 0; iterStck < stackSze; iterStck++) {
		// If there is no RF in bin, skip printing!! Reverse check here
		// For fun. Why? so you can see what inversion does.
		if (epicStack[iterStck][1] < 0) {
			
		}else {
			cout << "Epicentral bin: " << epicStack[iterStck][0] << endl;
			// Dump in File if Azimuth bin indicates hit
			for (int itercmp = 0; itercmp < 2; itercmp++) {
				
				if (itercmp == 0) {
					
					if (flagDepthMig) {
						
						// Normal Output With Time - RFamp - Distance Triplet
						for (int i = 0; i < (outLen); i++) {
							radialout << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgRadialRFAscii[iterStck][i]
							<< endl;
						}
						
						// plot uncertainties from jacknife analysis ...
						for (int i = 0; i < (outLen); i++) {
							radialOutMin << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgRadialRFAscii[iterStck][i] - devRadialRFAscii[iterStck][i]
							<< endl;
							
							
							radialOutMax << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgRadialRFAscii[iterStck][i] + devRadialRFAscii[iterStck][i]
							<< endl;
						}
						
						radialout  << ">" << endl;
						radialOutMin  << ">" << endl;
						radialOutMax  << ">" << endl;
						
					}else{
						
						// Normal Output With Time - RFamp - Distance Triplet
						for (int i = 0; i < (outLen); i++) {
							radialout << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgRadialRFAscii[iterStck][i]
							<< endl;
						}
						
						// plot uncertainties from jacknife analysis ...
						for (int i = 0; i < (outLen); i++) {
							radialOutMin << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgRadialRFAscii[iterStck][i] - devRadialRFAscii[iterStck][i]
							<< endl;
							
							
							radialOutMax << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgRadialRFAscii[iterStck][i] + devRadialRFAscii[iterStck][i]
							<< endl;
						}
						
						radialout  << ">" << endl;
						radialOutMin  << ">" << endl;
						radialOutMax  << ">" << endl;
						
					}
					
					
				}else {
					// I should either change name or find a way to make this
					// manageable ...
					
					if (flagDepthMig) {
						
						// Output File in Time Domain for Second File ****************************** 2.0 ....
						// Normal Output With Time - RFamp - Distance Triplet
						for (int i = 0; i < (outLen); i++) {
							transout << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgTransRFAscii[iterStck][i]
							<< endl;
						}
						
						// plot uncertainties from jacknife analysis ...
						for (int i = 0; i < (outLen); i++) {
							transOutMin << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgTransRFAscii[iterStck][i] - devTransRFAscii[iterStck][i]
							<< endl;
							
							
							transOutMax << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgTransRFAscii[iterStck][i] + devTransRFAscii[iterStck][i]
							<< endl;
						}
						
						transout  << ">" << endl;
						transOutMin  << ">" << endl;
						transOutMax  << ">" << endl;
						
						
						// Output File in Depth Domain for Second File ****************************** 2.0 ....
						for (int i = 0; i < (outLen); i++) {
							transMigOut << depthAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< radialRFAscii[iterStck][i]
							<< endl;
						}
						transMigOut  << ">" << endl;
						
					}else{
						
						// Normal Output With Time - RFamp - Distance Triplet
						for (int i = 0; i < (outLen); i++) {
							transout << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgTransRFAscii[iterStck][i]
							<< endl;
						}
						
						// plot uncertainties from jacknife analysis ...
						for (int i = 0; i < (outLen); i++) {
							transOutMin << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgTransRFAscii[iterStck][i] - devTransRFAscii[iterStck][i]
							<< endl;
							
							
							radialOutMax << timeShift + timeAscii[iterStck][i] << "\t"
							<< epicStack[iterStck][0] << "\t"
							<< avgTransRFAscii[iterStck][i] + devTransRFAscii[iterStck][i]
							<< endl;
						}
						
						transout  << ">" << endl;
						transOutMin  << ">" << endl;
						transOutMax  << ">" << endl;
						
					}
					
					
				}
				
			}
		}
		
		
	}
	outdata = false;  // mute after print once  ...

	if (outdata) {
		// Go through the stack and dump in file, first dimension
		int stackSze = radialRFAscii.nrows();   // Hopefully, this holds row size
		
		//Here is the dump code snippet . I use here again the file flag - now more than 2 options ? how handle others?
		if (fileFlag == 0) {
			// Azimuth Dump here ...
			for (int iterStck = 0; iterStck < stackSze; iterStck++) {
				// If there is no RF in bin, skip printing!! Reverse check here
				// For fun. Why? so you can see what inversion does.
				if (azimStack[iterStck][1] < 0) {
					
				}else {
					cout << "Azimuth bin " << azimStack[iterStck][0] << endl;
					// Dump in File if Azimuth bin indicates hit
					for (int itercmp = 0; itercmp < 2; itercmp++) {
						
						if (itercmp == 0) {
							
							
							if (flagDepthMig) {
								
								// Output File in Time Domain for First File ****************************** 1.0 ....
								for (int i = 0; i < (outLen); i++) {
									radialout << timeShift + timeAscii[iterStck][i] << "\t"
									<< azimStack[iterStck][0] << "\t"
									<< radialRFAscii[iterStck][i]
									<< endl;
								}
								radialout  << ">" << endl;
								
								// Output File in Depth Domain for First File ****************************** 2.0 ....
								for (int i = 0; i < (outLen); i++) {
									radialMigOut << depthAscii[iterStck][i] << "\t"
									<< azimStack[iterStck][0] << "\t"
									<< radialRFAscii[iterStck][i]
									<< endl;
								}
								radialMigOut  << ">" << endl;
								
							}else{
								
								// Normal Output With Time - RFamp - Distance Triplet
								for (int i = 0; i < (outLen); i++) {
									radialout << timeShift + timeAscii[iterStck][i] << "\t"
									<< azimStack[iterStck][0] << "\t"
									<< radialRFAscii[iterStck][i]
									<< endl;
								}
								radialout  << ">" << endl;
								
							}
							
							
						}else {
							// I should either change name or find a way to make this
							// manageable ...
							
							
							if (flagDepthMig) {
								
								// Output File in Time Domain for First File ****************************** 1.0 ....
								for (int i = 0; i < (outLen); i++) {
									transout  << timeShift + timeAscii[iterStck][i] << "\t"
									<< azimStack[iterStck][0] << "\t"
									<< transRFAscii[iterStck][i]
									<< endl;
								}
								transout  << ">" << endl;
								
								// Output File in Depth Domain for First File ****************************** 2.0 ....
								for (int i = 0; i < (outLen); i++) {
									transMigOut << depthAscii[iterStck][i] << "\t"
									<< azimStack[iterStck][0] << "\t"
									<< transRFAscii[iterStck][i]
									<< endl;
								}
								transMigOut  << ">" << endl;
								
							}else{
								
								// Normal Output With Time - RFamp - Distance Triplet
								for (int i = 0; i < (outLen); i++) {
									transout << timeShift + timeAscii[iterStck][i] << "\t"
									<< azimStack[iterStck][0] << "\t"
									<< transRFAscii[iterStck][i]
									<< endl;
								}
								transout  << ">" << endl;
								
							}
							
							
						}
						
					}
				}
				
				
			}
			/************************** END AZIM. DUMP HERE *******************************/
		} else {
			// Epicentral Dump here
			for (int iterStck = 0; iterStck < stackSze; iterStck++) {
				// If there is no RF in bin, skip printing!! Reverse check here
				// For fun. Why? so you can see what inversion does.
				if (epicStack[iterStck][1] < 0) {
					
				}else {
					cout << "Epicentral bin: " << epicStack[iterStck][0] << endl;
					// Dump in File if Azimuth bin indicates hit
					for (int itercmp = 0; itercmp < 2; itercmp++) {
						
						if (itercmp == 0) {
							
							if (flagDepthMig) {
								
								// Output File in Time Domain for First File ****************************** 1.0 ....
								for (int i = 0; i < (outLen); i++) {
									radialout << timeShift + timeAscii[iterStck][i] << "\t"
									<< epicStack[iterStck][0] << "\t"
									<< radialRFAscii[iterStck][i]
									<< endl;
								}
								radialout  << ">" << endl;
								
								
								// Output File in Depth Domain for Second File ****************************** 2.0 ....
								
								for (int i = 0; i < (outLen); i++) {
									radialMigOut << depthAscii[iterStck][i] << "\t"
									<< epicStack[iterStck][0] << "\t"
									<< radialRFAscii[iterStck][i]
									<< endl;
								}
								radialMigOut  << ">" << endl;
								
								
							}else{
								
								// Normal Output With Time - RFamp - Distance Triplet
								for (int i = 0; i < (outLen); i++) {
									radialout << timeShift + timeAscii[iterStck][i] << "\t"
									<< epicStack[iterStck][0] << "\t"
									<< radialRFAscii[iterStck][i]
									<< endl;
								}
								radialout  << ">" << endl;
								
							}
							
							
						}else {
							// I should either change name or find a way to make this
							// manageable ...
							
							if (flagDepthMig) {
								
								// Output File in Time Domain for Second File ****************************** 2.0 ....
								for (int i = 0; i < (outLen); i++) {
									transout << timeShift + timeAscii[iterStck][i] << "\t"
									<< epicStack[iterStck][0] << "\t"
									<< radialRFAscii[iterStck][i]
									<< endl;
								}
								transout  << ">" << endl;
								
								
								// Output File in Depth Domain for Second File ****************************** 2.0 ....
								for (int i = 0; i < (outLen); i++) {
									transMigOut << depthAscii[iterStck][i] << "\t"
									<< epicStack[iterStck][0] << "\t"
									<< radialRFAscii[iterStck][i]
									<< endl;
								}
								transMigOut  << ">" << endl;
								
							}else{
								
								// Normal Output With Time - RFamp - Distance Triplet
								for (int i = 0; i < (outLen); i++) {
									transout << timeShift + timeAscii[iterStck][i] << "\t"
									<< epicStack[iterStck][0] << "\t"
									<< radialRFAscii[iterStck][i]
									<< endl;
								}
								transout  << ">" << endl;
								
							}
							
							
						}
						
					}
				}
				
				
			}
			/************************** END EPIC. DUMP HERE *******************************/
		}

		
	}else {
		
	}
	
	
	// Print binary statistics here ...
	if(setPrintStats){
		
		filePrintStats.append("_BinStats.txt");
		ofstream statsOut(filePrintStats.c_str());
		statsOut << logReport << endl;
		
	}
	 
	
	
}

/* getter Routine. Returns timeVector for  use in external migration routine.. */

MatDoub TraceStack::getTimeAscii(){
	return timeAscii;
}

void TraceStack::setDepthAscii(MatDoub& inDepth){
	depthAscii = inDepth;
	
	//cout << "Depth Dimensions" <<  endl;
	flagDepthMig = true;
}

