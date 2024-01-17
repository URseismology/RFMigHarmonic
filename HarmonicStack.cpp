/**
 *  @file HarmonicStack.cpp
 *
 *
 *  @author Tolulope Olugboji
 *  @date Begin November 15, 2012.
 *  @version Updated on July 13, 2013
 *  @brief Computes Harmonic decompostion of RFs in Freq domain Cannibalized from TracesStack: 12/26/10
 *
 *
 *
 *  @copyright 2010 Yale University. All rights reserved.
 *
 *  This class is a specialist class that stacks the Reciever functions
 *  Using Harmonic Decomposition of the back Azimuths, while also using the coherence estimates
 *  to regularize the decomposition in the frequency domain.
 *
 *   It has an operator that allows the printing
 *  of the final stack as files or stdout? In the process of designing this.
 *
 *
 *  Updated on July 13, 2013 --- Updated print function to bias time axis with time delay
 *												for moving window migration display.
 *  @bug no known bugs
 */

#include "HarmonicStack.h"
#include "ludcmp.h"
#include "ran.h"
#include "MTCDriver.h"
#include <string>

HarmonicStack::HarmonicStack(int recLen, int nCmps, int nRecs,
					   VecDoub &A, VecDoub &E, float f, float delT, float pre, 
					   float post, int nBoot)
: Azim(A), Epic(E), zero(0.0,0.0), FreqMax(f), DeltaT(delT), preTime(pre), postTime(post), outLen(0),
binSze(3),ntimesBoot(nBoot), dimCmps(nCmps),dimTotRecs(nRecs)
{
	
	/* Initialize your 3Dimensional arrays. Extension: 
	 Write data structure driver in nr3.h
	 Here I allocate the same memory structure ..
	 Solution Point: Tentative problem: I can't initialize NRMat3D look into
	 this. For now put this initialization in stackAzim routine*/
	//Mat3DCmplx RFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	//Mat3DDoub sigmaRFTrace( RF.dim1(), RF.dim2(), RF.dim3() );
	
	
	FreqNyq = recLen/2 + 1;
	FreqRlg = 1.0 / (DeltaT * recLen);
	nFreqMax = (FreqMax / FreqRlg) + 1;
	Pih = M_PI/2.0;
	iZero = 0;
	nPad = recLen;
	padNyq = nPad/2 + 1;
	//outStream << iZero << endl;
	
	logReport.append("No of records stacked: ");
	int nRecords = Azim.size();
	logReport.append( to_string( nRecords ) );
	logReport.append( "\n" );
        
    
	// Initialize BootStrap Random Vector
	MatDoub LocalRand(nBoot,nRecords,0.0);
	RandomSet4Boot = LocalRand;
    
    // Initialize Back-Azimuth Functions, RecordWide
    
	dimBckAzimCoeff = 10;
	constBazCoeffBias = 3.0;
	
	/* BackAzim Function: Dynamically Populated By Recorded, 
	   Write A Helper Function That Uses:
						  1.Azim@RecIndex, 3. RadialTransFlag
	   And Populates the 
	*/
	BckAzimFuncQ.assign(dimBckAzimCoeff,0.0);
	BckAzimFuncT.assign(dimBckAzimCoeff,0.0);
	
	/* Initialize the Transform matrix with the right dimensions.
	   Consult the design docs, RFDemo to understand this.
	*/
	MatDoub LocalTrans(dimBckAzimCoeff,dimBckAzimCoeff,0.0);
	BckAzimTransMtrix = LocalTrans;	
	
	
	MatComplex LocalData(dimBckAzimCoeff,nPad, zero);
	BckAzimData = LocalData;  meanBckAzimCoeffCmplx = LocalData; 
	
	
	initializeBckAzimTransMtrix();
	initializeBckAzimData();
	
	if (ntimesBoot <= 0) {
		ntimesBoot = 1;
		flagBoot = false;
	}else{
		flagBoot = true;
	}
	
}




/* Use RFs and deltaRFs in Freq Domain to Regress for BackAzim Coefficients
 */
void HarmonicStack::regressNBootTimes(const Mat3DCmplx& RFTrace, const Mat3DDoub& deltaRFTrace) {
    logReport.append("Regressing Using RFs and RandomSet \n ");
	
	// Define the LocalBootStrap 3D Matrices. Do bounds checking though.
	Mat3DCmplx nBootBckAzimCoeffCmplx (ntimesBoot, dimBckAzimCoeff,nPad);
	int ntotRec = Azim.size();
    generateRandomSet(ntimesBoot, ntotRec);  // Compute BootStrap nBoot times.
	//dumpMatrix(RandomSet4Boot);
	
    for (int iBoot = 0; iBoot < ntimesBoot; iBoot++) {
        
		//cout << iBoot << endl; // noOfBootStrap Resampling to Be Done
		
		
        /* For each frequency, stack all records, and regress for bckAzim coefficients
		 This new update doesn't save the backAzim transform matrix since it can be 
		 reproduced very easily from the original RFinFreq data. Review this!!!!???!!
		 */
        for (int iFreq = 0; iFreq < nFreqMax; iFreq++) {
            
			// For each frequency re-initialize your transform Matrix, and redo-regression!!
			if (iFreq > 0) {
				initializeBckAzimTransMtrix();
				initializeBckAzimData();
			}
			
            for (int itercmp = 0; itercmp < 2; itercmp++) {
				for (int iRec = 0; iRec < ntotRec; iRec++) {
					
					int rndRec = RandomSet4Boot[iBoot][iRec];
					double rndAzim = Azim[rndRec];
					
					// Pick out reqularization value: sigmaRF[iFreq]
					//cout << "Record No. in Rand Set: " << rndRec << endl;
					
					double ithSigmaRF = pow(deltaRFTrace[rndRec][itercmp][iFreq], 2.0);
					Complex ithRF = RFTrace[rndRec][itercmp][iFreq];
					populateBckAzimFunc(rndAzim, itercmp);
					
					//Update each bck azimuth component for each record.
					updateBckAzimTransMtrix(ithSigmaRF,  itercmp);
					updateBckAzimData(ithRF,ithSigmaRF, iFreq, itercmp);
					
				}
                
            }
            
			VecDoub realData; getReal(BckAzimData, iFreq, realData);
			VecDoub imagData; getImag(BckAzimData, iFreq, imagData);
            
			VecDoub realCoeffVec(dimBckAzimCoeff, 0.0);
			VecDoub imagCoeffVec(dimBckAzimCoeff, 0.0);
			
			// do actual regression here by alu decomposition.
			// scope operator to free data vectors. 
			// update regression coefficients.
				
			//cout <<" Real & Imag. Data Vecs" << endl;
            //dumpVecDoub(realData); dumpVecDoub(imagData);
			//dumpMatrix(BckAzimTransMtrix);
			
				LUdcmp alu(BckAzimTransMtrix);
				alu.solve(realData, realCoeffVec);
				alu.solve(imagData, imagCoeffVec);
			
			//cout << " Done LUDCmp" << endl;
			//dumpVecDoub(realCoeffVec); dumpVecDoub(imagCoeffVec);
			storeiBootBazCoeffMatrix(realCoeffVec, imagCoeffVec, nBootBckAzimCoeffCmplx, iFreq, iBoot);
			updateMeanBazCoeffMatrix(nBootBckAzimCoeffCmplx, iFreq, iBoot);
			//cout << "Done Freq  " << iFreq << endl; 
			
			
        }
		
		
    }
	
	// Done with nBootstrap stacking of the regression coefficients
	// I use the mean in Freq domain to invert and then use this to find the variance.
	// at end of bootstrap, scale and square at the same time. Keeps things clean..
	
	int nBazComps = dimBckAzimCoeff/2;   
							// no. of baz. components = 1 for constant 4 for costheta and cos2theta
	
	int npre = preTime / DeltaT, npost = postTime / DeltaT;
	outLen = npre + npost;                 // Variable holds output legth
	meanModelledRFAscii.assign(nBazComps, outLen, 0.0);
	meanUnmodelledRFAscii.assign(nBazComps, outLen, 0.0);
	timeAscii.assign(dimBckAzimCoeff, outLen, 0.0);
	
	// Here is the Mean Bootstrap in time domain:
	buildTimeDomainBazCoeffs(meanBckAzimCoeffCmplx, meanModelledRFAscii, meanUnmodelledRFAscii);
		
	
	// Store sqrt of RF variance in time domain, using mean value above. If nBoot == 0 Skip Step.
	devModelledRFAscii.assign(nBazComps, outLen, 0.0);
	devUnmodelledRFAscii.assign(nBazComps, outLen, 0.0);
	ithModelledRFAscii.assign(nBazComps, outLen, 0.0);
	ithUnmodelledRFAscii.assign(nBazComps, outLen, 0.0);
	
	MatComplex iBootFreqBazCoeff(nBootBckAzimCoeffCmplx.dim2(), nBootBckAzimCoeffCmplx.dim3());
	
	for (int iBoot = 0; iBoot < ntimesBoot; iBoot++) {
		// Pick out the iBoot Complex Matrix ... 
		
		// ****** write member functioin that takes 3D Matrix and Initiliazes 2D Matrix:
		make2Dfrom3DMatrix(iBootFreqBazCoeff, nBootBckAzimCoeffCmplx, iBoot);
		// ***** nBoot is flattened to compute the variance ... interesting right? ***

		buildTimeDomainBazCoeffs(iBootFreqBazCoeff, ithModelledRFAscii, ithUnmodelledRFAscii);
		//dumpMatrix(ithModelledRFAscii);
		
		double cdev = 0;
		for (int iBazComp = 0; iBazComp < nBazComps; iBazComp++) {
			for (int ithLen = 0; ithLen < outLen; ithLen++) {
				
				cdev = ithModelledRFAscii[iBazComp][ithLen] -  meanModelledRFAscii[iBazComp][ithLen];
				
				//cout << "Dev:  " << cdev << endl;
				devModelledRFAscii[iBazComp][ithLen] = devModelledRFAscii[iBazComp][ithLen]
				+ pow(cdev, 2.0);
				
				cdev = ithUnmodelledRFAscii[iBazComp][ithLen] -  meanUnmodelledRFAscii[iBazComp][ithLen];
				devUnmodelledRFAscii[iBazComp][ithLen] = devUnmodelledRFAscii[iBazComp][ithLen]
				+ pow(cdev, 2.0);
				
			}
		}
		
	}

	

	
}

/* Populate ithRandomSet of NBootStrap Samples. Called by 'regress module' 
 nBoot Times
 */
void HarmonicStack::generateRandomSet(int nBoot, int ntotRec){
    // srand( time(NULL));
	
	for (int iBoot = 0; iBoot < nBoot; iBoot++) {
		Ran ranset(time(NULL) + iBoot);
		for (int iRand =0; iRand < ntotRec; iRand++) {
			//Random numbers from 1 to ntotRec
			RandomSet4Boot[iBoot][iRand] =  ranset.int64() % (ntotRec);
			//RandomSet4Boot[iRand] = rand() % ntotRec;
		}
	}

}	

void HarmonicStack::initializeBckAzimTransMtrix(){
	/*
	 first & second dimensions are the square matrix of the dimBckAzimCoeff.
	 third dimension is the number of frequencies.
	 */ 
	
	for (int first = 0; first < BckAzimTransMtrix.nrows(); first++) {
		for (int second = 0; second < BckAzimTransMtrix.ncols(); second++) {
				BckAzimTransMtrix[first][second] = 0.0;
			
		}
	}

	logReport.append("Transform Matrix Initialized. \n ");
}

void HarmonicStack::initializeBckAzimData(){
	/*
	 first & second dimensions are the square matrix of the dimBckAzimCoeff.
	 third dimension is the number of frequencies.
	 */
	//cout << "Rows and Cols: " << BckAzimData.nrows() << " " << BckAzimData.ncols() << endl;
	for (int first = 0; first < BckAzimData.nrows(); first++) {
		for (int second = 0; second < BckAzimData.ncols(); second++) {
			BckAzimData[first][second] = zero;
			
		}
	}
	logReport.append("Complex Data Matrix Initialized. \n ");
}

/* Populate The BackAzimuthFunction at a particular AzimbyRecIndx
 */
void HarmonicStack::populateBckAzimFunc(double ithAzim, int RadTransFlag) {
    double AzimRad = ithAzim * M_PI/180; // Convert Degree to Radians
    if (RadTransFlag == 0) {
        BckAzimFuncQ[0] = constBazCoeffBias;
        BckAzimFuncQ[1] = cos(AzimRad);
        BckAzimFuncQ[2] = sin(AzimRad);
        BckAzimFuncQ[3] = cos(2.0*AzimRad);
        BckAzimFuncQ[4] = sin(2.0*AzimRad);
        BckAzimFuncQ[5] = 0.0;
        BckAzimFuncQ[6] = cos(AzimRad);
        BckAzimFuncQ[7] = sin(AzimRad);
        BckAzimFuncQ[8] = cos(2*AzimRad);
        BckAzimFuncQ[9] = sin(2*AzimRad);
    }else{
        BckAzimFuncT[0] = 0.0;
        BckAzimFuncT[1] = sin(AzimRad);
        BckAzimFuncT[2] = -cos(AzimRad);
        BckAzimFuncT[3] = sin(2.0*AzimRad);
        BckAzimFuncT[4] = -cos(2.0*AzimRad);
        BckAzimFuncT[5] = constBazCoeffBias;
        BckAzimFuncT[6] = -sin(AzimRad);
        BckAzimFuncT[7] = cos(AzimRad);
        BckAzimFuncT[8] = -sin(2*AzimRad);
        BckAzimFuncT[9] = cos(2*AzimRad);
        
    }
	logReport.append("Transform Functions Re-Initialized. \n ");
}

void HarmonicStack::updateBckAzimTransMtrix(double ithSigmaRF, int flagQT){
	
	if (flagQT == 0) {
		for (int first = 0; first < dimBckAzimCoeff; first++) {
			for (int second = 0; second < dimBckAzimCoeff; second++) {
				//cout << "Q " << first << second  <<endl;
				//cout << BckAzimTransMtrix.ncols()<< endl;
				BckAzimTransMtrix[second][first] = 
				BckAzimTransMtrix[second][first] + BckAzimFuncQ[second]*BckAzimFuncQ[first]/(ithSigmaRF); 
			}
		}
	}else{
		for (int first = 0; first < dimBckAzimCoeff; first++) {
			for (int second = 0; second < dimBckAzimCoeff; second++) {
				BckAzimTransMtrix[second][first] = 
				BckAzimTransMtrix[second][first] + BckAzimFuncT[second]*BckAzimFuncT[first]/(ithSigmaRF); 
			}
		}
	}
	
	//logReport.append("BckAzimMatrix Transform Functions Updated. \n ");
	
}

void HarmonicStack::updateBckAzimData(Complex ithRF, double ithSigmaRF, int iFreq, int flagQT){
	
	if (flagQT == 0) {
		for (int first = 0; first < dimBckAzimCoeff; first++) {
			Complex storeLast = BckAzimData[first][iFreq];
			//cout << storeLast ;
			double ithBckFunc = BckAzimFuncQ[first]/ithSigmaRF;
			BckAzimData[first][iFreq] = storeLast + ithRF*ithBckFunc; 
			
		}
	}else{
		for (int first = 0; first < dimBckAzimCoeff; first++) {
			Complex storeLast = BckAzimData[first][iFreq];
			//cout << storeLast ;
			double ithBckFunc = BckAzimFuncT[first]/ithSigmaRF;
			BckAzimData[first][iFreq] = storeLast + ithRF*ithBckFunc;  
			
		}
	}
	//cout << endl;
	//logReport.append("Data Matrix Updated. \n ");
	
}

void HarmonicStack::getReal(const MatComplex& BckAzimData, int iFreq, VecDoub&outReal){
	VecDoub realData(dimBckAzimCoeff, 0.0);
	for (int iDim = 0; iDim < dimBckAzimCoeff; iDim++) {
		realData[iDim] = real(BckAzimData[iDim][iFreq]);
		//cout << realData[iDim] << endl;
	}
	outReal = realData;
}

void HarmonicStack::getImag(const MatComplex& BckAzimData, int iFreq, VecDoub& outImag){
	VecDoub imagData(dimBckAzimCoeff, 0.0);
	for (int iDim = 0; iDim < dimBckAzimCoeff; iDim++) {
		imagData[iDim] = imag(BckAzimData[iDim][iFreq]);
        //cout << imagData[iDim] << endl;
	}
	outImag = imagData;
}

void HarmonicStack::storeiBootBazCoeffMatrix(const VecDoub& real, const VecDoub& imag, Mat3DCmplx& nBootBckAzimCoeffCmplx, int iFreq, int iBoot){
	for (int iRow = 0; iRow < dimBckAzimCoeff; iRow++) {
		if ( iRow == 0 || iRow == 5 ) {
			nBootBckAzimCoeffCmplx[iBoot][iRow][iFreq] = Complex(real[iRow] * constBazCoeffBias, imag[iRow]*constBazCoeffBias);
		}else{
			nBootBckAzimCoeffCmplx[iBoot][iRow][iFreq] = Complex(real[iRow], imag[iRow]);			
		}
		
	}
}

void HarmonicStack::updateMeanBazCoeffMatrix(Mat3DCmplx& nBootBckAzimCoeffCmplx, int iFreq, int iBoot){
	double scaleMean = 1.0;

	for (int iRow = 0; iRow < dimBckAzimCoeff; iRow++) {
		
		meanBckAzimCoeffCmplx[iRow][iFreq] = meanBckAzimCoeffCmplx[iRow][iFreq] + nBootBckAzimCoeffCmplx[iBoot][iRow][iFreq];
	}
	
	// Check If you have summed up all the nBootStrap Realizations.
	// If so, scale the array by the right measure. Done only at the last call. 
	
	if (iBoot == (ntimesBoot - 1)) {
		scaleMean = 1.0 / double(ntimesBoot);
		
		for (int iRow = 0; iRow < dimBckAzimCoeff; iRow++) {
			meanBckAzimCoeffCmplx[iRow][iFreq] = meanBckAzimCoeffCmplx[iRow][iFreq] * scaleMean;
		}
	}
}

void HarmonicStack::make2Dfrom3DMatrix(MatComplex& twoD, Mat3DCmplx& threeD, int indx){
	int dim1 = threeD.dim2();
	int dim2 = threeD.dim3();
	
	for (int iDim1 = 0; iDim1 < dim1; iDim1++) {
		for (int iDim2 = 0; iDim2 < dim2; iDim2++) {
			twoD[iDim1][iDim2] = threeD[indx][iDim1][iDim2];
		}
	}
	
	
}

void HarmonicStack::TaperBazCoeffs(MatComplex& freqDomainBazCoeffs, int nFMax){
	Complex ithBazVal; // Store Before Taper
	
	
	for (int iterBazCmps = 0; iterBazCmps < dimBckAzimCoeff; iterBazCmps++) {
		//cout << "Begin Taper: " << endl;
		for (int iternf = 0; iternf < nFMax; iternf++) {
			//cout << "Taper Func "<< ( double(iternf)/double(nFMax) ) << endl;
			Hann = pow( cos( Pih * ( double(iternf)/double(nFMax) ) ), 2.0 );
			ithBazVal = freqDomainBazCoeffs[iterBazCmps][iternf];
			freqDomainBazCoeffs[iterBazCmps][iternf] = Hann * ithBazVal;
		}
	}
	
}


void HarmonicStack::buildTimeDomainBazCoeffs(MatComplex& freqDomainBazCoeffs, MatDoub& modelledRFAscii, MatDoub& unmodelledRFAscii){
	// This code is gonna be repeated often. Put it in a Method Function. 
	// ComputeIRF in Time Domain. OR something like that.
	TaperBazCoeffs( freqDomainBazCoeffs, nFreqMax);
	for (int itercmp = 0; itercmp < 2; itercmp++) {
		
		// Dump RF here.. big bug. Yish!!
		/*for (int i = 0; i < RFtimebyAzim.size(); i++) {
		 debugFile << abs(RFtimebyAzim[i]) << "  " << endl;
		 } */
		//Build Reciever Function in time domain? 
		//Causative & Non Causative Build
		if (itercmp == 0) {
			for (int iterBazCmps = 0; iterBazCmps < 5; iterBazCmps++) {
				
				VecComplex RFfreq(freqDomainBazCoeffs.ncols(), freqDomainBazCoeffs[iterBazCmps]);
				VecDoub_IO RFtimebyAzim;
				MTCDriver invertStack;
				
				invertStack.icmplxfft(RFfreq, RFtimebyAzim, padNyq);
				
				buildAsciiOut(preTime, postTime, RFtimebyAzim, modelledRFAscii, timeAscii,
							  iterBazCmps);
			}
			
		}else {
			for (int iterBazCmps = 5; iterBazCmps < 10; iterBazCmps++) {
				
				
				
				VecComplex RFfreq(freqDomainBazCoeffs.ncols(), freqDomainBazCoeffs[iterBazCmps]);
				VecDoub_IO RFtimebyAzim;
				MTCDriver invertStack;
				
				invertStack.icmplxfft(RFfreq, RFtimebyAzim, padNyq);
				
				int biasIndex = iterBazCmps - 5;  // bias Index of Unmodelled For Storage.
				buildAsciiOut(preTime, postTime, RFtimebyAzim, unmodelledRFAscii, timeAscii,
							  biasIndex);
			}
		}
	}
	
}


void HarmonicStack::dumpMatrix(const MatDoub& mat){
	for (int iRow = 0; iRow < mat.nrows(); iRow++) {
		for (int iCol = 0; iCol < mat.ncols(); iCol++) {
			cout << mat[iRow][iCol] << "  ";
		}
		cout << endl;
	}
}

void HarmonicStack::dumpVecDoub(const VecDoub& vec){
	for (int iRow = 0; iRow < vec.size(); iRow++) {
			cout << vec[iRow] << "  ";
		}
		cout << endl;
}


/* This member method helps build the ascii output before using the out function 
   All RFs in this section are in the time domain, what I need to do is get more
   printfns to print RFs in the frequency domain for test cases ...
 */
void HarmonicStack::buildAsciiOut(double pretime, double posttime, VecDoub_IO &RFin, 
							  MatDoub &RFout, MatDoub &timeOut, int ithBazCmps)
{
	int recLen = RFin.size();
	int npre = pretime / DeltaT, npost = posttime / DeltaT;
 
	double scaleRF = 2.0 * FreqNyq / (FreqMax/FreqRlg);
	
	for (int i = 1; i <= (npre+npost) ; i++) {
		timeOut[ithBazCmps][i-1] =  -npre*DeltaT + (i-1)*DeltaT ;
	}
	
	for (int i =1 ; i <= npost; i++) {
		RFout[ithBazCmps][npre+i-1] = scaleRF * RFin[i-1];
	}
	
	for (int i = 1; i <= npre; i++) {
		RFout[ithBazCmps][npre-i] = scaleRF * RFin[recLen - i-1];
	}
	
		
}

	/* Utility method for printing the ascii traces into a file .. 
	 New Adaptation: I need to convert the frequency decimation to time decimation
	 */
void HarmonicStack::print(char* outfn, int fileFlag, int timeShift){
	
	// Time - Constant And BAZ Components
	string fileMeanModelled(outfn);
	string  fileMeanUnModelled(outfn);
	
	string fileDevMaxModelled(outfn);
	string  fileDevMaxUnModelled(outfn);
	
	string fileDevMinModelled(outfn);
	string  fileDevMinUnModelled(outfn);
	
	// Time - Transverse Vector Functions
	string fileMeanModelledVector(outfn);
	string fileMeanUnModelledVector(outfn);
	
	string fileDevMaxModelledVector(outfn);
	
	string fileDevMinModelledVector(outfn);
	
	// Depth - Transverse Vector Function
	string fileMeanModelledVectorD(outfn);
	
	string fileDevMaxModelledVectorD(outfn);
	
	string fileDevMinModelledVectorD(outfn);
	
	
	if (fileFlag == 0) {
		
		// Mean ... Components
		fileMeanModelled.append(".Modelled.mean.xyz");
		fileMeanUnModelled.append(".UnModelled.mean.xyz");
		
		// Variances ...
		fileDevMaxModelled.append(".Modelled.devMax.xyz");
		fileDevMaxUnModelled.append(".UnModelled.devMax.xyz");
		
		fileDevMinModelled.append(".Modelled.devMin.xyz");
		fileDevMinUnModelled.append(".UnModelled.devMin.xyz");
		
		
		// Mean ... Transverse Vectors
		fileMeanModelledVector.append(".Modelled.mean.rt");
		fileMeanUnModelledVector.append(".UnModelled.mean.rt");
		
		// Mean ... Transverse Vectors - Depth
		fileMeanModelledVectorD.append(".Modelled.mean.rt_D");
		
		
		// Variance ... Transverse Vectors
		fileDevMaxModelledVector.append(".Modelled.devMax.rt");
		fileDevMinModelledVector.append(".Modelled.devMin.rt");
		
		// Variance ... Transverse Vectors - Depth
		fileDevMaxModelledVectorD.append(".Modelled.devMax.rt_D");
		fileDevMinModelledVectorD.append(".Modelled.devMin.rt_D");
		
	}else {

	}

	
	ofstream meanModelledOut(fileMeanModelled.c_str());
	ofstream meanUnModelledOut(fileMeanUnModelled.c_str());
	
	ofstream devMaxModelledOut(fileDevMaxModelled.c_str());
	ofstream devMaxUnModelledOut(fileDevMaxUnModelled.c_str());
	
	ofstream devMinModelledOut(fileDevMinModelled.c_str());
	ofstream devMinUnModelledOut(fileDevMinUnModelled.c_str());
	
	
	// Files For Vector Rose Plots - No Migration  ....
	ofstream meanVectorOut(fileMeanModelledVector.c_str());
	ofstream meanVectorUnModelledOut(fileMeanUnModelledVector.c_str());
	
	ofstream devMaxVectorOut(fileDevMaxModelledVector.c_str());
	ofstream devMinVectorOut(fileDevMinModelledVector.c_str());
	
	// Files For Vector Rose Plots - Migration ......
	ofstream meanVectorOutD(fileMeanModelledVectorD.c_str());
	ofstream devMaxVectorOutD(fileDevMaxModelledVectorD.c_str());
	ofstream devMinVectorOutD(fileDevMinModelledVectorD.c_str());

	
	cout << "Opened file: " << outfn << endl;
	logReport.append("Files Opened  \n ");
	
	bool outdata = true; // Check here for file, or outside the fucction?
	
	if (outdata) {
		// Go through the stack and dump in file, first dimension 
		int nBazCmps = meanModelledRFAscii.nrows();   // Hopefully, this holds no. of Coefficients
		
		//Here is the dump code snippet . I use here again the file flag
		if (fileFlag == 0) {
			
			// Dump Modelled RFs in time domain here ...
			for (int iterBazCmps = 0; iterBazCmps < nBazCmps; iterBazCmps++) {
				
				// Compute Axis Label. Arbitrary Units, Just For Display
				int cmpAxisUnit = (4-iterBazCmps)*10;
				
				for (int itercmp = 0; itercmp < 2; itercmp++) {
					
					if (itercmp == 0) {
						
						// Dump all modelled files. 1.Mean 2.DevMax and 3.DevMin
						for (int i = 0; i < (outLen); i++) {
							
							// Mean Here
							meanModelledOut << timeShift + timeAscii[iterBazCmps][i] << "\t"
							<< cmpAxisUnit << "\t"
							<< meanModelledRFAscii[iterBazCmps][i] 
							<< endl;
							
							
							//DevMax Here
							devMaxModelledOut << timeShift + timeAscii[iterBazCmps][i] << "\t"
							<< cmpAxisUnit << "\t"
							<< meanModelledRFAscii[iterBazCmps][i] + 
							                sqrt(devModelledRFAscii[iterBazCmps][i]/double(ntimesBoot)) << endl;
							
							//DevMin Here
							devMinModelledOut << timeShift + timeAscii[iterBazCmps][i] << "\t"
							<< cmpAxisUnit << "\t"
							<< meanModelledRFAscii[iterBazCmps][i] - 
									sqrt(devModelledRFAscii[iterBazCmps][i]/double(ntimesBoot)) << endl;
						}
						meanModelledOut  << ">" << endl;
						devMaxModelledOut << ">" << endl;
						devMinModelledOut << ">" << endl;
						
						
						
					}else {
						// I should either change name or find a way to make this 
						// manageable ...
						
						// Dump all unmodelled files. 1. Mean 2. DevMax and 3. DevMin
						for (int i = 0; i < (outLen); i++) {
							
							// Mean Here
							meanUnModelledOut << timeShift + timeAscii[iterBazCmps][i] << "\t"
							<< cmpAxisUnit << "\t"
							<< meanUnmodelledRFAscii[iterBazCmps][i] 
							<< endl;
							
							
							//DevMax Here
							devMaxUnModelledOut << timeShift + timeAscii[iterBazCmps][i] << "\t"
							<< cmpAxisUnit << "\t"
							<< meanUnmodelledRFAscii[iterBazCmps][i] + 
									sqrt(devUnmodelledRFAscii[iterBazCmps][i]/double(ntimesBoot))
							<< endl;
							
							//DevMin Here
							devMinUnModelledOut << timeShift + timeAscii[iterBazCmps][i] << "\t"
							<< cmpAxisUnit << "\t"
							<< meanUnmodelledRFAscii[iterBazCmps][i] - 
									sqrt(devUnmodelledRFAscii[iterBazCmps][i]/double(ntimesBoot))
							<< endl;

						}
						meanUnModelledOut  << ">" << endl;
						devMaxUnModelledOut  << ">" << endl;
						devMinUnModelledOut  << ">" << endl;
						
					}
					
				}
				}
			
			// Dump Vector RFs in time domain here ...
			//** Modelled First ...
			double dipX, dipY, horX, horY;
			double tempAngle, tempMag;
			
			vector<double> dipAngle, dipMag;
			vector<double> horAngle, horMag;
			
			//** umModelled Next ...
			double dipX_U, dipY_U, horX_U, horY_U;
			
			vector<double> dipAngle_U, dipMag_U;
			vector<double> horAngle_U, horMag_U;
			

			
			// Build vectors here at each time slice ..
			for (int i = 0; i < (outLen); i++) {
				
				//** Modelled First
				dipX = meanModelledRFAscii[2][i];
				dipY = meanModelledRFAscii[1][i];
				
				tempAngle = atan2( dipY, dipX ) * 180/M_PI;
				tempMag = sqrt(dipX*dipX + dipY*dipY);
				
				dipAngle.push_back(tempAngle);
				dipMag.push_back(tempMag);
				
				//** UnModelled Next ...
				dipX_U = meanUnmodelledRFAscii[2][i];
				dipY_U = meanUnmodelledRFAscii[1][i];
				
				tempAngle = atan2( dipY_U, dipX_U ) * 180/M_PI;
				tempMag = sqrt(dipX_U*dipX_U + dipY_U*dipY_U);
				
				dipAngle_U.push_back(tempAngle);
				dipMag_U.push_back(tempMag);
				
				
				//** Modelled First
				horX = meanModelledRFAscii[4][i];
				horY = meanModelledRFAscii[3][i];
				
				tempAngle = atan2( horY, horX ) * 180/M_PI;
				tempMag = sqrt(horX*horX + horY*horY);
				
				horAngle.push_back(tempAngle);
				horMag.push_back(tempMag);
				
				//** UnModelled Next ...
				horX_U = meanUnmodelledRFAscii[4][i];
				horY_U = meanUnmodelledRFAscii[3][i];
				
				tempAngle = atan2( horY_U, horX_U ) * 180/M_PI;
				tempMag = sqrt(horX_U*horX_U + horY_U*horY_U);
				
				horAngle_U.push_back(tempAngle);
				horMag_U.push_back(tempMag);
								
				
			}
			
			// Save Mean Vector Here
			for (int anisoType = 0; anisoType < 4; anisoType++) {
				
				switch (anisoType) {
					case 0:
						for (int i = 0; i < outLen; i++) {
							
							//** Modelled ...
							meanVectorOut << timeShift + timeAscii[0][i] << "\t"
							<< anisoType << "\t" << dipMag[i] << endl;
							
							//** UnModelled ...
							meanVectorUnModelledOut << timeShift + timeAscii[0][i] << "\t"
							<< anisoType << "\t" << dipMag_U[i] << endl;
							
						}
						break;
					case 1:
						for (int i = 0; i < outLen; i++) {
							
							//** Modelled ...
							meanVectorOut << timeShift + timeAscii[0][i] << "\t"
							<< anisoType << "\t" << dipAngle[i] << endl;
							
							//** UnModelled ...
							meanVectorUnModelledOut << timeShift + timeAscii[0][i] << "\t"
							<< anisoType << "\t" << dipAngle_U[i] << endl;
							
						}
						break;
					case 2:
						for (int i = 0; i < outLen; i++) {
							
							//** Modelled ...
							meanVectorOut << timeShift + timeAscii[0][i] << "\t"
							<< anisoType << "\t" << horMag[i] << endl;
							
							//** UnModelled ...
							meanVectorUnModelledOut << timeShift + timeAscii[0][i] << "\t"
							<< anisoType << "\t" << horMag_U[i] << endl;
						}
						break;
					case 3:
						for (int i = 0; i < outLen; i++) {
							
							//** Modelled ...
							meanVectorOut << timeShift + timeAscii[0][i] << "\t"
							<< anisoType << "\t" << horAngle[i] << endl;
							
							//** UnModelled ...
							meanVectorUnModelledOut << timeShift + timeAscii[0][i] << "\t"
							<< anisoType << "\t" << horAngle_U[i] << endl;
						}
						break;
					default:
						break;
				}
				meanVectorOut << ">" << endl;
				meanVectorUnModelledOut << ">" << endl;
			}

			
			// Dump Vector RFs in depth domain here ...
				
				
			
			/************************** END HARMONIC DUMP HERE *******************************/
		} else {
			/************************** PUT FULL EXPANSION HERE *******************************/
		}

		
	}else {
		
	}
	 
	
	
}

/* getter Routine. Returns timeVector for  use in external migration routine.. */

MatDoub HarmonicStack::getTimeAscii(){
	return timeAscii;
}

void HarmonicStack::setDepthAscii(MatDoub& inDepth){
	depthAscii = inDepth;
	
	//cout << "Depth Dimensions" <<  endl;
	flagDepthMig = true;
}

void HarmonicStack::saveSpecMeanBAZcmps(int cmps, char* outfn){
	// dump meanBaz component, use zero for icmp to print modelled constant term.
	string fileMeanCnstBAZcmp(outfn);
	fileMeanCnstBAZcmp.append(".constBAZspec.txt");
	ofstream specOut( fileMeanCnstBAZcmp.c_str() );
	float realOut, imagOut;
	
	// print to Nyquist or to Rayleigh ... Nyquist.
	//	int nFreqMax = (FreqNyq/ FreqRlg) + 1;
	
	for (int iFreq = 0; iFreq < nFreqMax; iFreq++){
		
		realOut = real( meanBckAzimCoeffCmplx[cmps][iFreq] );
		imagOut = imag( meanBckAzimCoeffCmplx[cmps][iFreq] );
		
		specOut << iFreq*FreqRlg << "\t"  << realOut  << " \t "
		<< imagOut << endl;
		
	}
	
}