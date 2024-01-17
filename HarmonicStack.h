/**
 *  @file HarmonicStack.h
 *  
 *
 *  @author Tolulope Olugboji   
 *  @date Begin November 15, 2012.
 *  @version 1.0 
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

#include "nr3.h"
#include <sstream>

class HarmonicStack{
private:
	VecDoub Azim, Epic;
	Mat3DCmplx stackRFTrace; 
	//Mat3DDoub  sigmaRFTrace;
	MatDoub  meanModelledRFAscii, meanUnmodelledRFAscii, timeAscii, devModelledRFAscii, devUnmodelledRFAscii, ithModelledRFAscii, ithUnmodelledRFAscii;
	
	MatDoub depthAscii; bool flagDepthMig; // This is value time transformed to depth using velocity migration.
	
	
	bool flagBoot;	// set false if you want to ignore Bootstrap resampling
   	
    MatDoub BckAzimTransMtrix;
    VecDoub BckAzimFuncQ, BckAzimFuncT;
    MatComplex  meanBckAzimCoeffCmplx; 
	MatComplex BckAzimData;
    
    // Mat3DCmplx nBootBckAzimCoeffCmplx; // I remove this matrix and I put it in
	// the regressNBootTimes. This is because I need it only temporarily to build
	// the mean and deviation RF Matrices. After I'm done creating that data,
	// I don't need it anymore. I can then send it's reference to helper functions.
	
    MatDoub RandomSet4Boot;
	double constBazCoeffBias;
    
    //float binmin, binmax, bininc;             Not useful in Harmonic stack, since we use all azimuths in the regression.
	float FreqNyq, FreqRlg, FreqMax, DeltaT;
	int nFreqMax;
	float preTime, postTime;
	int outLen, nPad, expo, iZero, padNyq, binSze, ntimesBoot;
    int dimBckAzimCoeff, dimCmps, dimTotRecs;
	Complex zero;
	double Hann;
	double Pih;
	
    

	
public:
	HarmonicStack(int RFLen, int nCmps, int nRecs, VecDoub &Azim, VecDoub &Epic,
			   float fMax, float delT, float pre, float post, int nBoot);
	
	void regressNBootTimes(const Mat3DCmplx& RFTrace, const Mat3DDoub& sigmaRFTrace);
	void generateRandomSet(int nBoot, int ntotRec);
	
	void initializeBckAzimTransMtrix();
	void initializeBckAzimData();
	void populateBckAzimFunc(double ithAzim, int RadTransFlag );
    void updateBckAzimTransMtrix(double ithDelaRF, int flagQT);
    void updateBckAzimData(Complex ithRF, double ithDelaRF, int iFreq, int flagQT);

	void  getReal(const MatComplex& BckAzimData, int iFreq, VecDoub& outReal);
	void  getImag(const MatComplex& BckAzimData, int iFreq, VecDoub& outImag);
	
	void storeiBootBazCoeffMatrix(const VecDoub& real, const VecDoub& imag, Mat3DCmplx& nBootMtrx, int iFreq, int iBoot);
	void updateMeanBazCoeffMatrix(Mat3DCmplx& nBootMtrx, int iFreq, int iBoot);
	void make2Dfrom3DMatrix(MatComplex& twoD, Mat3DCmplx& threeD, int indx);
	
	
	void buildTimeDomainBazCoeffs(MatComplex& freqDomainBazCoeffs, MatDoub& modelledRFAscii, MatDoub& unmodelledRFAscii);
	void TaperBazCoeffs( MatComplex& freqDomainBazCoeffs, int FreqMax);
	
	void dumpMatrix(const MatDoub& BckAzimTransMtrix);
	void dumpVecDoub(const VecDoub& Vector);
	
	void buildAsciiOut(double prelen, double postlen, VecDoub_IO &RFin,
					   MatDoub &RFout, MatDoub &timeOut, int ithAzim);
	void print(char*, int, int); // resolve stream type. Filestream? FileStream Always
	string logReport;

	// MODELLED AFTER TRACE-STACK ...
	/**
	 * @brief  returns a copy of the timeAscii vector for used by external migration routine
	 * @return timeAscii Array
	 * @bug still implementing
	 */
	MatDoub getTimeAscii(); // Return 2D Time Ascii ..
	
	/**
	 * @brief  updates the time to depth migrated data ..
	 * @param  MatDoub receive value ...
	 * @bug still implementing
	 */
	void setDepthAscii(MatDoub& inDepthAscii); // Set 2D Depth Ascii ...
	
	
	template <class T>
	inline std::string to_string (const T& t)
	{
		std::stringstream ss;
		ss << t;
		return ss.str();
	}
	
	// put utility codes here for external processing of spectrum in MATLAB. design code to do resonance muting by dividing out the resonance poles...
	void saveSpecMeanBAZcmps(int chckCmp, char* outfn);
	
};