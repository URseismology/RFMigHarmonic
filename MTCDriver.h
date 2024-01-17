/**
 *  @file MTCDriver.h
 *  
 *
 *  @author Tolulope Olugboji 
 *	@date  11/22/10.
 *  @copyright Yale University. All rights reserved.
 *
 *	@brief Driver class: it recieves Noise & Event (R(Q),T,Z(L)) arrays, and using Slepian tapers, it  computes the spectrum estimates.
 *
 *  The spectrum estimates computed are:
 *	1. Noise Spectrum
 *	2. Event Spectrum
 *
 *  Using this spectrum estimates, the squared coherence estimates  are used to compute the Receiver functions using the method of  J.Park and V.Levin 2000 (http://earth.geology.yale.edu/~jjpark/Park_Levin_2000.pdf)
 *
 */
#include "nr3.h"
#include "fourier.h"
#include "spectrum.h"


class MTCDriver  {
public:
	/* Constructor followed by helper functions. 
	   I overload constructor so that the TraceStack class can use the inversefft
	   Routine, Reuseability is good, but is the unity of MTCDriver respected
	   Also there is the potential unitialized objects might crash class.
	   Investigate.
	 */
	MTCDriver(); //Enables using helper functions 
	MTCDriver(int pos, Mat3DDoub &N, Mat3DDoub &E, Int p, Int k, float f, Doub dt, int nse, int nevt);
	void realSpectrum(VecDoub_IO &a, Int n,VecDoub_IO &spec );
	void cmplxfft(VecDoub_IO &a, Int n, VecComplex &specCmplx);
	void icmplxfft( VecComplex &specCmplx, VecDoub_IO &a, Int n);
	void printVecDoub(VecDoub_IO &a);
	void printVecCmplx(VecComplex &a);
	void TaperRF();
	void writeRFdata(double , double, VecDoub_IO &a);
	
	/* The Numerical workhorse. This routine computes the RFs and 
	 * Associated functions
	 */
	void computeRFs();

	// Updates RFs & Coherence by Record
	void addRF(Mat3DCmplx &, int);
	void addCoher(Mat3DDoub &, int);
	
	// Updates spectrum and coherencd
	void addCorrVH(Mat3DDoub &, int);
	void addSpecAll(Mat3DDoub &, int);
	
	//Write Function that dumps spectrum object into file as update
	void updateSpectrum(ofstream &, int, int);
	void updateCoher(ofstream &, int, int);
	
private:
	/* These are the object states */
	MatDoub NoiseTable, EventTable;		  // resident copy of noise & event
	
	/* Reciever Function[RF], Coherence, Spectrums and UncertaintyRF */
	MatDoub Coher, NoiseSpectrum, EventSpectrum, DeltaRF;
	MatComplex RF;
	
	enum RFcmp {
		R,
		T
	};
	
	enum Components {
		Vertical,
		Radial,
		Transverse
	};
	

	int P, K; // Slepian Parameters... changed only in the calling routine
	int recLen; // Record Length, useful for nyquist determination
	int ncmp, FreqNyq;
	double FreqMax, FreqRlg;
	double DeltaT, boostNoise;
	int nNoise, nEvent; // Size of noise & event, Tapers & Boosting.
	
};
