/**
 *  @file MTCDriver.cpp
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
#include "MTCDriver.h"
#define  VERTSPEC 0
#define  RADSPEC 1
#define  TRANSPEC 2
#define  NOISESPEC 3

MTCDriver::MTCDriver() :P(0), K(0), recLen(0), ncmp(0), nNoise(0), nEvent(0),
						FreqNyq(0), FreqMax(0.0), FreqRlg(0.0), DeltaT(0.0), 
						boostNoise(0.0)
{
	
}

MTCDriver::MTCDriver (int pos, Mat3DDoub &N, Mat3DDoub &E, Int p, Int k, float f, Doub  dt, int nse, int nevt) :
					P(p), K(k), ncmp(3), nNoise(nse), nEvent(nevt),
					FreqMax(f), DeltaT(dt)
			
{
	
	
	Complex zero(0.0,0.0);
	
	if (nNoise < nEvent){
		boostNoise = double (nEvent/nNoise);
	} else {
		boostNoise = 1.0;
	}
	//cout << boostNoise << endl;
	
	NoiseTable.assign(3, N.dim3(), 0);
	EventTable.assign(3, E.dim3(), 0);
	
	recLen = NoiseTable.ncols();
	FreqNyq = recLen/2 + 1;
	FreqRlg = 1.0 / (dt * recLen);
	
	for (int iterall = 0; iterall < N.dim3(); iterall++) {
		for (int itercmp = 0; itercmp < 3; itercmp++) {
			NoiseTable[itercmp][iterall] = N[pos][itercmp][iterall];
			EventTable[itercmp][iterall] = E[pos][itercmp][iterall];
		}
	}
	
	NoiseSpectrum.assign(ncmp, FreqNyq, 0);
	EventSpectrum.assign(ncmp, FreqNyq, 0);
	Coher.assign(ncmp - 1, FreqNyq, 0);
	RF.assign(ncmp-1, FreqNyq, zero);
	DeltaRF.assign(ncmp-1, FreqNyq, 0);
	
	//Compute RFs
	computeRFs();
	
};

void MTCDriver::realSpectrum(VecDoub_IO &a, Int n, VecDoub_IO &spec ){
	int nnf = 0, last = n/2 ;
	realft(a, 1);
	spec.assign(last + 1, 0);
	Int nf = n;
	spec[0] = pow(a[0], 2);
	spec[last] =  pow(a[1], 2);
	
	for (int iter = 2; iter < n-1; iter+=2) {
		nnf = nnf + 1;
		spec[nnf] = pow(a[iter],2) + pow(a[iter+1],2) ;
	}
	//cout << nnf<< endl;;
	
}

void MTCDriver::cmplxfft(VecDoub_IO &a, Int n, VecComplex &spec ){
	int nnf = 0, last = n/2 ;
	realft(a, 1);
	Complex zero(0.0,0.0);
	spec.assign(last + 1,zero);
	Int nf = n;
	spec[0] = Complex(a[0],0);
	spec[last] =  Complex(a[1],0);
	
	for (int iter = 2; iter < n-1; iter+=2) {
		nnf = nnf + 1;
		spec[nnf] = Complex(a[iter],a[iter+1]) ;
	}
	
}

/* I make the routine self contained here, so external classes can use */
void MTCDriver::icmplxfft(VecComplex &spec, VecDoub_IO &newa, Int Fn){
	//reconstruct Vector into realft data format
	// Take n, which is freq. domain length and unwrap it..
	int datLen = 2*(Fn-1);
	
	Int newLen = 2*datLen, nnf = 0;
	VecDoub_IO a;
	a.assign(newLen, 0);
	newa.assign(datLen, 0); // return only the useful values
	
	for (int iter = 0; iter < Fn; iter++) {
		a[nnf] = real(spec[iter]);
		nnf++;
		a[nnf] = imag(spec[iter]);
		nnf++;
	}
	
	for (int iter = Fn-2; iter > 0; iter--) {
		a[nnf] = real(spec[iter]);
		nnf++;
		a[nnf] = -1.0 * imag(spec[iter]);
		nnf++;
	}
	
	four1(a, -1);
	//printVecDoub(a);
	int map = 0;
	
	for (int iter = 0; iter < newLen; iter+=2) {
		map = iter/2;
		newa[map] = (1.0/datLen) * a[iter];      //scale inverse fourier by 1/n
	}
	
	
	
}

void MTCDriver::computeRFs(){
	
	RFcmp RFcmps[2] = {R, T};
	Components Z = Vertical;
	
	VecDoub_IO tempN(FreqNyq), tempE(FreqNyq);
	VecComplex tempEcmplx;
	

	Mat3DCmplx  KEigenSpectraEvent(ncmp, K, FreqNyq) ;
	Slepian pPiTaper(recLen/2, P, K);
	
	// Compute Spectrum for each seismic component ...
	for (int itercmp = 0; itercmp < ncmp ; itercmp++) {
		
		for (int iterK = 0; iterK < K; iterK++) {
			
			// Pick out the data each time before tapering ...
			VecDoub taperNoise(NoiseTable.ncols(), NoiseTable[itercmp]);
			VecDoub taperEvent(EventTable.ncols(), EventTable[itercmp]);
			
			// Taper the Noise and Event Table
			for (int itern = 0; itern < recLen ; itern++) {
				taperNoise[itern] *= pPiTaper.dpss[iterK][itern];
				taperEvent[itern] *= pPiTaper.dpss[iterK][itern];
				
			}
			// Sum spectrum over K
			//cout << recLen << endl;
			
			//Fixed potential bug here by refreshing taperEvent which, without this line, 
			// would be holding the fourier transform, instead of the data entries
			VecDoub taperEventCopy = taperEvent;
			
			realSpectrum(taperNoise, recLen, tempN);
			realSpectrum(taperEventCopy, recLen, tempE);
			
			taperEventCopy = taperEvent;
			cmplxfft(taperEventCopy, recLen, tempEcmplx);
			
			for (int iternf =0; iternf < FreqNyq; iternf++) {
				NoiseSpectrum[itercmp][iternf]+=tempN[iternf]*boostNoise;
				EventSpectrum[itercmp][iternf]+=tempE[iternf];
				KEigenSpectraEvent[itercmp][iterK][iternf] = tempEcmplx[iternf];
				//cout << abs(tempEcmplx[iternf]) << "  " << sqrt(tempE[iternf]) << endl;
				
			}
		}
		//End Kth Eigen estimate ..
		
	}
	
	// End EigenSpectrum estimate for each event and component
	
	//Compute Reciever Functions and their coherence estimates
	for (int itercmps = 0; itercmps < 2; itercmps++) {
		for (int iternf = 0; iternf < FreqNyq; iternf++) {
			Complex Zc(0.0,0.0);
			for (int iterK=0; iterK < K; iterK++) {
				Zc += conj(KEigenSpectraEvent[Z][iterK][iternf]) * 
				KEigenSpectraEvent[itercmps+1][iterK][iternf];
			}
			RF[itercmps][iternf] = Zc / (EventSpectrum[Z][iternf] + NoiseSpectrum[Z][iternf]);
			Coher[itercmps][iternf] = pow(abs(Zc),2.0) / 
			           (EventSpectrum[Z][iternf] * EventSpectrum[itercmps+1][iternf]) ;
			
			//DeltaRF[itercmps][iternf] = 1;

			DeltaRF[itercmps][iternf] = sqrt( (1.0 - Coher[itercmps][iternf]) / 
											 ( (K - 1)*Coher[itercmps][iternf] ) ) * abs(RF[itercmps][iternf]);
			
		}
	}
	
	//TaperRF();
	//VecComplex RadialRFf(RF.ncols(), RF[1]);
	//VecDoub_IO RadialRFt;
	//icmplxfft(RadialRFf, RadialRFt, FreqNyq);
	//printVecDoub(RadialRFt);
	//printVecCmplx(RadialRFf);
	 
	//Build Reciever Function in time domain? Causative & Non Causative Build
	//writeRFdata(10., 30., RadialRFt);
	
}

void MTCDriver::TaperRF(){
	double Hann, ithHann;
	double Pih = M_PI/2.0;
	int nFreqMax = ( FreqMax / FreqRlg) + 1;
	//cout << Pih;
	
	for (int iterF = 0; iterF < FreqNyq; iterF++) {
		if (iterF > FreqMax/FreqRlg + 1) {
			for (int itercmp = 0; itercmp < 2; itercmp++) {
				RF[itercmp][iterF] = 0.0;
				DeltaRF[itercmp][iterF] = 0.0;
			}
		}else {
			Hann = pow( cos( Pih * double(iterF/nFreqMax) ), 2.0);
			//cout << Hann << endl;
			for (int itercmp = 0; itercmp < 2; itercmp++) {
				RF[itercmp][iterF] = Hann * RF[itercmp][iterF];
				//DeltaRF[itercmp][iterF] = Hann * DeltaRF[itercmp][iterF];
			}
			
		}

		
	}
	
}

void MTCDriver::addRF(Mat3DCmplx &cpRF, int pos){
	for (int iterF = 0; iterF < FreqNyq; iterF++) {
		for (int itercmp = 0; itercmp < 2; itercmp++) {
			cpRF[pos][itercmp][iterF] = RF[itercmp][iterF];
		}
	}
}


// actually this saves variance .. just write two more to save spectrum and coherence

void MTCDriver::addCoher(Mat3DDoub &cpDeltaCoher, int pos){
	for (int iterF = 0; iterF < FreqNyq; iterF++) {
		for (int itercmp = 0; itercmp < 2; itercmp++) {
			// Check Coherence before dumping inverting ... 
			 cpDeltaCoher[pos][itercmp][iterF] = DeltaRF[itercmp][iterF];
			//cpDeltaCoher[pos][itercmp][iterF] = Coher[itercmp][iterF];

		}
	}
}


// save coherence here - vertical and horizontal coherence
void MTCDriver::addCorrVH(Mat3DDoub &cpDeltaCoher, int pos){
	for (int iterF = 0; iterF < FreqNyq; iterF++) {
		for (int itercmp = 0; itercmp < 2; itercmp++) {
			// Check Coherence before dumping inverting ...
			//cpDeltaCoher[pos][itercmp][iterF] = DeltaRF[itercmp][iterF];
			cpDeltaCoher[pos][itercmp][iterF] = Coher[itercmp][iterF];
			
		}
	}
}

// save 6 spectrum estimate 3 each for event and for noise
void MTCDriver::addSpecAll(Mat3DDoub &cpSpec, int pos){
	
	
	for (int iterF = 0; iterF < FreqNyq; iterF++) {
		for (int itercmp = 0; itercmp < 3; itercmp++) {
			// Check Coherence before dumping inverting ...
			cpSpec[pos][itercmp][iterF] = EventSpectrum[itercmp][iterF];
			cpSpec[pos][itercmp+3][iterF] = NoiseSpectrum[itercmp][iterF];
		}
		
	}
}




/* Function that saves noise & event spectrum to file. Called by RecFunc.cpp */
void MTCDriver::updateCoher(ofstream &spec, int datType, int irec )
{
	int nFreqMax = (FreqMax / FreqRlg) + 1;
	int FreqSkip = 10;
	
	switch (datType) {
		case VERTSPEC:
			for (int iterF =0 ; iterF < nFreqMax; iterF+=FreqSkip) { 
				spec << iterF*FreqRlg << "\t" <<
				irec << "\t" << Coher[VERTSPEC][iterF] << endl;
			}
			spec << "> " << endl;
			break;
		case RADSPEC:
			for (int iterF =0 ; iterF < nFreqMax; iterF+=FreqSkip) { 
				spec << iterF*FreqRlg << "\t" <<
				irec << "\t" << Coher[RADSPEC][iterF] << endl;
			}
			spec << "> " << endl;
			break;
		case TRANSPEC:
			for (int iterF =0 ; iterF < nFreqMax; iterF+=FreqSkip) { 
				spec << iterF*FreqRlg << "\t" <<
				irec << "\t" << pow(DeltaRF[TRANSPEC][iterF], 2.0) << endl;
			}
			spec << "> " << endl;
			break;
		case NOISESPEC:
			for (int iterF =0 ; iterF < nFreqMax; iterF+=FreqSkip) { 
				spec << iterF*FreqRlg << "\t" <<
				irec << "\t" << pow(DeltaRF[VERTSPEC][iterF], 2.0) << endl;
			}
			spec << "> " << endl;
			break;
		default:
			break;
	}
}

/* Function that saves Coherence & Variance? Called by RecFunc.cpp */
void MTCDriver::updateSpectrum(ofstream &spec, int datType, int irec )
{
	int nFreqMax = (FreqMax / FreqRlg) + 1;
	int FreqSkip = 10;
	
	switch (datType) {
		case VERTSPEC:
			for (int iterF =0 ; iterF < nFreqMax; iterF+=FreqSkip) { 
				spec << iterF*FreqRlg << "\t" <<
				irec << "\t" << log10(EventSpectrum[VERTSPEC][iterF]) << endl;
			}
			spec << "> " << endl;
			break;
		case RADSPEC:
			for (int iterF =0 ; iterF < nFreqMax; iterF+=FreqSkip) { 
				spec << iterF*FreqRlg << "\t" <<
				irec << "\t" << log10(EventSpectrum[RADSPEC][iterF]) << endl;
			}
			spec << "> " << endl;
			break;
		case TRANSPEC:
			for (int iterF =0 ; iterF < nFreqMax; iterF+=FreqSkip) { 
				spec << iterF*FreqRlg << "\t" <<
				irec << "\t" << log10(EventSpectrum[TRANSPEC][iterF]) << endl;
			}
			spec << "> " << endl;
			break;
		case NOISESPEC:
			for (int iterF =0 ; iterF < nFreqMax; iterF+=FreqSkip) { 
				spec << iterF*FreqRlg << "\t" <<
				irec << "\t" << log10(NoiseSpectrum[VERTSPEC][iterF]) << endl;
			}
			spec << "> " << endl;
			break;
		default:
			break;
	}
}



/* ------------Miscelleanuous Routines to print Vectors */
void MTCDriver::printVecDoub(VecDoub_IO &a){
	for (int it =0; it < a.size(); it++) {
		cout  << a[it] << " " << endl;
	}
}
void MTCDriver::printVecCmplx(VecComplex &a){
	for (int it =0; it < a.size(); it++) {
		cout << arg(a[it]) * (180.0/ M_PI) << endl;
	}
}
void MTCDriver::writeRFdata(double pretime, double posttime, VecDoub_IO &RF){
	int npre = pretime / DeltaT, npost = posttime / DeltaT;
	VecDoub time(npre+npost, 0.);
	VecDoub RFdump(npre+npost, 0.);
	
	bool outdata = true;
	
	double scaleRF = 2.0 * FreqNyq / (FreqMax/FreqRlg);
	//double scaleRF = 1.0;
	
	for (int i = 1; i <= (npre+npost) ; i++) {
		time[i-1] =  -npre*DeltaT + (i-1)*DeltaT ;
	}
	
	for (int i =1 ; i <= npost; i++) {
		RFdump[npre+i-1] = scaleRF * RF[i-1];
	}
	
	for (int i = 1; i <= npre; i++) {
		RFdump[npre-i] = scaleRF * RF[recLen - i-1];
	}
	
	if (outdata) {
		for (int i = 0; i < FreqNyq; i++) {
			for (int j = 0; j < 3; j++) {
				cout << EventSpectrum[j][i] << "  " << 
				NoiseSpectrum[j][i] << "   ";
			}
			cout << Coher[0][i] << "  " << Coher[1][i] << endl;
		}
	}else {
		for (int i = 0; i < (npre + npost); i++) {
			cout << time[i] << "  " << RFdump[i] << endl;
		}
	}
	
}

/* End Miscelleanuous Routines */
