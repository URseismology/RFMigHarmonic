/**
 *  @file TraceStack.h
 *  
 *
 *  @author Tolulope Olugboji 
 *	@date 12/26/10.
 *  @copyright Yale University. All rights reserved.
 *
 *  @brief This class is a specialist class that stacks the Reciever functions using the coherence estimates. 
 *
 *  TraceStack does azimuthal stacks, epicentral stacks.
 *	It has an operator that allows the printing
 *  of the final stack as files or stdout
 *	
 *	@warning Currently designing a method to do summary stacks for H-K grid stacks 
 *  
 *  Updated on July 13, 2013 --- 
 *  Updated print function to bias time axis with time delay for moving window migration display. Not useful
 * @bug no known bugs
 */
#include "nr3.h"
#include <sstream>

/** Azimuth Stack */
#define  STACKAZIM 0
/** Epicentral Stack */
#define  STACKEPIC 1
/** Grid Stack */
#define  STACKGRD  2
/** Epicentral Stack With Jacknife*/
#define  STACKEPICJCK 3


/** 
 * Stack Receiver Functions in Frequency Domain. Types of stack:
 * 1. Stack By Azimuth
 * 2. Stack By Epicentral Distance
 * 3. Grid Stack - Single Summary Stack for Set of Earthquakes.
 */
class TraceStack{
private:
	VecDoub Azim, Epic;
	VecDoub BinCntStat, BinCntrStat;    // statistics for  no. of records stacked in bin (epicentral or azimuth)
	VecDoub recIndxInBin;
	Mat3DCmplx RFTrace, stackRFTrace; 
	Mat3DDoub  sigmaRFTrace;
	MatDoub azimStack, epicStack, radialRFAscii, transRFAscii, timeAscii;
	MatDoub iJckRadialRFAscii, avgRadialRFAscii, devRadialRFAscii;
	MatDoub iJckTransRFAscii, avgTransRFAscii, devTransRFAscii; // Records are for jacknife matrices for each stack and for each record ...
	MatDoub depthAscii; bool flagDepthMig; // This is value time transformed to depth using velocity migration.
	bool setPrintStats; // set variable if code should print stats ....
	float binmin, binmax, bininc;
	float FreqNyq, FreqRlg, FreqMax, DeltaT;
	int nFreqMax;
	float preTime, postTime;
	int outLen, nPad, expo, iZero, padNyq, binSze;
	Complex zero;
	double Hann;
	double Pih;
	

	
public:
	/**
	 * @brief  TraceStack Constructor - Needed by operations
	 * @param RF holds radial and transverse receiver func for all records in freq. domain
	 * @param sigma is the estimate of coherence for all records
	 * @param Azim is the back-Azimuth for each record
	 * @param Epic is the epicentral for each record
	 * @param fMax is the cut off Frequency 
	 * @param delT is the inverse sampling rate
	 * @param pre is the total time before zeroTime for RF display in the time domain: ASCII print
	 * @param post is the total time after zeroTime for RF display in the time domain: ASCII print
	 * @param bSze epicentral or Azimuth decimation: set to single for summaryStacks.
	 * @bug No known bugs.
	 */
	TraceStack(Mat3DCmplx &RF, Mat3DDoub &sigma, VecDoub &Azim, VecDoub &Epic,
			   float fMax, float delT, float pre, float post, int bSze);
	
	/**
	 * @brief  Stack RF traces in freq. domain - ordering by Back-Azimuth
	 * @param max is the maximum back azimuth for arriving earthquake
	 * @param min is the minimum back azimuth for arriving earthquake
	 * @param inc is decimation in the back-Azimuth range: set by user.
	 * @param RF is the radial and transverse receiver func in freq. domain
	 * @param sigma is the estimate of coherence in the freq. domain
	 * @return no return. updates member stackRFTrace.
	 * @bug No known bugs.
	 */
	void stackAzim(float max, float min, float inc, 
				   Mat3DCmplx &RF, Mat3DDoub &sigma);
	
	/**
	 * @brief  Stack RF traces in freq. domain - ordering by Epicentral Distance
	 * @param epbazmax is the maximum back azimuth for arriving earthquake
	 * @param epbazmin is the minimum back azimuth for arriving earthquake
	 * @param max is the maximum distance for arriving earthquake
	 * @param min is the minimum distance for arriving earthquake
	 * @param inc is decimation in the distance range: set by user.
	 * @param RF is the radial and transverse receiver func in freq. domain
	 * @param sigma is the estimate of coherence in the freq. domain
	 * @return no return. updates member stackRFTrace.
	 * @bug No known bugs.
	 */
	void stackEpic(float epmax, float epmin, float max, float min, float inc,
				   Mat3DCmplx &RF, Mat3DDoub &sigma);
	
	
	/** 
	 * @brief  jacknife stacking of RF traces in freq. domain - ordering by Epicentral Distance
	 * @param epbazmax is the maximum back azimuth for arriving earthquake
	 * @param epbazmin is the minimum back azimuth for arriving earthquake
	 * @param max is the maximum distance for arriving earthquake
	 * @param min is the minimum distance for arriving earthquake
	 * @param inc is decimation in the distance range: set by user.
	 * @param RF is the radial and transverse receiver func in freq. domain
	 * @param sigma is the estimate of coherence in the freq. domain
	 * @return no return. updates member stackRFTrace.
	 * @bug No known bugs.
	 */
	void jacknifeStackEpic(float epmax, float epmin, float max, float min, float inc,
				   Mat3DCmplx &RF, Mat3DDoub &sigma);
	
	
	/**
	 * @brief Returns Zero time result for H-K Stack GrdSearch ... [RFs are parsed in sectors]
	 * @param epbazmax is the maximum back azimuth for arriving earthquake
	 * @param epbazmin is the minimum back azimuth for arriving earthquake
	 * @param max is the maximum distance for arriving earthquake
	 * @param min is the minimum distance for arriving earthquake
	 * @param inc is decimation in the distance range (not used):  set to max-min+1
	 * @param RF is the radial and transverse receiver func in freq. domain
	 * @param sigma is the estimate of coherence in the freq. domain
	 * @return return maximum value in stack at time zero...
	 * @warning Currently working on this. Notify when closed.
	*/
	double grdStack(float epbazmax, float epbazmin, float max, float min, float inc,
				   Mat3DCmplx &RF, Mat3DDoub &sigma);
	
	/**
	 * @brief Finds the maximum amplitude for Radial RF(time) at zero time
	 * @param RFout radial receiver function in the time domain
	 * @param timeOut time vector showing the timing of RF
	 * @return maxOut is the maximum distance for arriving earthquake
	 * @return minInd is the minimum distance for arriving earthquake
	 * @warning Currently working on this. future tweak?
	 */
	void findMAX(MatDoub &RFout, MatDoub &timeOut, double& maxOut, double& maxInd);
	
	/**
	 * @brief Builds RF in time domain for ASCII display using pre and post time information
	 * @param prelen length of time for pre-time informatin for RF display
	 * @param postlen length of time for pre-time informatin for RF display
	 * @param RFin radial or transverse receiver function in the time domain
	 * @param RFout the truncated receiver function in the time domain ready for printing
	 * @param timeOut the timing information for printing
	 * @param ithStck the stack index - Azimuth, Epicentral Distance or SingleStack
	 * @return maxOut is the maximum distance for arriving earthquake
	 * @return minInd is the minimum distance for arriving earthquake
	 * @warning Currently working on this. Notify when closed.
	 */
	void buildAsciiOut(double prelen, double postlen, VecDoub_IO &RFin,
					   MatDoub &RFout, MatDoub &timeOut, int ithStck);
	
	/**
	 * @brief  prints RF traces to file using macros to guide display
	 * @param fName is the file name for output. function appends extra text depending on stack type
	 * @param stckFlag is integer code for printing uses MACROS STACKAZIM STACKEPIC and STACKGRD
	 * @param timeShift is the time shift bias. Set to Zero by default.
	 * @bug No known bugs.
	 */
	void print(char* fName, int stckFlag, int timeShift); // resolve stream type. Filestream? FileStream Always
	
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
	string logReport;

	
	template <class T>
	inline std::string to_string (const T& t)
	{
		std::stringstream ss;
		ss << t;
		return ss.str();
	}
};