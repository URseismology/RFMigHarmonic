/**
 *  @file RecordList.h
 *
 *
 *  @author Tolulope Olugboji 
 *	@date   12/12/10.
 *  @copyright Yale University. All rights reserved.
 *
 *  @brief This class builds the list of records used to load the seismograms. Using the input file option. 
 *	
 *	 Currently this class expect single station records, maintains list of 
 *	 header information and records that pass the QC check.
 *	 Future versions(extensions) of this class should be able to load seismograms
 *   and store them in data structures based on their locations in a grid. a line. or by
 *	 common conversion points.
 *
 * @bug no known bugs
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "../nr_c304_bl/code/nr3.h"

class RecordList{
public:
	RecordList(const char*, float, float); //file input
	void BuildRecords();
	void printVec(vector<string> &a);
	void passNew(int);
	void setAzim(); void setEpic(); void setEvLat(); void setEvLong(); void setEvMag();
	
	void initHeaders();
	void setHeaders(double Az, double Ep, double lng, double lat, double mag, int irec);
	
	// Added June 26 - 2013 --- Method to remove defective
	// records. i.e. earthquakes that can't be migrated
	// because they are trapped in layer ...
	void removeRecord(int);
	
	
	// Data for tracking records
	int nTotrec;
	int nGoodrec;
	
	// I need to extend this for variable timing information - e.g. BOREHOLE KAWAKATSU.
	float tEventTrace, tNoiseTrace, Fcutoff;
	int nNoiseTrace, nEventTrace, nMax;
	
	vector<string> recordname; 
	vector<string> goodRecordname;              //All records that pass tag Test
	
	std::string line;
	VecDoub Azim, Epic, EvLong, EvLat, EvMag;
	Mat3DDoub noiseTrace, postEventTrace, CoherTrace;
	Mat3DCmplx RFTrace;
	std::string temp;
};
