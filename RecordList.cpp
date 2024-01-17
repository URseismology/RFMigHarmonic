/**
 *  @file RecordList.cpp
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

#include "RecordList.h"
#include "sacread.h"

using namespace std;

RecordList::RecordList(const char* fn, float time, float Freq){
	nTotrec = 0;
	nGoodrec = 0;
	nNoiseTrace = 0; nEventTrace = 0; nMax = 0;
	Fcutoff = Freq;
	
	/* 
	 Set time length of noise and event here,
	 adaptation should read from console
	 */
	tEventTrace = time; // in sec
	
	ifstream s(fn, std::ios_base::in);
	if (!s) {
		cerr << "RecordList: Cannot open " << fn << "\n";
		exit(1);
	}

	line = "";
	while ( !s.eof() ) {
		s >> line ;
		recordname.push_back(line);
		
	}
	s.close();
	// printVec(recordname); // Diagnostics
	//Check Array list and create a vector for Azimuth & Epicentre
	nTotrec = recordname.size();


	
};	

/** When this is called, It pushes new record into */
void RecordList::passNew(int recPos){
	goodRecordname.push_back(recordname[recPos]);
	nGoodrec++;
	
}


/** Read Header entry for design requirements ... */
void RecordList::removeRecord(int vecIndx){
	
	// Erase Record and Headers....
	goodRecordname.erase(goodRecordname.begin() + vecIndx);
	//Azim.erase(Azim.begin() + vecIndx);		// 1. Azimuth
	//Epic.erase(Epic.begin() + vecIndx);		// 2. Distance
	//EvLong.erase(EvLong.begin() + vecIndx);     // 3. Equake Longitude
	//EvLat.erase(EvLat.begin() + vecIndx);		// 4. Equake Latitude
	//EvMag.erase(EvMag.begin() + vecIndx);		// 5. Equake Magnitude
	
	nGoodrec--;
}
 

/** Important Headers stored for statistics and selective RF computation ****/
void RecordList::setAzim(){
	Azim.assign(nGoodrec,0);
}

void RecordList::setEpic(){
	Epic.assign(nGoodrec,0);
}

void RecordList::setEvLong(){
	EvLong.assign(nGoodrec,0);
}

void RecordList::setEvLat(){
	EvLat.assign(nGoodrec,0);
}
void RecordList::setEvMag(){
	EvMag.assign(nGoodrec,0);
}

void RecordList::initHeaders(){
	Azim.assign(nGoodrec,0);		// 1. Azimuth
	Epic.assign(nGoodrec,0);		// 2. Distance
	EvLong.assign(nGoodrec,0);      // 3. Equake Longitude
	EvLat.assign(nGoodrec,0);		// 4. Equake Latitude
	EvMag.assign(nGoodrec,0);		// 5. Equake Magnitude
	

}


void RecordList::setHeaders(double Az, double Ep, double lng, double lat, double mag,
							int irec){
	Azim[irec] = Az;		// 1. Azimuth
	Epic[irec] = Ep;		// 2. Distance
	EvLong[irec] = lng;      // 3. Equake Longitude
	EvLat[irec] = lat;		// 4. Equake Latitude
	EvMag[irec] = mag;		// 5. Equake Magnitude
	
	// Add headers for timing information??
	
	
}

//Utility for debugig..
void RecordList::printVec(vector<string> &a){
	for (int i = 0; i < a.size(); i++) {
		cout << a[i] << endl;
	}
		
}