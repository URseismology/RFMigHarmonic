/**
 *  @file sacread.cpp
 *
 *
 *  @author Tolulope Olugboji
 *	@date 11/22/10.
 *
 *  @copyright 2010 Yale University. All rights reserved.
 *
 *  @brief File reads SAC triplets using recrord names.
 *	Added LQT rotation on May 23, 2013
 *
 *	@bug no known bugs
 */


#include "sacread.h"
#include <fstream>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
using namespace std;

/* Basic constructor to do just lqt rotation
 */
Sacread::Sacread( string name, int hdrField, float lqtNo ) {
	
	//cout << "File Name is: " << name << endl;
	fName.assign(name);
	
	if (read_sac(name) != -1) {
		// cout << "Sacread: success!" << endl;
		hdrNumber = hdrField;
	
	}else {
		cerr << "Sacread: cannot open all or some of file(s) " << name << endl;
		exit(1);
	}
	getSACfields();
	
	// Check If LQT rotation is to be done - If lqtNo is greater than zero then  rotate	
	if (lqtNo >= 0) {
		isLQTset = 1;
		rotateToLQT(lqtNo);
	}else {
		 cout << "No rotation: LQT argument " << lqtNo << endl;
	}
	

	
};

/* Higher level constructor that does both horizontal and lqt rotation
 */
Sacread::Sacread( string name, int hdrField, float lqtNo, float rotAngle, int sense) {
	
	//cout << "File Name is: " << name << endl;
	fName.assign(name);
	
	if (read_sac(name) != -1) {
		// cout << "Sacread: success!" << endl;
		hdrNumber = hdrField;
		
	}else {
		cerr << "Sacread: cannot open all or some of file(s) " << name << endl;
		exit(1);
	}
	getSACfields();
	
	// Do horizontal rotation before lqtRotation
	if (rotAngle > 0.0) {
		int sense = 1;
		rotateHorizontals(rotAngle, sense);
	}
	

	// Check If LQT rotation is to be done - If lqtNo is greater than zero then  rotate
	if (lqtNo >= 0) {
		isLQTset = 1;
		rotateToLQT(lqtNo);
	}else {
		cout << "No rotation: LQT argument " << lqtNo << endl;
	}
	
	
	
};

/***********************************************************
 Store relevant header variables needed to pick out chunk of
 data from  sac variable.
 
 ** FOR EarthQuake Stats. I include:
		1. EventLong.
		2. EventLat.
		3. Magnitude of Earthquake
  ************************************************************/

void Sacread::getSACfields(){
	float inverseRadiusEarth = 1.0 / (6371.0);
	datLen = hd.npts; 
	deltaT = hd.delta;
	startTime = hd.b;
	
	// Recompute Back Azimuth here using the trig formulae
	recBaz = hd.baz;
	if (1) {
		float lat2 = (hd.evla) * M_PI/180;
        float lat1 = (hd.stla) * M_PI/180;
        float dLong = (hd.evlo - hd.stlo) * M_PI/180;
        
        float bazRad = atan2( (sin(dLong)*cos(lat2) ) , 
							 (cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dLong)) );
        float num = (bazRad * 180/M_PI) + 360.0;
		float den = 360.0;
        recBaz = fmod( num , den );
        
	}
	
	recEpic = hd.gcarc;
	timeStart = hd.b;
	evLon = hd.evlo;
	evLat = hd.evla;
	evMag = hd.mag;
	phaseName.assign(hd.kt0, hd.kt0+8);
	
	
	// Olugboji2016 store event statistics for reporting ..
	if (1) {
		eventStats.append(fName);
		eventStats.append(  "  ");
		eventStats.append( to_string(hd.nzyear));
		eventStats.append(  "  ");
		eventStats.append( to_string(hd.nzjday) );
		eventStats.append(  "  ");
		eventStats.append( to_string(hd.nzhour) );
		eventStats.append(  "-");
		eventStats.append( to_string(hd.nzmin) );
		eventStats.append(  ":");
		eventStats.append( to_string(hd.nzsec) );
		eventStats.append(  ".");
		eventStats.append( to_string(hd.nzmsec) );
		eventStats.append(  "  ");
		eventStats.append( to_string(evLat) );
		eventStats.append(  "  ");
		eventStats.append( to_string(evLon) );
		eventStats.append(  "  ");
		eventStats.append( to_string(hd.evdp) );
		eventStats.append(  "  ");
		eventStats.append( to_string(evMag) );
		eventStats.append(  "  ");
		eventStats.append( to_string(hd.imagtyp) );
		eventStats.append(  "  \n");
	}

	
	// save statistics field in string format ...
	
	if (hd.user0 < 30.0) {
		// unit is in sec/deg so use this conversion
		raySlowness = hd.user0 * (180.0/M_PI)* (inverseRadiusEarth);  //now in sec/km
		//cout << "IN sec/deg. conversion " << raySlowness << endl;
		//raySlowness = hd.user0 * (inverseRadiusEarth);  //now in sec/km
	}else{
		// unit is in sec/rad so use this conversion
		raySlowness = hd.user0 * (inverseRadiusEarth);  //now in sec/km
	}
	
	
	
	timeTag = timeMark_stats();
	
	
	
}

/***********************************************************
 Added on the 23rd of April, 2014
 
 HACK TO DO SIMPLE ROTATIONS OF THE HORIZONTALS AROUND NORTH
 
 Useful for faulty horizontal orientation in ocean bottom data
 to rotate all RT by a specified angle ??
 
 ************************************************************/

void Sacread::rotateHorizontals(float rotAngle, int sense){

	float calcAngle = rotAngle * sense;
	//cout << "Wave incindence angle: " << incAngle << endl;
	
	
	// Now rotate The R and T components.
	// Component[1] - Radial, Component[2] - Transverse
	
	// Transform From ZRT to ZR_newT_new Here.
	double cosInc =  cos(calcAngle * M_PI/180); // Convert Degree to Radians
	double sinInc = sin(calcAngle * M_PI/180); // Convert Degree to Radians
	
	double ithRnval, ithTnval;
	
	for (int iter = 0; iter < hd.npts; iter++) {
		ithRnval = DataCmps[1][iter];
		ithTnval = DataCmps[2][iter];
		DataCmps[1][iter] =	cosInc*ithRnval	+ sinInc*ithTnval;  // Radial now rotated into Radial_new coordinates
		DataCmps[2][iter] =	cosInc*ithTnval	- sinInc*ithRnval;  // Transverse now rotated into Trans_new coordinates
	}
	/* 3 Component Data has been Transformed from ZRT to ZRnewTnew */
	
	
	
	
}

/***********************************************************
 Added on the 23rd of May, 2013
 Use Empirical relationship in Park and Levin, 2013, in review,
 to rotate Vertical (Z) and Radial Components into L(P)-Q(SV) coords.
 using the piecewise incidence angle function.
 
 ** FOR EarthQuake Stats. I include:
 1. EventLong.
 2. EventLat.
 3. Magnitude of Earthquake
 ************************************************************/

void Sacread::rotateToLQT(float lqtNo){
	/* First determine what phase type it is, then select the empirical
	   incidence angle. With incidence angle found, rotate.
	   A future modification of this routine, would perform local 
	   determination of the incidence angle, or would require the incidence
	   angle to be specified in the SAC header.
	   This empirical relationship was determined using station ARU, so
	 mod */
	
	
	bool isP =  (phaseName.compare(0,1,"P") == 0 ? TRUE : FALSE);
	bool isPP = (phaseName.compare(0,2,"PP") == 0 ? TRUE : FALSE);
	bool isPdiff = (phaseName.compare(0,5,"Pdiff") == 0 ? TRUE : FALSE);
	bool isPKP = (phaseName.compare(0,3,"PKP") == 0 ? TRUE : FALSE);
	bool isPKIKP = (phaseName.compare(0,5,"PKIKP") == 0 ? TRUE : FALSE);
	
	bool isPhaseSet = false;
	float ithEpic = recEpic;
	int branch = 0;                         // determine which epicentral branch to use to determine
											// incidence angle.
	
	float incAngle = 9;                     // Incidence angle
	
	if (isPKIKP ) {
		isPhaseSet = TRUE;
	} else if (isPKP && !isPhaseSet) {
		isPhaseSet = TRUE;
	} else if ( isPdiff && !isPhaseSet) {
		isPhaseSet = TRUE;
	} else if ( isPP && !isPhaseSet ) {
		ithEpic = recEpic / 2.0;
		isPhaseSet = TRUE;
	} else if ( isP && !isPhaseSet ) {
		isPhaseSet = TRUE;
	}
	
	
	// Determine the function to use to determine the incidence  angle.
	if (ithEpic <= 80.0) {
		branch = 1;
	} else if ( (ithEpic >= 80.0 && ithEpic <= 118.0) || isPdiff ){
		branch = 2;
	} else if ( ithEpic >= 118.0 || isPKP || isPKIKP) {
		branch = 3;
	}
	
	
	switch (branch) {
		case 1:
			incAngle = 42.0 + (20.0 - recEpic)/3.0 ;
			break;
		case 2:
			incAngle = 16.0 + ( pow((110.0 - recEpic), 2.0)/150.0 ) ;
			break;
		case 3:
			incAngle = 9.0 ;
			break;
		default:
			break;
	}
	
	incAngle = lqtNo * incAngle;
	//cout << "Wave incindence angle: " << incAngle << endl;
	
	
	// Now rotate The Z and R components.
	// Component[0] - Vertical, Component[1] - Radial
	
	// Transform From ZRT to LQT Here.
	double cosInc =  cos(incAngle * M_PI/180); // Convert Degree to Radians
	double sinInc = sin(incAngle * M_PI/180); // Convert Degree to Radians
	
	double ithZval, ithRval;
	
	for (int iter = 0; iter < hd.npts; iter++) {
		ithZval = DataCmps[0][iter]; 
		ithRval = DataCmps[1][iter];
		DataCmps[0][iter] =	cosInc*ithZval	+ sinInc*ithRval;  // Vertical now in L coordinates
		DataCmps[1][iter] =	cosInc*ithRval	- sinInc*ithZval;  // Radial now in Q coordinates
	}
	/* 3 Component Data has been Transformed from ZRT to LQT */
	
	
	
	
}

/***********************************************************
Clean function to return headertime marker with a count of 
 how many time tags were marked
 ************************************************************/

int  Sacread::timeMark_stats(){
	/* Design specification requires using t3. Use just t1 */
	int countTags = 0;
	float timeMark = 0;
	
	switch (hdrNumber) {
		case 0:
			timeMark = hd.t0;
			countTags++;
			break;
		case 1:
			timeMark = hd.t1;
			countTags++;
			break;
		case 2:
			timeMark = hd.t2;
			countTags++;
			break;
		case 3:
			timeMark = hd.t3;
			countTags++;
			break;
		default:
			break;
	}
	
	noTags = countTags;
	return timeMark;
	
}

/***********************************************************
 
 read_sac
 
 Description:	read binary SAC header from file.
 
 Original Author:	Lupei Zhu
 Adapted for c++ functions by: Tolulope Olugboji
 now used to read SAC header as well as Data pts
 It is used to read a 3 component Seismogram
 
 Arguments:	const char *name 	file name
 SACHEAD *hd		SAC header to be filled
 
 Return:	0 if success, -1 if failed
 
 Modify history:
 05/29/97	Lupei Zhu	Initial coding
 11/22/10.  Tolulope Olugboji Adaptations for c++
 ************************************************************/

int Sacread::read_sac(std::string name ) {
	
	int dataLenBytes, startpos, ndata=1;
	Components Component[3] = {Vertical, Radial, Transverse};
	
	for (int cmp = 0; cmp < 3; cmp++) {
		
		// nameNdComp : name appended with component, neccessary for ARRAYS
		string nameNdComp = name;			
		switch (Component[cmp]) {
			case 0:
				nameNdComp.append("Z");
				break;
			case 1:
				nameNdComp.append("R");
				break;
			case 2:
				nameNdComp.append("T");
				break;
			default:
				cout << "Encoding wrong" << endl;
				break;
		}
						
				
		
		ifstream strm( nameNdComp.c_str(), ios::binary);
		
		
		if (!strm) {
					cout <<  "Unable to open " << nameNdComp << endl;
			return -1;
		}
		
		//cout << "Attempting to read ... " << nameNdComp << endl;
		strm.read((char *)(&hd), sizeof(SACHEAD)); // read header here ...
		strm.close();
		//cout << "Read ... " << nameNdComp << endl;
		
		if(sac_byte_order()) {
			swab4((char *)(&hd), SAC_HEADER_SIZE_NUMBERS);
			if(sac_byte_order()) {
				cerr << "sacread: Error determining SAC header" << nameNdComp << endl;
				return(-1);
			}
		}
		ifstream nstrm(nameNdComp.c_str(), ios::binary);
		dataLenBytes = sizeof(float)*(hd.npts) ;
		ndata = (int)hd.npts;
		startpos = sizeof(float)*158;
		nstrm.seekg(startpos);
		//cout << hd.npts << endl;
		if (cmp == 0){
			DataCmps.assign(3, ndata,0);
			DataCmpsRw.assign(3, ndata,0);  // Olugboji2016 added raw, left untouched.
		}
		
		if ( dataLenBytes > 0) {
			 // cout << "Byte size determined!" << endl;
			nstrm.read( (char *)(&Data), dataLenBytes); // read in data here
		}
		nstrm.close();
		
		for (int iter = 0; iter < hd.npts; iter++) {
			DataCmps[Component[cmp]][iter] = Doub (Data.data[iter]);
			DataCmpsRw[Component[cmp]][iter] = Doub (Data.data[iter]); //Rw.
			//cout << Data.data[iter];
		}
		
	}
		return 0;
}


/***********************************************************
 
 sac_byte_order (hd)
 
 Description:  Determine if the header variable (internal4)
 which defines the header version is between 0 and 6
 do_swap  = 1 Try and byte swap
 = 0 Do not byte swap
 Returns:
 0 - No swapping needs to be done
 1 - Swapping is needed
 */
int Sacread::sac_byte_order() {
	if( ( hd.nvhdr > 0  ) && ( hd.nvhdr <= 6 ) ) 
		return(0);
	return(1);
}


/*****************************************************
 
 swab4
 
 Description:	reverse byte order for float/integer
 
 Author:	Lupei Zhu
 
 Arguments:	char *pt	pointer to byte array
 int    n	number of bytes
 
 Return:	none
 
 Modify history:
 12/03/96	Lupei Zhu	Initial coding
 
 ************************************************************/

void Sacread::swab4( char *pt, int n ) {
	int i;
	char temp;
	for(i=0;i<n;i+=4) {
		temp = pt[i+3];
		pt[i+3] = pt[i];
		pt[i] = temp;
		temp = pt[i+2];
		pt[i+2] = pt[i+1];
		pt[i+1] = temp;
	}
}
