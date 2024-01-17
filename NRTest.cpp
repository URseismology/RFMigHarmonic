/*
 *  NRTest.cpp
 *  
 *
 *  Created by Tolulope Olugboji on 11/22/10.
 *  Copyright 2010 Yale University. All rights reserved.
 *
 */


#include "nr3.h"
#include "calendar.h"
#include "moment.h"
#include "fourier.h"

//#include "sacread.h"

//SACHEAD hd;
//char ko[16];

Int main(void) {
	const  Int NTOT=20;
	Int i,jd,nph=2;
	Doub frac,ave,vrnce;
	//VecDoub data(NTOT);
	MatDoub mdata(2,NTOT);
	/*
	for (i=0;i<NTOT;i++){
		flmoon(i,nph,jd,frac);
		data[i]=jd;
		
	}
	avevar(data,ave,vrnce);
	cout << "Average = " << setw(12) << ave;
	cout << "Variance = " << setw(13) << vrnce << endl;
	cout << "The size of sacheader is " << sizeof(SACHEAD) << '\n';
	cout << "The size of int is " << sizeof(int) << '\n';
	cout << "The size of float is " << sizeof(float) << '\n';
	 */
	/*
	enum Components {
		Vertical,
		Radial,
		Transverse
	};
	 Components Component[3] = {Vertical,Radial,Transverse};
	complex Zc(3,4);
	cout << Zc.conj() << endl;
	
	for (int cmp = 0; cmp < 3; cmp++) {
		switch (Component[cmp]) {
			case 0:
				cout << "Vertical" << endl;
				break;
			case 1:
				cout << "Radial" << endl;
				break;
			case 2:
				cout << "Transverse" << endl;
				break;
			default:
				cout << "Encoding wrong" << endl;
				break;
		}
	}
	mdata[1] = data;
	
	*/
	
	//Test the fourier routine, realft.
	// Create test data
	Int len = 16;
	VecDoub data(len);
	for (int i = 0; i < len; i++) {
		data[i] = sin(2*i);
		cout << data[i] << endl;
	}
	realft(data,1);
	cout << "After fft" << endl;
	for (int i = 0; i < len; i++) {
		cout << data[i] << endl;
	}
	realft(data, -1);
	cout << "After fft" << endl;
	for (int i = 0; i < len; i++) {
		//cout << 2/len << endl;
		cout << (2.0/len) * (data[i]) << endl;
	}
						 
	return 0;
}
