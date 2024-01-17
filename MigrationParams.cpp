/**
 * @file MigrationParams.cpp
 *
 *
 *  @date Updated on July 16, 2013
 *  @brief Helper function to calculate delay times for migration. Reads Velocity file.
 *
 *  @author Tolulope  Olugboji on 6/10/13.
 *  @copyright 2013 Yale University. All rights reserved.
 *
 *	Tweaked getTimeDelayPs to bias time plots for MWM. Code also does direct migration for grid search
 *  @bug no known bugs
 */


#include "MigrationParams.h"


using namespace std;

MigrationParams::MigrationParams(const char* velfn){
	std:string line, title;
	int nLayers;
	double iTheta, iPhi;
	double z, vp, vp2, vp4, vs, vs2, rho;
	
	
	if (velfn != NULL) {
		ifstream stream(velfn, std::ios_base::in);
		
		if (!stream) {
			cerr << "MigrationParams: Cannot open ;" << velfn << "\n";
			exit(1);
		}
		
		getline(stream, title);
		stream >> noLayers;
		
		//cout << "Title:  " << title << "No of Layers: "  << noLayers << endl;
		
		
		int iLayer = 1;
		while (!stream.eof()) {
			
			
			stream >> iTheta >> iPhi;
			stream >> z >> vp >> vp2 >> vp4 >> vs >> vs2 >> rho;
			updateLayerParams(z, vp, vp2, vp4, vs, vs2, rho);
			//cout << "Layer " << iLayer << endl;
			//cout << " Azimuth and Tilt: " << iTheta  << " " << iPhi << endl;
			//cout << " Isotropic Vp Vs " << vp  << " " << vs << endl;
			
			if (iLayer >= noLayers+ 1) {
				break;
			}
			iLayer++;
		}
		stream.close();
		isRayGood = true;
	}else {
		cout << "No Velocity File Passed for Migration." << endl;
		noLayers = -1;
	}
	directMig = false;
		

	
	
};

//Grid search constructor - Used by Single Layer Grid Search .
MigrationParams::MigrationParams(double ithH, double ithVp, double ithVs, double ithRho){
	
	noLayers = 1;
	directMig = true;
	updateLayerParams(ithH, ithVp, 0.0, 0.0, ithVs, 0.0, 0.0);
	
}

//Grid search constructor - Used by Multi-Layer Grid Search .
MigrationParams::MigrationParams(Layer ithLayer, vecLayers prevLayers, int nlayersAbove){
	
	if (nlayersAbove == 0) {
		
		noLayers = 1;
		directMig = true;
		updateLayerParams(ithLayer.H, ithLayer.Vp, 0.0, 0.0, ithLayer.Vs, 0.0, 0.0);
	}else{
		
		//cout << "Layers in mig. " << nlayersAbove << endl;
		noLayers = nlayersAbove + 1;
		directMig = true;
		
		for (int iterLayers = 0; iterLayers < nlayersAbove; iterLayers++) {
			// Do I need directMig? Explore ..
			updateLayerParams(prevLayers[iterLayers].H, prevLayers[iterLayers].Vp, 0.0, 0.0, prevLayers[iterLayers].Vs, 0.0, 0.0);
			//cout << "test layers: " << prevLayers[iterLayers].H <<  "  "<< prevLayers[iterLayers].Vp << "  " << prevLayers[iterLayers].Vs << endl;
		}	
		// Current grid parameters ...
		updateLayerParams(ithLayer.H, ithLayer.Vp, 0.0, 0.0, ithLayer.Vs, 0.0, 0.0);
		
		//cout << "test layers: " << ithLayer.H << " " << ithLayer.Vp << "  " <<ithLayer.Vs <<endl;
		
	}

	
}


void MigrationParams::updateLayerParams(double z, double vp, double vp2,
										double vp4, double vs, double vs2, double rho){
	
	//Update orientation of layer... not used. is it? I ignore.
	
	//Update the thickness and depth
	depthStep.push_back(z); density.push_back(rho);
	
	if (depthStep.size() ==  1) {
		thicknessLayer.push_back(depthStep[0]);

	}else{
		int currLayer = depthStep.size()-1;
		int prevLayer = currLayer - 1;
		thicknessLayer.push_back(depthStep[currLayer] - depthStep[prevLayer]);
	}
	
	// Update all velocities ...
	pVelIsotropic.push_back(vp); pVelCos2Theta.push_back(vp2); pVelCos4Theta.push_back(vp4);
	sVelIsotropic.push_back(vs); sVelCos2Theta.push_back(vs2);
	

}

double MigrationParams::getVertTimeDelay(double targetDepthKm, int phaseFlag){

	// Update time delay here ... The time delay is updated incrementally.
	double tDelayP = 0.0; double tDelayS = 0.0;
	for (int iLayer = 0; iLayer < noLayers+1; iLayer++) {
		
		tDelayS += thicknessLayer[iLayer] / sVelIsotropic[iLayer];
		tDelayP += thicknessLayer[iLayer] / pVelIsotropic[iLayer];
	
		double layerDelay = 0.0;
		switch (phaseFlag) {
			case PS:
				layerDelay = tDelayS - tDelayP;
				timeDelayLayer.push_back(layerDelay);		// Ps vertical time delay
				break;
			case PPSMS:
				layerDelay = 2.0 * tDelayS;
				timeDelayLayer.push_back(layerDelay);		// PpSms vertical time delay
				break;
			case PPPMS:
				layerDelay = tDelayS + tDelayP;
				timeDelayLayer.push_back(layerDelay);		// PpPms vertical time delay
				break;
			default:
				break;
		}
		//cout << "!!!!DebugTimeDelay: " << sVelIsotropic[iLayer] << endl;
		
		
	}
	
	double vertTimeDelay = 0.0;								
	
	double tgtDepMetre = targetDepthKm * 1000.0;
	bool inLayer = false;
	
	// search stack of layers and compute vertical time delay through layers.
	// Improvement... add ray parameter correction.. no?
	
	int iLayer = noLayers-1;							// start above halfspace
	
	while ( (inLayer == false) && (iLayer >=  0) ) {
		
		// Put test here to shortcircuit for grd- search... In this case no need to search 
		// for layer depth right?? Indeed - Extend to Sequential Search ...
		
		if (directMig) {
			inLayer = true;
		} else {
			inLayer = (tgtDepMetre > depthStep[iLayer]);
		}
		
		
		if (inLayer || iLayer == 0) {
			
			double extraDepth;
			double extraTimeDelay = 0.0;
			
			if (directMig) {
				extraDepth = 0.0;
			}else{
				extraDepth = tgtDepMetre - depthStep[iLayer];
				
				switch (phaseFlag) {
					case PS:
						extraTimeDelay = extraDepth * 
						(1.0/sVelIsotropic[iLayer+1] - 1.0/pVelIsotropic[iLayer+1]);	// Ps vertical time delay
						break;
					case PPSMS:
						extraTimeDelay = extraDepth * (2.0/sVelIsotropic[iLayer+1]);		// PpSms vertical time delay
						break;
					case PPPMS:
						extraTimeDelay = extraDepth * 
						(1.0/sVelIsotropic[iLayer+1] + 1.0/pVelIsotropic[iLayer+1]);		// PpPms vertical time delay
						break;
					default:
						break;
				}
			}
			
			
			vertTimeDelay = timeDelayLayer[iLayer] + extraTimeDelay;
			
			break;
		}
		
		iLayer--;
	}
	
	if (iLayer == 0) {

		
		switch (phaseFlag) {
			case PS:
				vertTimeDelay = tgtDepMetre * 
				(1.0/sVelIsotropic[iLayer] - 1.0/pVelIsotropic[iLayer]);	// Ps vertical time delay
				break;
			case PPSMS:
				vertTimeDelay = tgtDepMetre * (2.0/sVelIsotropic[iLayer]);		// PpSms vertical time delay
				break;
			case PPPMS:
				vertTimeDelay = tgtDepMetre * 
				(1.0/sVelIsotropic[iLayer] + 1.0/pVelIsotropic[iLayer]);		// PpPms vertical time delay
				break;
			default:
				break;
		}
		
	}
	
	int foundLayer = iLayer;			// C++ style layer index
	int reportLayer = foundLayer + 1;		// Actual visual interpretation of layer
	
	//cout << "Layer for target depth:" << targetDepthKm << " is " << reportLayer	<< " Vertical Time Delay: " << vertTimeDelay << endl;
	return vertTimeDelay;  // no. return slowness dependent timedelay if rayParam = 0.
}

double MigrationParams::getTimeDelayPs(double targetDepthKm, double rayParam){
	
	// Update time delay here ... The time delay is updated incrementally.
	rayParam = rayParam / 1000.0;			// now in sec per metre.
	
	double timeDelay = 0.0;			// time delay through bottom of stack
	double tgtDepMetre = targetDepthKm * 1000.0;
	bool inLayer = false;
	double extraDepth;
	
	// search stack of layers and compute vertical time delay through layers.
	// Improvement... add ray parameter correction.. no?
	
	int iLayer = noLayers-1;				// start above halfspace
	while ( (inLayer == false) && (iLayer >=  0) ) {
		inLayer = (tgtDepMetre >= depthStep[iLayer]);
		
		if (inLayer || iLayer == 0) {
			
			if (directMig) {
				extraDepth = 0.0;
				break;
			}else{
				extraDepth = tgtDepMetre - depthStep[iLayer];
				break;
			}
		}
		
		iLayer--;
	}
	
	
	int foundLayer = iLayer;			// C++ style layer index
	int reportLayer = foundLayer + 1;		// Actual visual interpretation of layer
	
	// Check if ray impinges here ...
	isRayGood = ( (rayParam * pVelIsotropic[foundLayer]) < 1.0 );
	if (!isRayGood) {
		//cout << "Layer for target depth:" << targetDepthKm << " is "<< reportLayer << " Ray Does Not Impinge!! " << isRayGood << endl;
		return -1.0;			// if evanescent: fail & return negative time shift
	}
	
	
	// Now run ray through entire stack and calculate time Delay (orTimeShift)
	double sVelCosine, pVelCosine, sVelRayP, pVelRayP;
	if (foundLayer == 0) {
		// Within First Layer ...
		sVelRayP = sVelIsotropic[foundLayer] * rayParam;
		pVelRayP = pVelIsotropic[foundLayer] * rayParam;
		
		sVelCosine = sqrt(1.0 - pow(sVelRayP, 2.0));
		pVelCosine = sqrt(1.0 - pow(pVelRayP, 2.0));
		
		timeDelay += tgtDepMetre *
		(sVelCosine/sVelIsotropic[foundLayer] - pVelCosine/pVelIsotropic[foundLayer]);
		
	}else{
		// accumulate time delay above the foundLayerth Layer ...
		for (int currLayer = 0; currLayer <= foundLayer; currLayer++) {
			sVelRayP = sVelIsotropic[foundLayer] * rayParam;
			pVelRayP = pVelIsotropic[foundLayer] * rayParam;
			
			sVelCosine = sqrt(1.0 - pow(sVelRayP, 2.0));
			pVelCosine = sqrt(1.0 - pow(pVelRayP, 2.0));
			
			timeDelay += thicknessLayer[currLayer] *
			(sVelCosine/sVelIsotropic[currLayer] - pVelCosine/pVelIsotropic[currLayer]);
		}
		
		// use depth below the stack to accumulate time delay at targetDepth
		// Don't look below if directMig! for grid stack ...
		if (directMig) {
			
		}else{
			sVelRayP = sVelIsotropic[foundLayer+1] * rayParam;
			pVelRayP = pVelIsotropic[foundLayer+1] * rayParam;
			
			sVelCosine = sqrt(1.0 - pow(sVelRayP, 2.0));
			pVelCosine = sqrt(1.0 - pow(pVelRayP, 2.0));
			
			timeDelay += extraDepth *
			(sVelCosine/sVelIsotropic[foundLayer+1] - pVelCosine/pVelIsotropic[foundLayer+1]);
		}

	}

	

	
	// cout << "Layer for target depth:" << targetDepthKm << "(km) is " << reportLayer << " Ray Impinges? " << isRayGood << "Timing: " << timeDelay << " ExtraDepth: " << extraDepth <<  " Depth: " << depthStep[iLayer] << " Mig? " << directMig <<  endl;
	return timeDelay;  // no. return slowness dependent timedelay if rayParam = 0.
}

// Return Depth, when given time
double MigrationParams::getVertDepthPs(double PsTime, int phaseFlag){

	
	
	accumTimeDelayLayer.assign(noLayers+1, 0.0);
	timeDelayLayer.assign(noLayers+1, 0.0);
	// Update time delay for EACH LAYER! ... INITIALIZE? Check Make sure memory leak not happening here.
	double tDelayP = 0.0; double tDelayS = 0.0;
	for (int iLayer = 0; iLayer < noLayers+1; iLayer++) {
		
		tDelayS += thicknessLayer[iLayer] / sVelIsotropic[iLayer];
		tDelayP += thicknessLayer[iLayer] / pVelIsotropic[iLayer];
		
		double layerDelay = 0.0;
		switch (phaseFlag) {
			case PS:
				layerDelay = tDelayS - tDelayP;
				timeDelayLayer[iLayer] = (layerDelay);		// Ps vertical time delay
				break;
			case PPSMS:
				layerDelay = 2.0 * tDelayS;
				timeDelayLayer[iLayer] = (layerDelay);		// PpSms vertical time delay
				break;
			case PPPMS:
				layerDelay = tDelayS + tDelayP;
				timeDelayLayer[iLayer] = (layerDelay);		// PpPms vertical time delay
				break;
			default:
				break;
		}
		//cout << "!!!!DebugTimeDelay: " << sVelIsotropic[iLayer] << endl;
		
		
	}
	
	// Accumulate total time delay INCREMENTALLY for search - WRONG! I THINK??
	/*
	for (int iLayer = 0; iLayer < noLayers+1; iLayer++) {
		
		int currLayer = iLayer;
		int topLayer = iLayer;
		
		if (currLayer > 0) {
			
			while (topLayer >= 0) {
				accumTimeDelayLayer[iLayer] += timeDelayLayer[topLayer];
				topLayer--;
			}
			
		}else{
			// This is the first layer so just push on top.
			accumTimeDelayLayer[iLayer] = timeDelayLayer[iLayer];
			
		}
		
		cout << "\n Acum. TDelay: " << timeDelayLayer[iLayer] << endl;

		
		
		
	}
	*/
	
	// Algorithm to find layer using timing information
	// search from top ...
	
	int iLayer = 0; // start above halfspace
	bool inLayer = false;
	
	while ( (inLayer == false) && (iLayer <  noLayers + 1) ) {
		inLayer = (PsTime <= timeDelayLayer[iLayer]);
		
		if (inLayer || (iLayer == noLayers) ) {
			break;
		}else{
			iLayer++;
		}
		
		
	}
	
	
	// Layer already found so just calculate depth iteratively one after the other till you hit the actual layer ..
	double migDepthInMetre = 0.0;		// Depth for Ps time delay
	double thicknessTopLayers = 0.0;
	double extraTimeDelay = 0.0;
	
		
	if (iLayer == 0) {
		thicknessTopLayers = 0.0;
		extraTimeDelay = PsTime;
	}else if ( iLayer > (noLayers+1) ){
		thicknessTopLayers = depthStep[noLayers];
		extraTimeDelay = PsTime - timeDelayLayer[noLayers];
	}else if ( (iLayer > 0) && (iLayer < (noLayers + 1) ) ){
		thicknessTopLayers = depthStep[iLayer-1];
		extraTimeDelay = PsTime - timeDelayLayer[iLayer-1];
	}

	switch (phaseFlag) {
		case PS:
			migDepthInMetre = extraTimeDelay /
			(1.0/sVelIsotropic[iLayer] - 1.0/pVelIsotropic[iLayer]);	// Ps vertical time delay
			break;
		case PPSMS:
			migDepthInMetre = extraTimeDelay / (2.0/sVelIsotropic[iLayer]);		// PpSms vertical time delay
			break;
		case PPPMS:
			migDepthInMetre = extraTimeDelay /
			(1.0/sVelIsotropic[iLayer] + 1.0/pVelIsotropic[iLayer]);		// PpPms vertical time delay
			break;
		default:
			break;
	}
	
		
	int foundLayer = iLayer;			// C++ style layer index
	int reportLayer = foundLayer + 1;		// Actual visual interpretation of layer
	//cout << "Migrated depth is in Layer: " << reportLayer << " of " << noLayers + 1 << endl;
	//cout << "MigDepth: " << migDepthInMetre /1000.0 << "Top Layer: " << thicknessTopLayers/1000 << endl;
	return (migDepthInMetre + thicknessTopLayers);
			

	
	

	
	
}

double MigrationParams::getBotLayer(double targetDepthKm, double rayParam){
	
	// Update time delay here ... The time delay is updated incrementally.
	double tDelayP = 0.0; double tDelayS = 0.0;
	for (int iLayer = 0; iLayer < noLayers+1; iLayer++) {
		
		tDelayS += thicknessLayer[iLayer] / sVelIsotropic[iLayer];
		tDelayP += thicknessLayer[iLayer] / pVelIsotropic[iLayer];
		
		double layerDelay = tDelayS - tDelayP;
		timeDelayLayer.push_back(layerDelay);
	}
	
	double vertTimeDelay = 0.0;
	double tgtDepMetre = targetDepthKm * 1000.0;
	bool inLayer = false;
	
	// search stack of layers and compute vertical time delay through layers.
	// Improvement... add ray parameter correction.. no?
	
	int iLayer = noLayers;							// start at halfspace
	
	while ( (inLayer == false) && (iLayer >=  0) ) {
		inLayer = (tgtDepMetre > depthStep[iLayer]);
		
		if (inLayer) {
			
			double extraDepth = tgtDepMetre - depthStep[iLayer];
			double extraTimeDelay = extraDepth *
			(1.0/sVelIsotropic[iLayer+1] - 1.0/pVelIsotropic[iLayer+1]);
			vertTimeDelay = timeDelayLayer[iLayer] + extraTimeDelay;
			
			break;
		}
		
		iLayer--;
	}
	
	if (iLayer < 0) {
		vertTimeDelay = tgtDepMetre *
		(1.0/sVelIsotropic[iLayer+1] - 1.0/pVelIsotropic[iLayer+1]);
		
	}
	
	int foundLayer = iLayer + 1;			// C++ style layer index
	int reportLayer = foundLayer + 1;		// Actual visual interpretation of layer
	
	//cout << "Layer for target depth:" << targetDepthKm << " is " << reportLayer << endl;
	return vertTimeDelay;  // no. return slowness dependent timedelay if rayParam = 0.
}

double MigrationParams::getTimeDelayPpSms(double targetDepthKm, double rayParam){
	
	// Update time delay here ... The time delay is updated incrementally.
	rayParam = rayParam / 1000.0;			// now in sec per metre.
	
	double timeDelay = 0.0;			// time delay through bottom of stack
	double tgtDepMetre = targetDepthKm * 1000.0;
	bool inLayer = false;
	double extraDepth;
	
	// search stack of layers and compute vertical time delay through layers.
	// Improvement... add ray parameter correction.. no?
	
	int iLayer = noLayers-1;				// start above halfspace
	while ( (inLayer == false) && (iLayer >=  0) ) {
		inLayer = (tgtDepMetre >= depthStep[iLayer]);
		
		if (inLayer || iLayer == 0) {
			extraDepth = tgtDepMetre - depthStep[iLayer];
			break;
		}
		
		iLayer--;
	}
	
	
	int foundLayer = iLayer;			// C++ style layer index
	int reportLayer = foundLayer + 1;		// Actual visual interpretation of layer
	
	// Check if ray impinges here ...
	isRayGood = ( (rayParam * pVelIsotropic[foundLayer]) < 1.0 );
	if (!isRayGood) {
		//cout << "Layer for target depth:" << targetDepthKm << " is "<< reportLayer << " Ray Does Not Impinge!! " << isRayGood << endl;
		return -1.0;			// if evanescent: fail & return negative time shift
	}
	
	
	// Now run ray through entire stack and calculate time Delay (orTimeShift)
	double sVelCosine, pVelCosine, sVelRayP, pVelRayP;
	if (foundLayer == 0) {
		// Within First Layer ...
		sVelRayP = sVelIsotropic[foundLayer] * rayParam;
		pVelRayP = pVelIsotropic[foundLayer] * rayParam;
		
		sVelCosine = sqrt(1.0 - pow(sVelRayP, 2.0));
		pVelCosine = sqrt(1.0 - pow(pVelRayP, 2.0));
		
		// Two S and 1 P phase = [P + 2S] - P = 2S... Flat layer PpSms
		timeDelay += tgtDepMetre * (2.0 * sVelCosine/sVelIsotropic[foundLayer]);
		
	}else{
		// accumulate time delay above the foundLayerth Layer ...
		for (int currLayer = 0; currLayer <= foundLayer; currLayer++) {
			sVelRayP = sVelIsotropic[foundLayer] * rayParam;
			pVelRayP = pVelIsotropic[foundLayer] * rayParam;
			
			sVelCosine = sqrt(1.0 - pow(sVelRayP, 2.0));
			pVelCosine = sqrt(1.0 - pow(pVelRayP, 2.0));
			
			// Two S and 1 P phase = [P + 2S] - P = 2S... Flat layer PpSms
			timeDelay += thicknessLayer[currLayer] * 
					(2.0 * sVelCosine/sVelIsotropic[foundLayer]);
		}
		
		// use depth below the stack to accumulate time delay at targetDepth
		// Don't look below if directMig! Grid Stack has no halfspace!
		
		if (directMig) {
			
		}else{
			sVelRayP = sVelIsotropic[foundLayer+1] * rayParam;
			pVelRayP = pVelIsotropic[foundLayer+1] * rayParam;
			
			sVelCosine = sqrt(1.0 - pow(sVelRayP, 2.0));
			pVelCosine = sqrt(1.0 - pow(pVelRayP, 2.0));
			
			// Two S and 1 P phase = [P + 2S] - P = 2S... Flat layer PpSms
			timeDelay += extraDepth *(2.0 * sVelCosine/sVelIsotropic[foundLayer]);
		}

	}
	
	
	
	
	// cout << "Layer for target depth:" << targetDepthKm << "(km) is " << reportLayer << " Ray Impinges? " << isRayGood << endl;
	return timeDelay;  // no. return slowness dependent timedelay if rayParam = 0.
}

double MigrationParams::getTimeDelayPpPms(double targetDepthKm, double rayParam){
	
	// Update time delay here ... The time delay is updated incrementally.
	rayParam = rayParam / 1000.0;			// now in sec per metre.
	
	double timeDelay = 0.0;			// time delay through bottom of stack
	double tgtDepMetre = targetDepthKm * 1000.0;
	bool inLayer = false;
	double extraDepth;
	
	// search stack of layers and compute vertical time delay through layers.
	// Improvement... add ray parameter correction.. no?
	
	int iLayer = noLayers-1;				// start above halfspace
	while ( (inLayer == false) && (iLayer >=  0) ) {
		inLayer = (tgtDepMetre >= depthStep[iLayer]);
		
		if (inLayer || iLayer == 0) {
			extraDepth = tgtDepMetre - depthStep[iLayer];
			break;
		}
		
		iLayer--;
	}
	
	
	int foundLayer = iLayer;			// C++ style layer index
	int reportLayer = foundLayer + 1;		// Actual visual interpretation of layer
	
	// Check if ray impinges here ...
	isRayGood = ( (rayParam * pVelIsotropic[foundLayer]) < 1.0 );
	if (!isRayGood) {
		
		// cout << "Layer for target depth:" << targetDepthKm << " is " << reportLayer << " Ray Does Not Impinge!! " << isRayGood << endl;
		return -1.0;			// if evanescent: fail & return negative time shift
	}
	
	
	// Now run ray through entire stack and calculate time Delay (orTimeShift)
	double sVelCosine, pVelCosine, sVelRayP, pVelRayP;
	if (foundLayer == 0) {
		// Within First Layer ...
		sVelRayP = sVelIsotropic[foundLayer] * rayParam;
		pVelRayP = pVelIsotropic[foundLayer] * rayParam;
		
		sVelCosine = sqrt(1.0 - pow(sVelRayP, 2.0));
		pVelCosine = sqrt(1.0 - pow(pVelRayP, 2.0));
		
		// Two P and 1 S phase = [2P + S] - P = S + P
		timeDelay += tgtDepMetre *
		(sVelCosine/sVelIsotropic[foundLayer] + pVelCosine/pVelIsotropic[foundLayer]);
		
	}else{
		// accumulate time delay above the foundLayerth Layer ...
		for (int currLayer = 0; currLayer <= foundLayer; currLayer++) {
			sVelRayP = sVelIsotropic[foundLayer] * rayParam;
			pVelRayP = pVelIsotropic[foundLayer] * rayParam;
			
			sVelCosine = sqrt(1.0 - pow(sVelRayP, 2.0));
			pVelCosine = sqrt(1.0 - pow(pVelRayP, 2.0));
			
			// Two P and 1 S phase = [2P + S] - P = S + P
			timeDelay += thicknessLayer[currLayer] *
			(sVelCosine/sVelIsotropic[currLayer] + pVelCosine/pVelIsotropic[currLayer]);
		}
		
		// use depth below the stack to accumulate time delay at targetDepth
		// Don't look below if directMig! grid stack has no halfspace
		
		if (directMig) {
		}else{
			sVelRayP = sVelIsotropic[foundLayer+1] * rayParam;
			pVelRayP = pVelIsotropic[foundLayer+1] * rayParam;
			
			sVelCosine = sqrt(1.0 - pow(sVelRayP, 2.0));
			pVelCosine = sqrt(1.0 - pow(pVelRayP, 2.0));
			
			// Two P and 1 S phase = [2P + S] - P = S + P
			timeDelay += extraDepth *
			(sVelCosine/sVelIsotropic[foundLayer+1] + pVelCosine/pVelIsotropic[foundLayer+1]);
			
		}

	}
	
	
	
	
	//cout << "Layer for target depth:" << targetDepthKm << "(km) is "<< reportLayer << " Ray Impinges? " << isRayGood << endl;
	return timeDelay;  // no. return slowness dependent timedelay if rayParam = 0.
}

int MigrationParams::getNoLayers(){
	return noLayers;
}

