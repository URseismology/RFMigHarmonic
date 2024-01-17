/**
 * @file MigrationParams.h
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "nr3.h"


/** Phase Flag for Direct PS phase */
#define PS 0  
/** Phase Flag for Converted Phase PPSMS phase */
#define PPSMS 1 
/** Phase Flag for Converted Phase PPPMS phase */
#define PPPMS 2 
/** Phase Flag for ALL phases: used in H-K stack for full summary stacks. */
#define FULL 3

/** think of the grid parameters as a single node point of triplets - h, vp, & vs */
struct layerParams{
	double H;		//Thickness
	double Vp;		// p velocity
	double Vs;		// s velocity
};

typedef vector<layerParams> vecLayers;
typedef layerParams  Layer ;

class MigrationParams{
public:
	
	/**
	 * Constructor builds migration parameters using velocify file.
	 */
	MigrationParams(const char* velfn);
	
	/**
	 * Constructor builds migration parameters using ith Grid parameters - single layer
	 */
	MigrationParams(double ithH, double ithVp, double ithVs, double ithRho); 
	
	/**
	 * Constructor builds migration parameters using ith Grid parameters - multiple layers
	 */
	MigrationParams(Layer ithLayer, vecLayers prevLayers, int nlayersAbove); 
	
	/**
	 * no of layers used in migration
	 */
	int noLayers;
	
	/**
	 * time delay for each layer
	 */
	vector<double> timeDelayLayer;
	
	/**
	 * total time delay for each layer, accumulated
	 * from the top, used for PsTime search
	 * and depth migration ...
	 */
	vector<double> accumTimeDelayLayer;
	
	/**
	 * Calculate time delay for migration. Add a phase flag for reverberated phases.
	 * @param targetDepthKm is the migration depth in km
	 * @param phaseFlag is the code for phase @see PS
	 * @return vertical time delay. 
	 */
	double getVertTimeDelay(double targetDepthKm, int phaseFlag);   
			
	
	/**
	 * Calculate time delay for migration for PS phase @see PS  phaseFlag
	 * @param targetDepthKm is the migration depth in km
	 * @param rayParam is the ray parameter in sec/km
	 * @return time delay for earthquake with ray parameter
	 */
	double getTimeDelayPs(double targetDepthKm, double rayParam);
	
	/**
	 * Does reverse calculation for depth @see PS  phaseFlag
	 * @param PS time of ray through vel model. Must be +ve
	 * @return depth value for the particular PS time.
	 */
	double getVertDepthPs(double PsTime, int phaseFlag);
	
	/**
	 * Constructor builds migration parameters using velocify file.
	 */
	// retrieve timeDelay
	double getBotLayer(double targetDepthKm, double rayParam); // retrieve bottom of stack
	
	/**
	 * @return no. of layers in migration model.
	 */
	int getNoLayers();

	/**
	 * Calculate time delay for migration for PpSmS phase @see PPSMS  phaseFlag
	 * @param targetDepthKm is the migration depth in km
	 * @param rayParam is the ray parameter in sec/km
	 * @return time delay for earthquake with ray parameter
	 */
	double getTimeDelayPpSms(double mohoDepth, double rayParam);
	
	/**
	 * Calculate time delay for migration for PpPmS phase @see PPPMS  phaseFlag
	 * @param targetDepthKm is the migration depth in km
	 * @param rayParam is the ray parameter in sec/km
	 * @return time delay for earthquake with ray parameter
	 */
	double getTimeDelayPpPms(double mohoDepth, double rayParam);
	

	
	
private:
	VecDoub theta, phi;				// Azimuth, and Dip
	vector<double> depthStep, thicknessLayer, density;
	bool isRayGood;	 // set to false if ray does not impinge bottom of velocity stack...
	bool directMig;	// set to true if used by Grd search for single layer migration ..
	
	// Depth, Layer thickness and Density
	vector<double> pVelIsotropic, pVelCos2Theta, pVelCos4Theta; 
	// P Velocity with anisotropic perturbation
	vector<double> sVelIsotropic, sVelCos2Theta;
	// S Velocity with anisotropic perturbation
	
	void updateLayerParams(double z, double vp, double vp2,
						   double vp4, double vs, double vs2, double rho);
	
	
};

