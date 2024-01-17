/**
 *  @file sacread.h
 *  
 *
 *  @author Tolulope Olugboji 
 *	@date 11/22/10.
 *
 *  @copyright 2010 Yale University. All rights reserved.
 *  
 *  @brief File reads SAC triplets [R, T, Z] using recrord names. 
 *	Added LQT rotation on May 23, 2013
 *  Added Horizontal rotation on April 23, 2014 
 *				  -
 *				  -
 *                - Useful for faulty horizontals in ocean bottom (OBS) data
 *
 *	@bug no known bugs
 */


#include "sachead.h"
#include "nr3.h"
#include <string>


class Sacread  {
public:
	Sacread(std::string, int hdr, float lqtNo);
	
	/* Add a new constructor that conducts horizontal rotation before lqt rotation */
	Sacread(std::string, int hdr, float lqtNo, float rotAngle, int sense);

		
	/* get function to collect relevant data from header */
	void  getSACfields();
	
	/* Rotate Records into the LQT coordinates if flags sete */
	void rotateToLQT(float lqtNo);
	
	/* Rotate Horizontals - Hook this up with new constructor for backward compatibility */
	// In constructor, indicate a new argument - 1 extra for hor.angle, and sense maybe.
	void rotateHorizontals(float rotAngle, int sense);
	
	
	// SAC Data retrieved from file and accessible through header
	SACHEAD hd; float deltaT, preData, postData, startTime;
	int noTags, datLen, hdrNumber;
	
	bool isLQTset;
	float recBaz, recEpic, evLon, evLat, evMag;
	float raySlowness;			// ray slowness in seconds/km - converted to sec/km...
	std::string phaseName;	    // used to pick out phase identifier
							// useful for identifying PP phases during LQT rotation
	
	std::string eventStats; // load relevant event stats useful in publications
	
	std::string fName; // use to store file name
	
	// Should I place extra data for storing and displaying waveforms intermediate steps?
	
	SACDATA Data;
	MatDoub DataCmps;
	float timeTag, timeStart; // time tag
	
	//OLUGBOJI2016 -- also save raw waveform traces...
	// think useful for debugging, but also for demonstrating rotations ..
	MatDoub DataCmpsRw;
	
	template<class T>
	char *as_bytes(T& i){ //treat a T as a sequence of bytes
		void* addr = i; //get the address of the first byte
		return static_cast<char *> (addr); //treat that memory as bytes
		
	}
	
	// string concatenation class..
	template <class T>
	inline std::string to_string (const T& t)
	{
		std::stringstream ss;
		ss << t;
		return ss.str();
	}
	
private:
	/* prototype for SACIO functions */
	void	swab4(char *, int);
	int sac_byte_order() ;
	
	/* duty horse function that does all the work */
	int	read_sac(std::string);
	int  timeMark_stats();
	enum Components {
		Vertical,
		Radial,
		Transverse
	};
};
