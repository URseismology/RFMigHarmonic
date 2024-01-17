/**
 *  @file sachead.h
 *  
 *
 *  @author Tolulope Olugboji 
 *	@date 11/22/10.
 *  @copyright Yale University. All rights reserved.
 *
 *  @brief Much of the details of this header template is adapted from the sac.h header file in the utils directory of the sac distribution
 */

/* True/false definitions */
#ifndef TRUE
#define FALSE	0
#define TRUE	1
#endif

#define SAC_HEADER_FIELDS          133
#define SAC_HEADER_SIZE_NUMBERS    440
#define SAC_HEADER_SIZE            632 
#define SAC_HEADER_TMARK_POSITION  10
#define SAC_HEADER_USERN_POSITION  40

#define SAC_HEADER_FLOAT_MIN       0
#define SAC_HEADER_FLOAT_MAX       69
#define SAC_HEADER_INT_MIN         70
#define SAC_HEADER_INT_MAX         84
#define SAC_HEADER_ENUM_MIN        85
#define SAC_HEADER_ENUM_MAX        104
#define SAC_HEADER_LOGICAL_MIN     105
#define SAC_HEADER_LOGICAL_MAX     109
#define SAC_HEADER_CHAR_MIN        110
#define SAC_HEADER_CHAR_MAX        133
#define SAC_HEADER_CHAR_DOUBLE     111
#define SAC_HEADER_CHAR_DOUBLE_END 112

//char *SacHeaderNameNull = "";

#define SAC_HEADER_FLOAT_UNDEFINED (-12345.0)
#define SAC_HEADER_INT_UNDEFINED   (-12345)
#define SAC_HEADER_CHAR_UNDEFINED  ("-12345  ")
#define SAC_HEADER_UNDEFINED       ("UNDEFINED")

typedef struct sac_head
{
	float	delta;			/** RF time increment, sec    */
	float	depmin;			/**    minimum amplitude      */
	float	depmax;			/**    maximum amplitude      */
	float	scale;			/**    amplitude scale factor */
	float	odelta;			/**    observed time inc      */
	float	b;			/** RD initial time - wrt nz* */
	float	e;			/** RD end time               */
	float	o;			/**    event start            */
	float	a;			/**    1st arrival time       */
	float	fmt;		        /**    internal use           */
	float	t0;			/**    user-defined time pick */
	float	t1;			/**    user-defined time pick */
	float	t2;			/**    user-defined time pick */
	float	t3;			/**    user-defined time pick */
	float	t4;			/**    user-defined time pick */
	float	t5;			/**    user-defined time pick */
	float	t6;			/**    user-defined time pick */
	float	t7;			/**    user-defined time pick */
	float	t8;			/**    user-defined time pick */
	float	t9;			/**    user-defined time pick */
	float	f;			/**    event end, sec > 0     */
	float	resp0;			/**    instrument respnse parm*/
	float	resp1;			/**    instrument respnse parm*/
	float	resp2;			/**   instrument respnse parm*/
	float	resp3;			/**    instrument respnse parm*/
	float	resp4;			/**    instrument respnse parm*/
	float	resp5;			/**    instrument respnse parm*/
	float	resp6;			/**    instrument respnse parm*/
	float	resp7;			/**    instrument respnse parm*/
	float	resp8;			/**    instrument respnse parm*/
	float	resp9;			/**    instrument respnse parm*/
	float	stla;			/**  T station latititude     */
	float	stlo;			/**  T station longitude      */
	float	stel;			/**  T station elevation, m   */
	float	stdp;			/**  T station depth, m       */
	float	evla;			/**    event latitude         */
	float	evlo;			/**    event longitude        */
	float	evel;			/**    event elevation        */
	float	evdp;			/**    event depth            */
	float	mag;    		/**    magnitude value        */
	float	user0;			/**    available to user      */
	float	user1;			/**    available to user      */
	float	user2;			/**    available to user      */
	float	user3;			/**    available to user      */
	float	user4;			/**    available to user      */
	float	user5;			/**    available to user      */
	float	user6;			/**    available to user      */
	float	user7;			/**    available to user      */
	float	user8;			/**    available to user      */
	float	user9;			/**    available to user      */
	float	dist;			/**    stn-event distance, km */
	float	az;			/**    event-stn azimuth      */
	float	baz;			/**    stn-event azimuth      */
	float	gcarc;			/**    stn-event dist, degrees*/
	float	sb;     		/**    saved b value          */
	float	sdelta; 		/**    saved delta value      */
	float	depmen;			/**    mean value, amplitude  */
	float	cmpaz;			/**  T component azimuth      */
	float	cmpinc;			/**  T component inclination  */
	float	xminimum;		/**    XYZ X minimum value    */
	float	xmaximum;		/**    XYZ X maximum value    */
	float	yminimum;		/**    XYZ Y minimum value    */
	float	ymaximum;		/**    XYZ Y maximum value    */
	float	unused6;		/**    reserved for future use*/
	float	unused7;		/**    reserved for future use*/
	float	unused8;		/**    reserved for future use*/
	float	unused9;		/**    reserved for future use*/
	float	unused10;		/**    reserved for future use*/
	float	unused11;		/**    reserved for future use*/
	float	unused12;		/**    reserved for future use*/
	int	nzyear;			/**  F zero time of file, yr  */
	int	nzjday;			/**  F zero time of file, day */
	int	nzhour;			/**  F zero time of file, hr  */
	int	nzmin;			/**  F zero time of file, min */
	int	nzsec;			/**  F zero time of file, sec */
	int	nzmsec;			/**  F zero time of file, msec*/
	int	nvhdr;  		/**  R header version number  */
	int	norid;  		/**    Origin ID              */
	int	nevid;  		/**    Event ID               */
	int	npts;			/** RF number of samples      */
	int	nsnpts; 		/**    saved npts             */
	int	nwfid; 		        /**    Waveform ID            */
	int	nxsize;	                /**    XYZ X size             */
	int	nysize;   		/**    XYZ Y size             */
	int	unused15;		/**    reserved for future use*/
	int	iftype;			/** RA type of file          */
	int	idep;			/**    type of amplitude      */
	int	iztype;			/**    zero time equivalence  */
	int	unused16;		/**    reserved for future use*/
	int	iinst;			/**    recording instrument   */
	int	istreg;			/**    stn geographic region  */
	int	ievreg;			/**    event geographic region*/
	int	ievtyp;			/**    event type             */
	int	iqual;			/**    quality of data        */
	int	isynth;			/**    synthetic data flag    */
	int	imagtyp;		/**    magnitude type         */
	int	imagsrc;		/**    magnitude source       */
	int	unused19;		/**    reserved for future use*/
	int	unused20;		/**    reserved for future use*/
	int	unused21;		/**    reserved for future use*/
	int	unused22;		/**    reserved for future use*/
	int	unused23;		/**    reserved for future use*/
	int	unused24;		/**    reserved for future use*/
	int	unused25;		/**    reserved for future use*/
	int	unused26;		/**    reserved for future use*/
	int	leven;			/** RA data-evenly-spaced flag*/
	int	lpspol;			/**    station polarity flag  */
	int	lovrok;			/**    overwrite permission   */
	int	lcalda;			/**    calc distance, azimuth */
	int	unused27;		/**    reserved for future use*/
	char	kstnm[8];		/**  F station name           */
	char	kevnm[16];		/**    event name             */
	char	khole[8];		/**    man-made event name    */
	char	ko[8];			/**    event origin time id   */
	char	ka[8];			/**    1st arrival time ident */
	char	kt0[8];			/**    time pick 0 ident      */
	char	kt1[8];			/**    time pick 1 ident      */
	char	kt2[8];			/**    time pick 2 ident      */
	char	kt3[8];			/**    time pick 3 ident      */
	char	kt4[8];			/**    time pick 4 ident      */
	char	kt5[8];			/**    time pick 5 ident      */
	char	kt6[8];			/**    time pick 6 ident      */
	char	kt7[8];			/**    time pick 7 ident      */
	char	kt8[8];			/**    time pick 8 ident      */
	char	kt9[8];			/**    time pick 9 ident      */
	char	kf[8];			/**    end of event ident     */
	char	kuser0[8];		/**    available to user      */
	char	kuser1[8];		/**    available to user      */
	char	kuser2[8];		/**    available to user      */
	char	kcmpnm[8];		/**  F component name         */
	char	knetwk[8];		/**    network name           */
	char	kdatrd[8];		/**    date data read         */
	char	kinst[8];		/**    instrument name        */
} SACHEAD;

typedef struct sac_data{
	float data[100000];
} SACDATA;
