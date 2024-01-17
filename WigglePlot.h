/*
 *  WigglePlot.h
 *  
 *
 *  Created by Tolulope Olugboji on June 17 2011.
 *  Copyright Yale University. All rights reserved.
 *
 *  Test plplot utility in this class driver. See if you can display seismograms with the 
 *  xwindow utility
 */

#include "plc++demos.h"

#include "plevent.h"
#include <cctype>

#ifdef PL_HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef PL_USE_NAMESPACE
using namespace std;
#endif

static PLGraphicsIn gin;

static int          locate_mode;
static int          test_xor;
static int          fontset;
static char         *f_name;

static const char   *notes[] = { "Make sure you get it right!", "" };

/* Options data structure definition. */
static PLOptionTable options[] = {
    {
        "locate",               /* Turns on test of API locate function */
        NULL,
        NULL,
        &locate_mode,
        PL_OPT_BOOL,
        "-locate",
        "Turns on test of API locate function"
    },
    {
        "xor",                  /* Turns on test of xor function */
        NULL,
        NULL,
        &test_xor,
        PL_OPT_BOOL,
        "-xor",
        "Turns on test of XOR"
    },
    {
        "font",                 /* For switching between font set 1 & 2 */
        NULL,
        NULL,
        &fontset,
        PL_OPT_INT,
        "-font number",
        "Selects stroke font set (0 or 1, def:1)"
    },
    {
        "save",                 /* For saving in postscript */
        NULL,
        NULL,
        &f_name,
        PL_OPT_STRING,
        "-save filename",
        "Save plot in color postscript `file'"
    },
    {
        NULL,                   /* option */
        NULL,                   /* handler */
        NULL,                   /* client data */
        NULL,                   /* address of variable to set */
        0,                      /* mode flag */
        NULL,                   /* short syntax */
        NULL
    }                           /* long syntax */
};


class WigglePlot {
public:
    WigglePlot( int, const char ** );
	
    // Stack Varying plot interfaces... sweet structure.
	// void plot1( int );
    // void plot2();
    // void plot3();
	
	
private:
    // Class data
    PLFLT    xscale, yscale, xoff, yoff;
    plstream *pls;
};