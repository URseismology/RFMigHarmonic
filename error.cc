/**
 * @file error.cc
 *
 * @author Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include <iostream>
#include <cstdlib>
#include "error.h"

void error(string s)
{
    cerr << s << '\n';
    exit(1);
}

