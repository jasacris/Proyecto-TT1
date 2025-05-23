#ifndef _IERS_
#define _IERS_

#include <tuple>
#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\global.hpp"

using namespace std;

tuple<double,double,double,double,double,double,double,double,double> IERS(double Mjd_UTC, char interp = 'n');

#endif