#ifndef _AZELPA_
#define _AZELPA_

#define _USE_MATH_DEFINES

#include <math.h>

#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <tuple>

using namespace std;

tuple<double,double,Matrix&,Matrix&> AzElPa(Matrix &s);

#endif