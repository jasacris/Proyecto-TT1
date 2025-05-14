#ifndef _ACCELHARMONIC_
#define _ACCELHARMONIC_

#include "..\include\matrix.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\global.hpp"
#include <cmath>

using namespace std;
Matrix& AccelHarmonic(Matrix &r, Matrix &E, int n_max, int m_max);

#endif