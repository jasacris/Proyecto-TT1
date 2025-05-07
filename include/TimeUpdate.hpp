#ifndef _TIMEUPDATE_
#define _TIMEUPDATE_

#include "..\include\matrix.hpp"
#include <cmath>

Matrix TimeUpdate(Matrix &P, Matrix &Phi, Matrix &Qdt = 0.0);

#endif