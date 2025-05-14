#ifndef _GLOBAL_
#define _GLOBAL_

#include <cmath>
#include "..\include\matrix.hpp"

typedef struct{
	double Mjd_UTC, Mjd_TT;
	int n, m, sun, moon, planets;
} Param;

extern Param AuxParam;
extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;

void eop19620101(int c);

void GGM03S(int n);

void DE430Coeff(int f, int c);

#endif