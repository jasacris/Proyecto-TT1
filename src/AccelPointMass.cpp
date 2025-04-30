#include "..\include\AccelPointMass.hpp"

double AccelPointMass(Matrix r, Matrix s, double GM){
	
	double d = r - s;
	
	return -GM * ((d/pow(norm(d),3))+(s/pow(norm(s),3));
}