#include "..\include\gast.hpp"

double gast(double Mjd_UT1){

    double gstime = fmod (gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), pi2);

    return gstime;
}