/*#include "..\include\MeanObliquity.hpp"

double MeanObliquity(double Mjd_TT){
    Constants constants;

    double T = (Mjd_TT - constants.MJD_J2000) / 36525.0;
    
    double MOblq = constants.Rad * (84381.448 / 3600.0 - (46.8150 + (0.00059 - 0.001813 * T) * T) * T / 3600.0);
    return MOblq;
}*/