#include "..\include\IERS.hpp"

tuple<double,double,double,double,double,double,double,double,double> IERS(Matrix &eop, double Mjd_UTC, char interp){
    Constants constants;

    double x_pole;
    double y_pole;
    double UT1_UTC;
    double LOD;
    double dpsi;
    double deps;
    double dx_pole;
    double dy_pole;
    double TAI_UTC;

    if (interp =='l'){
        double mjd = (floor(Mjd_UTC));
        int i;

        for (int col = 1; col <= eop.n_column; col++) {
            if (mjd == eop(4, col)) {
                i = col;
                break;
            }
        }

        Matrix preeop = eop.extract_column(i);
        Matrix nexteop = eop.extract_column(i+1);

        double mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
        double = mfme/1440;

        x_pole  = preeop(5)+(nexteop(5)-preeop(5))*fixf;
        y_pole  = preeop(6)+(nexteop(6)-preeop(6))*fixf;
	    UT1_UTC = preeop(7)+(nexteop(7)-preeop(7))*fixf;
        LOD     = preeop(8)+(nexteop(8)-preeop(8))*fixf;
        dpsi    = preeop(9)+(nexteop(9)-preeop(9))*fixf;
        deps    = preeop(10)+(nexteop(10)-preeop(10))*fixf;
        dx_pole = preeop(11)+(nexteop(11)-preeop(11))*fixf;
        dy_pole = preeop(12)+(nexteop(12)-preeop(12))*fixf;
        TAI_UTC = preeop(13);
        
        x_pole  = x_pole/constants.Arcs;  % Pole coordinate [rad]
        y_pole  = y_pole/constants.Arcs;  % Pole coordinate [rad]
        dpsi    = dpsi/constants.Arcs;
        deps    = deps/constants.Arcs;
        dx_pole = dx_pole/constants.Arcs; % Pole coordinate [rad]
        dy_pole = dy_pole/constants.Arcs; % Pole coordinate [rad]
    }else if (interp =='n') {
        double mjd = (floor(Mjd_UTC));
        int i;

        for (int col = 1; col <= eop.n_column; col++) {
            if (mjd == eop(4, col)) {
                i = col;
                break;
            }
        }
        eop = eop.extract_column(i);
        
        x_pole  = eop(5)/constants.Arcs; 
        y_pole  = eop(6)/constants.Arcs;  
	    UT1_UTC = eop(7);          
        LOD     = eop(8);             
        dpsi    = eop(9)/constants.Arcs;
        deps    = eop(10)/constants.Arcs;
        dx_pole = eop(11)/constants.Arcs; 
        dy_pole = eop(12)/constants.Arcs; 
	    TAI_UTC = eop(13);          
    }

    return tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);
}