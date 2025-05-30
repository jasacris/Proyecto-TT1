#include "..\include\matrix.hpp"
#include "..\include\global.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Position.hpp"
#include "..\include\DEInteg.hpp"
#include "..\include\Accel.hpp"
#include "..\include\R_z.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\VarEqn.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\LTC.hpp"
#include <iostream>
#include <tuple>
#include <string.h>

using namespace std;

int main() {

    DE430Coeff(2285, 1020);
    
    GGM03S(181);

    eop19620101(21413);

    int nobs = 46;
    Matrix& obs = zeros(nobs,4);

    FILE *fid = fopen("../data/GEOS3.txt","r");

    if(fid == NULL){
        printf("Fail open GEOS3.txt file\n");
        exit(EXIT_FAILURE);
    }

    char tline[256];
    char sub[10];
    double YY,M,D,h,m,ss,az,el,Dist;

    for (int i = 1; i <= nobs; i++){
        fgets(tline, sizeof(tline), fid);

        if (strlen(tline) < 4){
            break;
        }

        strncpy(sub, &tline[0], 4);
        sub[4] = '\0';
        YY = atof(sub);

        strncpy(sub, &tline[5], 2);
        sub[2] = '\0';
        M = atof(sub);

        strncpy(sub, &tline[8], 2);
        sub[2] = '\0';
        D = atof(sub);

        strncpy(sub, &tline[12], 2);
        sub[2] = '\0';
        h = atof(sub);

        strncpy(sub, &tline[15], 2);
        sub[2] = '\0';
        m = atof(sub);

        strncpy(sub, &tline[18], 6);
        sub[6] = '\0';
        ss = atof(sub);

        strncpy(sub, &tline[25], 8);
        sub[8] = '\0';
        az = atof(sub);

        strncpy(sub, &tline[35], 7);
        sub[7] = '\0';
        el = atof(sub);

        strncpy(sub, &tline[44], 10);
        sub[10] = '\0';
        Dist = atof(sub);

        obs(i,1) = Mjday(YY,M,D,h,m,ss);
        obs(i,2) = Rad * az;
        obs(i,3) = Rad * el;
        obs(i,4) = 1e3 * Dist;
    }

    fclose(fid);

    double sigma_range = 92.5;
    double sigma_az = 0.0224 * Rad;
    double sigma_el = 0.0139 * Rad;

    double lat = Rad * 21.5748;
    double lon = Rad * (-158.2706);
    double alt = 300.20;

    Matrix& Rs = transponse(Position(lon, lat, alt));

    double Mjd1 = obs(1,1);
    double Mjd2 = obs(9,1);
    double Mjd3 = obs(18,1);

    Matrix &r2 = zeros(3), &v2 = zeros(3);

    r2(1) = 6221397.62857869; r2(2) = 2867713.77965738; r2(3) = 3006155.98509949;
    v2(1) = 4645.04725161806; v2(2) = -2752.21591588204; v2(3) = -7507.99940987031;
    
    //[r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
    // [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
    //                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

    Matrix& Y0_apr = transponse(r2.union_vector(v2));

    double Mjd0 = Mjday(1995,1,29,02,38,0);

    double Mjd_UTC = obs(9,1);

    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.n      = 20;
    AuxParam.m      = 20;
    AuxParam.sun     = 1;
    AuxParam.moon    = 1;
    AuxParam.planets = 1;

    int n_eqn  = 6;

    Matrix& Y = DEInteg(Accel, 0, -(obs(9,1)-Mjd0) * 86400.0, 1e-13, 1e-6, 6, transponse(Y0_apr));

    Matrix P = zeros(6,6);
    
    for (int i = 1; i <= 3; i++){
        P(i,i) = 1e8;
    }

    for (int i = 4; i <= 6; i++){
        P(i,i) = 1e3;
    }

    Matrix& LT = LTC(lon, lat);

    Matrix& yPhi = zeros(42);
    Matrix& Phi  = zeros(6,6);

    double t = 0.0;

    double t_old;

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC, Mjd_TT, Mjd_UT1, theta, Azim, Elev;

    Matrix& Y_old = zeros(1);
    Matrix& U = zeros(1);
    Matrix& r = zeros(1);
    Matrix& s = zeros(1);
    Matrix& dAds = zeros(1);
    Matrix& dEds = zeros(1);
    Matrix& dAdY = zeros(1);
    Matrix& K = zeros(1);
    Matrix& dEdY = zeros(1);
    Matrix& dDds = zeros(1);
    Matrix& dDdY = zeros(1);
    
    for (int i = 1; i <= nobs; i++){
        t_old = t;
        Y_old = Y;
        
        Mjd_UTC = obs(i,1);
        t       = (Mjd_UTC - Mjd0) * 86400.0;
        
        tie(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC) = IERS(Mjd_UTC, 'l');
        tie(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC) = timediff(UT1_UTC, TAI_UTC);
        Mjd_TT = Mjd_UTC + TT_UTC / 86400.0;
        Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
        AuxParam.Mjd_UTC = Mjd_UTC;
        AuxParam.Mjd_TT = Mjd_TT;
            
        for (int ii = 1; ii <= 6; ii++){
            yPhi(ii) = Y_old(ii);

            for (int j = 1; j <= 6; j++){
                if (ii == j){
                    yPhi(6 * j + ii) = 1.0; 
                }else{
                    yPhi(6 * j + ii) = 0.0;
                }
            }
        }
        
        yPhi = DEInteg (VarEqn, 0, t - t_old, 1e-13, 1e-6, 42, yPhi);
        
        for (int j = 1; j <= 6; j++){
            Phi.assign_column(yPhi.extract_vector(6 * j + 1, 6 * j + 6), j);
        }

        Y = DEInteg (Accel, 0.0, t - t_old, 1e-13, 1e-6, 6, Y_old);
        
        theta = gmst(Mjd_UT1);
        U = R_z(theta);
        r = transponse(Y.extract_vector(1,3));

        s = LT * (U * r - Rs);
        
        P = TimeUpdate(P, Phi);
            
        tie(Azim, Elev, dAds, dEds) = AzElPa(transponse(s));
        dAdY = (dAds * LT * U).union_vector(zeros(3));
        
        tie(K, Y, P) = MeasUpdate(Y, obs(i,2), Azim, sigma_az, dAdY, P, 6);
        
        r = transponse(Y.extract_vector(1,3));
        s = LT * (U * r - Rs);
        tie(Azim, Elev, dAds, dEds) = AzElPa(transponse(s));
        dEdY = (dEds * LT * U).union_vector(zeros(3));
        
        tie(K, Y, P) = MeasUpdate(Y, obs(i,3), Elev, sigma_el, dEdY, P, 6);
        
        r = transponse(Y.extract_vector(1,3));
        s = LT * (U * r - Rs);
        Dist = norm(transponse(s));
        dDds = transponse(s / Dist);
        dDdY = (dDds * LT * U).union_vector(zeros(3));
        
        tie(K, Y, P) = MeasUpdate(Y, obs(i,4), Dist, sigma_range, dDdY, P, 6);
    }

    tie(x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC) = IERS(obs(46,1), 'l');
    tie(UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC) = timediff(UT1_UTC, TAI_UTC);
    Mjd_TT = Mjd_UTC + TT_UTC / 86400;
    AuxParam.Mjd_UTC = Mjd_UTC;
    AuxParam.Mjd_TT = Mjd_TT;

    Matrix& Y0 = DEInteg (Accel, 0, -(obs(46,1) - obs(1,1)) * 86400.0, 1e-13, 1e-6, 6, Y);

    Matrix& Y_true = zeros(6,1);
    
    Y_true(1) = 5753.173e3;
    Y_true(2) = 2673.361e3;
    Y_true(3) = 3440.304e3;
    Y_true(4) = 4.324207e3;
    Y_true(5) = -1.924299e3;
    Y_true(6) = -5.728216e3;

    printf("\nError of Position Estimation\n");
    printf("dX%10.1f [m]\n",Y0(1)-Y_true(1));
    printf("dY%10.1f [m]\n",Y0(2)-Y_true(2));
    printf("dZ%10.1f [m]\n",Y0(3)-Y_true(3));
    printf("\nError of Velocity Estimation\n");
    printf("dVx%8.1f [m/s]\n",Y0(4)-Y_true(4));
    printf("dVy%8.1f [m/s]\n",Y0(5)-Y_true(5));
    printf("dVz%8.1f [m/s]\n",Y0(6)-Y_true(6));

    return 0;
}