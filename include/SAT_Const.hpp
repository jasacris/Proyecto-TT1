#ifndef _SAT_CONST_
#define _SAT_CONST_

#include "math.h"

struct Constants {
    // Constantes matemáticas
    constexpr static const double pi2 = 2 * M_PI;                    // 2pi
    const double Rad = M_PI / 180;                  // Radianes por grado
    const double Deg = 180 / M_PI;                  // Grados por radian
    const double Arcs = 3600 * 180 / M_PI;          // Segundos de arco por radian

    // Generales
    const double MJD_J2000 = 51544.5;               // Fecha Juliana Modificada de J2000
    const double T_B1950 = -0.500002108;            // Época B1950
    const double c_light = 299792458.000000000;     // Velocidad de la luz [m/s]
    const double AU = 149597870700.000000;          // Unidad astronómica [m]

    // Parámetros físicos de la Tierra, el Sol y la Luna
    const double R_Earth = 6378.1363e3;             // Radio de la Tierra [m]
    const double f_Earth = 1 / 298.257223563;       // Achatamiento de la Tierra
    const double R_Sun = 696000e3;                  // Radio del Sol [m]
    const double R_Moon = 1738e3;                   // Radio de la Luna [m]

    // Rotación de la Tierra (derivada de GMST en J2000; difiere del período inercial por precesión)
    const double omega_Earth = 15.04106717866910 / 3600 * Rad;  // [rad/s]

    // Coeficientes gravitacionales
    const double GM_Earth = 398600.435436e9;                   // [m^3/s^2]
    const double GM_Sun = 132712440041.939400e9;               // [m^3/s^2]
    const double GM_Moon = GM_Earth / 81.30056907419062;       // [m^3/s^2]
    const double GM_Mercury = 22031.780000e9;                  // [m^3/s^2]
    const double GM_Venus = 324858.592000e9;                   // [m^3/s^2]
    const double GM_Mars = 42828.375214e9;                     // [m^3/s^2]
    const double GM_Jupiter = 126712764.800000e9;              // [m^3/s^2]
    const double GM_Saturn = 37940585.200000e9;                // [m^3/s^2]
    const double GM_Uranus = 5794548.600000e9;                 // [m^3/s^2]
    const double GM_Neptune = 6836527.100580e9;                // [m^3/s^2]
    const double GM_Pluto = 977.0000000000009e9;               // [m^3/s^2]

    // Presión de radiación solar a 1 UA
    const double P_Sol = 1367 / c_light; // [N/m^2] (~1367 W/m^2)
};
#endif
