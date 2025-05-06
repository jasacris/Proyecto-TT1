#include "..\include\EccAnom.hpp"

double EccAnom(double M, double e){
    int maxit = 15;
    int i = 1;

    M = M % (2.0*M_PI);

    double E;
    if(e < 0.8){
        E = M;
    }else{
        E = M_PI;
    }

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    while(abs(f) > 1e2*eps){
        f = E - e * sin(E) - M;
        E = E - f / (1.0 - e * cos(E));
        i++;
        if (i == maxit) {
            cout << "Problemas de convergencia en EccAnom\n";
            exit(EXIT_FAILURE);
        }
    }

    return E;
}