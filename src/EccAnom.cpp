#include "..\include\EccAnom.hpp"

double EccAnom(double M, double e){
    int maxit = 15;
    int i = 1;

    M = fmod(M, pi2);

    double E;
    if(e < 0.8){
        E = M;
    }else{
        E = pi;
    }

    double f = E - e * sin(E) - M;
    E = E - f / (1.0 - e * cos(E));

    while(fabs(f) > 1e2*std::numeric_limits<double>::epsilon()){
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