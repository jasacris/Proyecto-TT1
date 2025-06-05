#include "..\include\MeasUpdate.hpp"

tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix& x, double z, double g, double s, Matrix& G, Matrix& P, int n){

    double Inv_W = s * s; 

    Matrix& K = transponse(P * transponse(G) * inv(G * P * transponse(G) + Inv_W));

    x = x + (K * (z - g));

    P = (eye(n) - transponse(K) * G) * P;

    return tie(K,x,P);
}