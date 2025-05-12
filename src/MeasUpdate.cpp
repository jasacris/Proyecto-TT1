#include "..\include\MeasUpdate.hpp"

tuple<Matrix&, Matrix&, Matrix&> MeasUpdate(Matrix& x, double z, double g, double s, Matrix& G, Matrix& P, int n){

    int m = 1;
    double Inv_W = s*s; 

    Matrix& K = P * transpose(G) * inv(G * P * transpose(G) + Inv_W);

    x = x + (K * (z - g));

    P = (eye(n) - K * G) * P;

    return tie(K,x,P);
}