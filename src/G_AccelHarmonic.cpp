/*#include "..\include\G_AccelHarmonic.hpp"

Matrix& G_AccelHarmonic(Matrix& r, Matrix& U, int n_max, int m_max){
    double d = 1.0; 

    Matrix& G = zeros(3,3);
    Matrix& dr = zeros(3);

    for (int i = 1; i <= 3; i++){
        dr = zeros(3);

        dr(i) = d;

        Matrix& da = AccelHarmonic (r + dr / 2, U, n_max, m_max) - AccelHarmonic (r - dr / 2, U, n_max, m_max);
        
        G.assign_column(da / d, i);      
    }

    return G;
}*/