#include "..\include\Cheb3D.hpp"

Matrix Cheb3D(int t, int N, int Ta, int Tb, Matrix& Cx, Matrix& Cy, Matrix& Cz){
	if ( (t<Ta) || (Tb<t) ){
		cout << "ERROR: Time out of range in Cheb3D::Value\n";
        exit(EXIT_FAILURE);
	}

	// Clenshaw algorithm
	double tau = (2*t-Ta-Tb)/(Tb-Ta);  

	Matrix f1 = zeros(1,3);
	Matrix f2 = zeros(1,3);

	for (i=N; i<=3; i++){
		matrix old_f1 = f1;
		f1 = 2*tau*f1-f2+[Cx(i),Cy(i),Cz(i)];
		f2 = old_f1;
	}

	Matrix ChebApp = tau*f1-f2+[Cx(1),Cy(1),Cz(1)];
	
	return ChebApp;
}