#include "..\include\AzElPa.hpp"

tuple<double,double,Matrix&,Matrix&> AzElPa(Matrix &s){
	double pi2 = 2.0*M_PI;
	
	double rho = sqrt(s(1)*s(1)+s(2)*s(2));
	
	double Az = atan2(s(1),s(2));
	
	if(Az < 0){
		Az += pi2;
	}
	
	double El = atan ( s(3) / rho );
	
	Matrix& dAds = zeros(3);
	Matrix& dEds = zeros(3);
	
	dAds(1) = s(2)/(rho*rho);
	dAds(2) = -s(1)/(rho*rho);
	dAds(3) = 0.0;
	
	dEds(1) = (-s(1)*s(3)/rho)/dot(s,s);
	dEds(2) = (-s(2)*s(3)/rho)/dot(s,s);
	dEds(3) = rho/dot(s,s);
	
	return tie(Az, El, dAds, dEds);
}