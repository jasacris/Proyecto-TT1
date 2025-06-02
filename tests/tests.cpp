#include "..\include\matrix.hpp"
#include "..\include\global.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\Frac.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\Position.hpp"
#include "..\include\sign_.hpp"
#include "..\include\timediff.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\IERS.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\LTC.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\gmst.hpp"
#include "..\include\gast.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\Accel.hpp"
#include "..\include\VarEqn.hpp"
#include "..\include\DEInteg.hpp"

#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sum_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  4; C(1,3) = 10; C(1,4) = 2;
	C(2,1) = 3; C(2,2) = 1; C(2,3) = 2; C(2,4) = 2;
	C(3,1) = 2; C(3,2) = 3; C(3,3) = 2; C(3,4) = 7;
	
	Matrix R = A + 2.0;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) = 0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_sub_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) = 0; C(1,3) = 6; C(1,4) = -2;
	C(2,1) = -1; C(2,2) = -3; C(2,3) = -2; C(2,4) = -2;
	C(3,1) = -2; C(3,2) = -1; C(3,3) = -2; C(3,4) = 3;
	
	Matrix R = A - 2.0;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_prod_01() {

    Matrix A(2, 3);
    A(1, 1) = 1; A(1, 2) = 2; A(1, 3) = 3;
    A(2, 1) = 4; A(2, 2) = 5; A(2, 3) = 6;

    Matrix B(3, 2);
    B(1, 1) = 7; B(1, 2) = 8;
    B(2, 1) = 9; B(2, 2) = 10;
    B(3, 1) = 11; B(3, 2) = 12;

	Matrix C(2, 2);
	C(1,1) = 58; C(1,2) = 64;
	C(2,1) = 139; C(2,2) = 154;

    Matrix R = A * B;
	
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_prod_02() {

    Matrix A(2, 3);
    A(1, 1) = 1; A(1, 2) = 2; A(1, 3) = 3;
    A(2, 1) = 4; A(2, 2) = 5; A(2, 3) = 6;

	Matrix C(2, 3);
    C(1, 1) = 2; C(1, 2) = 4; C(1, 3) = 6;
    C(2, 1) = 8; C(2, 2) = 10; C(2, 3) = 12;

    Matrix R = A * 2.0;
	
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_div_01() {

    Matrix A(3, 3);
    A(1, 1) = 0; A(1, 2) = 2; A(1, 3) = 1;
    A(2, 1) = 1; A(2, 2) = -1; A(2, 3) = 0;
    A(3, 1) = 0; A(3, 2) = 1; A(3, 3) = 0;

    Matrix B(3, 3);
    B(1, 1) = 0; B(1, 2) = 2; B(1, 3) = 1;
    B(2, 1) = 1; B(2, 2) = -1; B(2, 3) = 0;
    B(3, 1) = 0; B(3, 2) = 1; B(3, 3) = 0;

	Matrix C(3, 3);
    C(1, 1) = 1; C(1, 2) = 0; C(1,3) = 0;
    C(2, 1) = 0; C(2, 2) = 1; C(2,3) = 0;
    C(3, 1) = 0; C(3, 2) = 0; C(3,3) = 1;

    Matrix R = A / B;
	
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_div_02() {
	
	Matrix A(2, 3);
    A(1, 1) = 2; A(1, 2) = 4; A(1, 3) = 6;
    A(2, 1) = 8; A(2, 2) = 10; A(2, 3) = 12;

    Matrix B(2, 3);
    B(1, 1) = 1; B(1, 2) = 2; B(1, 3) = 3;
    B(2, 1) = 4; B(2, 2) = 5; B(2, 3) = 6;


    Matrix R = A / 2.0;
	
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_equ_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 0; B(1,2) =  2; B(1,3) = 8; B(1,4) = 0;
	B(2,1) = 1; B(2,2) = -1; B(2,3) = 0; B(2,4) = 0;
	B(3,1) = 0; B(3,2) =  1; B(3,3) = 0; B(3,4) = 5;
	
	
	Matrix R = A;
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_inv_01() {
    int f = 3;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 1;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0;
	
	Matrix B(f, c);
	B(1,1) = 0; B(1,2) =  1; B(1,3) = 1;
	B(2,1) = 0; B(2,2) =  0; B(2,3) = 1;
	B(3,1) = 1; B(3,2) =  0; B(3,3) = -2;
	
	Matrix R = inv(A);
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_trans_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(c, f);
	B(1,1) = 0; B(1,2) =  1; B(1,3) = 0;
	B(2,1) = 2; B(2,2) = -1; B(2,3) = 1;
	B(3,1) = 8; B(3,2) =  0; B(3,3) = 0;
	B(4,1) = 0; B(4,2) =  0; B(4,3) = 5;
	
	Matrix R = transponse(A);
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}

int m_eye_01() {
	
	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = 1; A(2,3) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 1;
	
	Matrix R = eye(3);
    
    _assert(m_equals(A, R, 1e-10));
    
    return 0;
}

int m_norm_01() {
	
	Matrix A(3);
	A(1) = 4; A(2) = 0; A(3) = 3;
	
	double ans = 5;
	
	double R = norm(A);
    
    _assert(fabs(R - ans) < 1e-10);
    
    return 0;
}

int m_dot_01() {
	
	Matrix A(3);
	A(1) = 1; A(2) = 2; A(3) = 3;
	
	Matrix B(3);
	B(1) = 3; B(2) = 2; B(3) = 1;
	
	double ans = 10;
	
	double R = dot(A, B);
    
    _assert(fabs(R - ans) < 1e-10);
    
    return 0;
}

int m_cross_01() {
	
	Matrix A(1, 3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	
	Matrix B(1, 3);
	B(1,1) = 3; B(1,2) = 2; B(1,3) = 1;
	
	Matrix C(1, 3);
	C(1,1) = -4; C(1,2) = 8; C(1,3) = -4;
	
	Matrix R = cross(A, B);
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_assign_column_01() {

	Matrix A(1, 3);
	A(1,1) = 1; A(1, 2) = 2; A(1,3) = 3;
	
	Matrix B(3, 1);
	B(1,1) = 4; B(2,1) = 5; B(3,1) = 6;
	
	Matrix R = transponse(B.assign_column(A, 1));
	
	_assert(m_equals(A, R, 1e-10));
	
	return 0;
}

int m_assign_row_01() {
	
	Matrix A(1, 3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	
	Matrix B(1, 3);
	B(1,1) = 4; B(1,2) = 5; B(1,3) = 6;
	
	Matrix R = B.assign_row(A, 1);
    
    _assert(m_equals(A, R, 1e-10));
    
    return 0;
}

int m_extract_column_01() {

	Matrix A(3, 1);
	A(1,1) = 1; A(2,1) = 2; A(3,1) = 3;
	
	Matrix R = transponse(A.extract_column(1));
	
	_assert(m_equals(A, R, 1e-10));
	
	return 0;
}

int m_extract_row_01() {
	
	Matrix A(1, 3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	
	Matrix R = A.extract_row(1);
    
    _assert(m_equals(A, R, 1e-10));
    
    return 0;
}

int m_union_vector_01() {
	
	Matrix A(1, 3);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	
	Matrix B(1, 3);
	B(1,1) = 4; B(1,2) = 5; B(1,3) = 6;

	Matrix C(1, 6);
	C(1,1) = 1; C(1,2) = 2; C(1,3) = 3;
	C(1,4) = 4; C(1,5) = 5; C(1,6) = 6;
	
	Matrix R = A.union_vector(B);
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_extract_vector_01() {
	
	Matrix A(1, 5);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3; A(1,4) = 4; A(1, 5) = 5;
	
	Matrix B(1, 3);
	B(1,1) = 2; B(1,2) = 3; B(1,3) = 4;
	
	Matrix R = A.extract_vector(2, 4);
    
    _assert(m_equals(B, R, 1e-10));
    
    return 0;
}


int m_r_x_01() {
	
	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 0; A(1,3) = 0;
	A(2,1) = 0; A(2,2) = -0.989992496600445; A(2,3) =  0.141120008059867;
	A(3,1) = 0; A(3,2) = -0.141120008059867; A(3,3) = -0.989992496600445;
	
	Matrix R = R_x(3);
    
    _assert(m_equals(A, R, 1e-10));
    
    return 0;
}

int m_r_y_01() {
	
	Matrix A(3, 3);
	A(1,1) = -0.989992496600445; A(1,2) = 0; A(1,3) = -0.141120008059867;
	A(2,1) =  0; 				 A(2,2) = 1; A(2,3) = 0;
	A(3,1) =  0.141120008059867; A(3,2) = 0; A(3,3) = -0.989992496600445;
	
	Matrix R = R_y(3);
    
    _assert(m_equals(A, R, 1e-10));
    
    return 0;
}

int m_r_z_01() {
	
	Matrix A(3, 3);
	A(1,1) = -0.989992496600445; A(1,2) =  0.141120008059867; A(1,3) = 0;
	A(2,1) = -0.141120008059867; A(2,2) = -0.989992496600445; A(2,3) = 0;
	A(3,1) =  0; 				 A(3,2) =  0; 				  A(3,3) = 1;
	
	Matrix R = R_z(3);
    
    _assert(m_equals(A, R, 1e-10));
    
    return 0;
}

int m_AccelPointMass_01(){
    Matrix r(3);
    r(1)=1; r(2)=2; r(3)=3;

    Matrix s(3);
    s(1)=4; s(2)=5; s(3)=6;

	Matrix C(3);
	C(1) = 0.0463899400720928; C(2) =  0.0419499176126264; C(3) = 0.0375098951531599;

    Matrix R = AccelPointMass(r,s,3);

    _assert(m_equals(C, R, 1e-10));
    return 0;
}

int m_Cheb_01() {
	
	int N = 3;

	double t = 0.5;
	double Ta = 0;
	double Tb = 1;

	Matrix Cx(3);
	Cx(1) = 1; Cx(2) =  2; Cx(3) = 2;

	Matrix Cy(3);
	Cy(1) = 1; Cy(2) =  2; Cy(3) = 2;

	Matrix Cz(3);
	Cz(1) = 1; Cz(2) =  2; Cz(3) = 2;

	Matrix C(3);
	C(1) = -1; C(2) =  -1; C(3) = -1;
	
	Matrix R = Cheb3D(t, N, Ta, Tb, Cx, Cy, Cz);
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_EccAnom_01() {

	double M = 2.5;
	double e = 3.5;

	double ans = 2.99863819328014;
	
	double R = EccAnom(M, e);
    
    _assert(fabs(R - ans) < 1e-10);
    
    return 0;
}

int m_Frac_01() {

	double x = 5.6;

	double ans = 0.6;

	double R = Frac(x);
    
    _assert(fabs(R - ans) < 1e-10);
    
    return 0;
}

int m_MeanObli_01() {

	double x = 5;

	double ans = 0.409413039119721;

	double R = MeanObliquity(x);
    
    _assert(fabs(R - ans) < 1e-10);
    
    return 0;
}

int m_Mjday_01() {

	double ans = -674450.5762731482;

	double R = Mjday(12, 05, 02, 10, 10, 10);
    
    _assert(fabs(R - ans) < 1e-9);
    
    return 0;
}

int m_Mjday_02() {

	double ans = -674451;

	double R = Mjday(12, 05, 02);
    
    _assert(fabs(R - ans) < 1e-9);
    
    return 0;
}

int m_MjdayTDB_01() {

	double ans = 4.999999987541869;

	double R = Mjday_TDB(5);
    
    _assert(fabs(R - ans) < 1e-9);
    
    return 0;
}

int m_Position_01() {
	
	double lon = 0.5;
	double lat = 0;
	double h = 1;

	Matrix C(3);
	C(1) = 5597342.07182254; C(2) =  3057841.91034406; C(3) = 0;
	
	Matrix R = Position(lon, lat, h);
    
    _assert(m_equals(C, R, 1e-7));
    
    return 0;
}

int m_sign__01() {

	double ans = 2;

	double R = sign_(2, 2);
    
    _assert(fabs(R - ans) < 1e-9);
    
    return 0;
}

int m_sign__02() {

	double ans = -2;

	double R = sign_(2, -2);
    
    _assert(fabs(R - ans) < 1e-9);
    
    return 0;
}

int m_timediff_01(){

    double UT1_UTC = 3;
    double TAI_UTC = 5;

    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);

    double ans_UT1_TAI = -2;
    double ans_UTC_GPS = 14;
    double ans_UT1_GPS = 17;
    double ans_TT_UTC = 37.184;
    double ans_GPS_UTC = -14;

    _assert(fabs(UT1_TAI - ans_UT1_TAI)< 1e-10);
    _assert(fabs(UTC_GPS - ans_UTC_GPS)< 1e-10);
    _assert(fabs(UT1_GPS - ans_UT1_GPS)< 1e-10);
    _assert(fabs(TT_UTC - ans_TT_UTC)< 1e-10);
    _assert(fabs(GPS_UTC - ans_GPS_UTC)< 1e-10);
	
    return 0;
}

int m_AzElPa_01(){

    Matrix s(3);
	
	s(1) = 1; s(2) = 2; s(3) = 3;

    auto [Az, El, dAds, dEds] = AzElPa(s);

    double ans_Az = 0.463647609000806;
    double ans_El = 0.930274014115472;

    Matrix ans_dAds(3);

	ans_dAds(1) = 0.4; ans_dAds(2) = -0.2; ans_dAds(3) = 0.0;
	
    Matrix ans_dEds(3);

	ans_dEds(1) = -0.095831484749991; ans_dEds(2) = -0.191662969499982; ans_dEds(3) = 0.159719141249985;

    _assert(fabs(Az - ans_Az)< 1e-10);
    _assert(fabs(El - ans_El)< 1e-10);
    _assert(m_equals(dAds, ans_dAds, 1e-10));
    _assert(m_equals(dEds, ans_dEds, 1e-10));
	
    return 0;
}

int m_IERS_01(){
    
    double Mjd_UTC = 49746.1163541665;
    char interp = 'l';

    auto [x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC] = IERS(Mjd_UTC, interp);

    double ans_x_pole = -5.5937872420407e-07;
    double ans_y_pole = 2.33559834147197e-06;
    double ans_UT1_UTC = 0.325747632958709;
    double ans_LOD = 0.00272698971874332;
    double ans_dpsi = -1.16882953161744e-07;
    double ans_deps = -2.4783506198648e-08;
    double ans_dx_pole = -8.43027359626024e-10;
    double ans_dy_pole = -1.56811369105037e-09;
    double ans_TAI_UTC = 29;

    _assert(fabs(x_pole - ans_x_pole)< 1e-10);
    _assert(fabs(y_pole - ans_y_pole)< 1e-10);
    _assert(fabs(UT1_UTC - ans_UT1_UTC)< 1e-10);
    _assert(fabs(LOD - ans_LOD)< 1e-10);
    _assert(fabs(dpsi - ans_dpsi)< 1e-10);
    _assert(fabs(deps - ans_deps)< 1e-10);
    _assert(fabs(dx_pole - ans_dx_pole)< 1e-10);
    _assert(fabs(dy_pole - ans_dy_pole)< 1e-10);
    _assert(fabs(TAI_UTC - ans_TAI_UTC)< 1e-10);

    return 0;
}

int m_Legendre_01(){

    auto [pnm, dpnm] = Legendre(1, 2, 3);

    Matrix ans_pnm(2,3);

	ans_pnm(1,1) = 1.0; ans_pnm(1,2) = 0.0; ans_pnm(1,3) = 0.0;
	ans_pnm(2,1) = 0.244427023924219; ans_pnm(2,2) = -1.71471730322393; ans_pnm(2,3) = 0.0;
	
    Matrix ans_dpnm(2,3);

	ans_dpnm(1,1) = 0.0; ans_dpnm(1,2) = 0.0; ans_dpnm(1,3) = 0.0;
	ans_dpnm(2,1) = -1.71471730322393; ans_dpnm(2,2) = -0.244427023924219; ans_dpnm(2,3) = 0.0;

    _assert(m_equals(pnm, ans_pnm, 1e-10));
    _assert(m_equals(dpnm, ans_dpnm, 1e-10));
	
    return 0;
}

int m_NutAngles_01(){

    auto [dpsi, deps] = NutAngles(2);

    double ans_dpsi = 2.7179807523643e-05;
    double ans_deps = 3.91872311875582e-05;

    _assert(fabs(dpsi - ans_dpsi)< 1e-10);
    _assert(fabs(deps - ans_deps)< 1e-10);
	
    return 0;
}

int m_TimeUpdate_01() {

	Matrix P(2,2);
    P(1,1) = 1; P(1,2) = 2;
    P(2,1) = 3; P(2,2) = 4;

    Matrix Phi(2,2);
    Phi(1,1) = 5; Phi(1,2) = 6;
    Phi(2,1) = 7; Phi(2,2) = 8;

    double Qdt = 5.5;
    
    Matrix R = TimeUpdate(P, Phi, Qdt);

	Matrix ans(2,2);
    ans(1,1) = 324.5; ans(1,2) = 438.5;
    ans(2,1) = 436.5; ans(2,2) = 590.5;
    
    _assert(m_equals(R, ans, 1e-10));
    
    return 0;
}

int m_TimeUpdate_02() {

	Matrix P(2,2);
    P(1,1) = 1; P(1,2) = 2;
    P(2,1) = 3; P(2,2) = 4;

    Matrix Phi(2,2);
    Phi(1,1) = 5; Phi(1,2) = 6;
    Phi(2,1) = 7; Phi(2,2) = 8;
    
    Matrix R = TimeUpdate(P, Phi);

	Matrix ans(2,2);
    ans(1,1) = 319; ans(1,2) = 433;
    ans(2,1) = 431; ans(2,2) = 585;
    
    _assert(m_equals(R, ans, 1e-10));
    
    return 0;
}

int m_AccelHarmonic_01() {

	Matrix r(3);
    r(1) = 1; r(2) = 2; r(3) = 0;

    Matrix E(3,3);
    E(1,1) = 4; E(1,2) = 0; E(1,3) = 5;
    E(2,1) = 6; E(2,2) = 7; E(2,3) = 0;
    E(3,1) = 5; E(3,2) = 1; E(3,3) = 0;
    
    Matrix R = AccelHarmonic(r, E, 2, 2);

	Matrix ans(3);
    ans(1) = -8.43253296816925e+20;
    ans(2) = -4.64625239595692e+20;
	ans(3) = -5.2988491468403e+19;
	
    _assert(m_equals(R, ans, 1e+10));
    
    return 0;
}

int m_EqnEquinox_01(){

    double R = EqnEquinox(2);

    double ans = 2.49335221515174e-05;

    _assert(fabs(R - ans)< 1e-10);
	
    return 0;
}

int m_JPL_Eph_DE430_01(){

    Matrix x(2,2);
	
	x(1,1) = 1; x(1,2) = 2;
    x(2,1) = 3; x(2,2) = 4; 
    
    Matrix G(2,2);
	
	G(1,1) = 4; G(1,2) = 5;
    G(2,1) = 6; G(2,2) = 7; 

    Matrix P(2,2);
	
	P(1,1) = 2; P(1,2) = 6;
    P(2,1) = 4; P(2,2) = 1; 

    auto [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE430(Mjday_TDB(AuxParam.Mjd_TT));

    Matrix ans_Mercury(3);

	ans_Mercury(1) = 83760631595.4308; ans_Mercury(2) = -65287281685.7782; ans_Mercury(3) = -23387398482.2596;
	
    Matrix ans_Venus(3);

	ans_Venus(1) = -15219008641.7309; ans_Venus(2) = -110141489204.256; ans_Venus(3) = -41024866021.001;

    Matrix ans_Earth(3);

	ans_Earth(1) = -92481354617.5946; ans_Earth(2) = 106387374396.777; ans_Earth(3) = 46126868979.9791;
	
    Matrix ans_Mars(3);

	ans_Mars(1) = -88274887128.7535; ans_Mars(2) = 46966034914.9639; ans_Mars(3) = 29071603267.4181;

	Matrix ans_Jupiter(3);
	
	ans_Jupiter(1) = -298370618068.54; ans_Jupiter(2) = -754492984903.812; ans_Jupiter(3) = -314408340637.005;
	
	Matrix ans_Saturn(3);
	
	ans_Saturn(1) = 1482045333721; ans_Saturn(2) = -453861584762.029; ans_Saturn(3) = -249397462197.266;
	
	Matrix ans_Uranus(3);
	
	ans_Uranus(1) = 1412381021799.02; ans_Uranus(2) = -2511346390560.05; ans_Uranus(3) = -1118104931691.82;
	
	Matrix ans_Neptune(3);
	
	ans_Neptune(1) = 1871263335928.14; ans_Neptune(2) = -3928967868616.56; ans_Neptune(3) = -1655016891088.18;
	
	Matrix ans_Pluto(3);
	
	ans_Pluto(1) = -2171402287905.05; ans_Pluto(2) = -3915426845649.37; ans_Pluto(3) = -552714009658.768;
	
	Matrix ans_Moon(3);
	
	ans_Moon(1) = 89840553.9111792; ans_Moon(2) = -336510464.818452; ans_Moon(3) = -114587871.494221;
	
	Matrix ans_Sun(3);
	
	ans_Sun(1) = 92308639955.5137; ans_Sun(2) = -105367651030.449; ans_Sun(3) = -45683095695.9811;
	
    _assert(m_equals(r_Mercury, ans_Mercury, 1e+10));
	_assert(m_equals(r_Venus, ans_Venus, 1e+10));
	_assert(m_equals(r_Earth, ans_Earth, 1e+10));
	_assert(m_equals(r_Mars, ans_Mars, 1e+10));
	_assert(m_equals(r_Jupiter, ans_Jupiter, 1e+10));
	_assert(m_equals(r_Saturn, ans_Saturn, 1e+10));
	_assert(m_equals(r_Uranus, ans_Uranus, 1e+10));
	_assert(m_equals(r_Neptune, ans_Neptune, 1e+10));
	_assert(m_equals(r_Pluto, ans_Pluto, 1e+10));
	_assert(m_equals(r_Moon, ans_Moon, 1e+10));
	_assert(m_equals(r_Sun, ans_Sun, 1e+10));
	
    return 0;
}

int m_LTC_01(){

    Matrix R = LTC(2,4);

    Matrix ans(3,3);    
	ans(1,1) = -0.909297426825682; ans(1,2) =  -0.416146836547142; ans(1,3) = 0;
	ans(2,1) = -0.314940964313378; ans(2,2) = 0.688158561598754; ans(2,3) = -0.653643620863612;
	ans(3,1) = 0.272011725051612; ans(3,2) =  -0.594356462512304; ans(3,3) = -0.756802495307928;

    _assert(m_equals(R, ans, 1e-10));
	
    return 0;
}

int m_NutMatrix_01(){

    Matrix R = NutMatrix(5);

    Matrix ans(3,3);    
	ans(1,1) = 0.999999999596936; ans(1,2) =  -2.60458970189426e-05; ans(1,3) = -1.13021883531116e-05;
	ans(2,1) = 2.6045466894546e-05; ans(2,2) = 0.999999998936717; ans(2,3) = -3.80552143565138e-05;
	ans(3,1) = 1.13031795232884e-05; ans(3,2) =  3.80549199703872e-05; ans(3,3) = 0.999999999212031;
	
    _assert(m_equals(R, ans, 1e-10));
	
    return 0;
}

int m_PoleMatrix_01(){

    Matrix R = PoleMatrix(3,5);

    Matrix ans(3,3);    
	ans(1,1) = -0.989992496600445; ans(1,2) =  -0.135323401369264; ans(1,3) = 0.04003040989885;
	ans(2,1) = 0; ans(2,2) = 0.283662185463226; ans(2,3) = 0.958924274663138;
	ans(3,1) = -0.141120008059867; ans(3,2) =  0.949327836724532; ans(3,3) = -0.280823435177878;

    _assert(m_equals(R, ans, 1e-10));
	
    return 0;
}

int m_PrecMatrix_01(){

    Matrix R = PrecMatrix(2,4);

    Matrix ans(3,3);    
	ans(1,1) = 0.99999999999911; ans(1,2) =  -1.22341466610288e-06; ans(1,3) = -5.32402961405738e-07;
	ans(2,1) = 1.22341466610288e-06; ans(2,2) = 0.999999999999252; ans(2,3) = -3.25674798695904e-13;
	ans(3,1) = 5.32402961405738e-07; ans(3,2) =  -3.25674792564771e-13; ans(3,3) = 0.999999999999858;

    _assert(m_equals(R, ans, 1e-10));
	
    return 0;
}

int m_gmst_01(){

    double R = gmst(5);

    double ans = 1.05922210457995;

    _assert(fabs(R - ans)< 1e-10);
	
    return 0;
}

int m_gast_01(){

    double R = gast(3);

    double ans = 1.02484149759777;

    _assert(fabs(R - ans)< 1e-10);
	
    return 0;
}

int m_MeasUpdate_01(){

    Matrix x(2);
	
	x(1) = 1; x(2) = 2;
    
    Matrix G(2);
	
	G(1) = 4; G(2) = 5;

    Matrix P(2,2);
	
	P(1,1) = 2; P(1,2) = 6;
    P(2,1) = 4; P(2,2) = 1; 

    auto [K, x1, P1] = MeasUpdate(x, 3, 4, 2, G, P, 2);
	
    Matrix ans_K(2);

	ans_K(1) = 0.145593869731801; ans_K(2) = 0.0804597701149425;
	
    Matrix ans_x(2);

	ans_x(1) = 0.854406130268199; ans_x(2) = 1.91954022988506;

    Matrix ans_P(2,2);

	ans_P(1,1) = -2.07662835249042; ans_P(1,2) = 1.77777777777778;
    ans_P(2,1) = 1.74712643678161; ans_P(2,2) = -1.33333333333333; 

    _assert(m_equals(K, ans_K, 1e-10));
    _assert(m_equals(x1, ans_x, 1e-10));
    _assert(m_equals(P1, ans_P, 1e-10));
	
    return 0;
}

int m_GAccelHarmonic_01() {

	Matrix r(3);
    r(1,1) = 1; r(1,2) = 2; r(1,3) = 0;

    Matrix U(3,3);
    U(1,1) = 4; U(1,2) = 0; U(1,3) = 5;
    U(2,1) = 6; U(2,2) = 7; U(2,3) = 0;
    U(3,1) = 5; U(3,2) = 1; U(3,3) = 0;
    
    Matrix R = G_AccelHarmonic(r, U, 2, 2);

	Matrix ans(3,3);
    ans(1,1) = 1.64484907894009e+21; ans(1,2) = 1.06288145342028e+21; ans(1,3) = 8.37443192155128e+19;
    ans(2,1) = 1.21743011956612e+21; ans(2,2) = 4.50817680497753e+20; ans(2,3) = 6.13698828978822e+19;
    ans(3,1) = 8.19164736011343e+19; ans(3,2) = 6.42949163098142e+19; ans(3,3) = -5.65166746051633e+19;
	
    _assert(m_equals(R, ans, 1e+10));
    
    return 0;
}

int m_GHAMatrix_01() {
    
    Matrix R = GHAMatrix(3);

	Matrix ans(3,3);
    ans(1,1) = 0.519234354564389; ans(1,2) = 0.854631900317384; ans(1,3) = 0;
    ans(2,1) = -0.854631900317384; ans(2,2) = 0.519234354564389; ans(2,3) = 0;
    ans(3,1) = 0; ans(3,2) = 0; ans(3,3) = 1;
    
    _assert(m_equals(R, ans, 1e-10));
    
    return 0;
}

int m_Accel_01() {
    
    Matrix Y = zeros(6);
	
	Y(1) = 1111111; Y(2) = 2222222; Y(3) = 3333333;
	Y(4) = 4444444; Y(5) = 5555555; Y(6) = 6666666;
	
    Matrix R = Accel(777777, Y);
	
	Matrix ans =zeros(6);
    ans(1) = 4444444.0; ans(2) = 5555555.0; ans(3) = 6666666.0;
    ans(4) = -6.2043905967211; ans(5) = -12.0025757141433; ans(6) = -18.5369187965181;
    
    _assert(m_equals(R, ans, 1e-8));
    
    return 0;
}

int m_VarEqn_01() {
    
    Matrix& yPhi = zeros(42);
	
    yPhi(1) = 7101576.98990384;
    yPhi(2) = 1295199.87127754;
    yPhi(3) = 12739.2823333892;
    yPhi(4) = 576.004651193009;
    yPhi(5) = -3084.62203617271;
    yPhi(6) = -6736.02594582756;
    yPhi(7) = 1.00002525535511;
    yPhi(8) = 7.08259815373561e-06;
    yPhi(9) = 1.91608861002907e-07;
    yPhi(10) = 1.01043851887223e-05;
    yPhi(11) = 2.82768336557965e-06;
    yPhi(12) = 6.44131451075285e-08;
    yPhi(13) = 7.08259834024473e-06;
    yPhi(14) = 0.999988040046622;
    yPhi(15) = 3.53015288644891e-08;
    yPhi(16) = 2.82768357826951e-06;
    yPhi(17) = -4.78603729288896e-06;
    yPhi(18) = 1.18527461137171e-08;
    yPhi(19) = 1.9160935046062e-07;
    yPhi(20) = 3.53016114843062e-08;
    yPhi(21) = 0.999986704774626;
    yPhi(22) = 6.44136325079115e-08;
    yPhi(23) = 1.18528331537947e-08;
    yPhi(24) = -5.31820682446032e-06;
    yPhi(25) = 5.00001498082565;
    yPhi(26) = 1.1781862826826e-05;
    yPhi(27) = 2.68389762645616e-07;
    yPhi(28) = 1.00002526606744;
    yPhi(29) = 7.05571100144785e-06;
    yPhi(30) = 1.30455137405173e-07;
    yPhi(31) = 1.17818628919961e-05;
    yPhi(32) = 4.99995293819715;
    yPhi(33) = 4.93630678596596e-08;
    yPhi(34) = 7.05571117883108e-06;
    yPhi(35) = 0.999988029832331;
    yPhi(36) = 2.39618837211068e-08;
    yPhi(37) = 2.68390168073246e-07;
    yPhi(38) = 4.93631303180711e-08;
    yPhi(39) = 4.99995072081276;
    yPhi(40) = 1.30455621823661e-07;
    yPhi(41) = 2.39619698989173e-08;
    yPhi(42) = 0.999986704276552;
    
    Matrix& R = VarEqn(2.5, yPhi);

    Matrix& ans = zeros(42);
	
    ans(1) = 576.004651193009;
    ans(2) = -3084.62203617271;
    ans(3) = -6736.02594582756;
    ans(4) = -7.53466223591456;
    ans(5) = -1.37422019436673;
    ans(6) = -0.0135523187801885;
    ans(7) = 1.01043851887223e-05;
    ans(8) = 2.82768336557965e-06;
    ans(9) = 6.44131451075285e-08;
    ans(10) = 2.02219654159749e-06;
    ans(11) = 5.62315203773639e-07;
    ans(12) = 5.54306320850871e-09;
    ans(13) = 2.82768357826951e-06;
    ans(14) = -4.78603729288896e-06;
    ans(15) = 1.18527461137171e-08;
    ans(16) = 5.62315384036549e-07;
    ans(17) = -9.58426103462621e-07;
    ans(18) = 1.01508481640735e-09;
    ans(19) = 6.44136325079115e-08;
    ans(20) = 1.18528331537947e-08;
    ans(21) = -5.31820682446032e-06;
    ans(22) = 5.54345246738689e-09;
    ans(23) = 1.0151546438867e-09;
    ans(24) = -1.06368579634921e-06;
    ans(25) = 1.00002526606744;
    ans(26) = 7.05571100144785e-06;
    ans(27) = 1.30455137405173e-07;
    ans(28) = 1.01107443600396e-05;
    ans(29) = 2.81153608464331e-06;
    ans(30) = 2.77154087134049e-08;
    ans(31) = 7.05571117883108e-06;
    ans(32) = 0.999988029832331;
    ans(33) = 2.39618837211068e-08;
    ans(34) = 2.81153629664837e-06;
    ans(35) = -4.79215601400399e-06;
    ans(36) = 5.07544126578285e-09;
    ans(37) = 1.30455621823661e-07;
    ans(38) = 2.39619698989173e-08;
    ans(39) = 0.999986704276552;
    ans(40) = 2.77158915838407e-08;
    ans(41) = 5.07552473765848e-09;
    ans(42) = -5.3184472780816e-06;

    _assert(m_equals(R, ans, 1e-10));
    return 0;
}

int m_DEInteg_01() {
    
    int n_eqn = 6;
    double t = 0.0, relerr = 1e-13, abserr = 1e-6, tout = -134.999991953373;

    Matrix Y = zeros(6);

    Y(1) = 6221397.62857869; Y(2) = 2867713.77965738; Y(3) = 3006155.98509949;
    Y(4) = 4645.04725161806; Y(5) = -2752.21591588204; Y(6) = -7507.99940987031;
    
    Matrix& R = DEInteg(Accel, t, tout, relerr, abserr, n_eqn, Y);

    Matrix& ans = zeros(6);

    ans(1) = 5542555.89427451; ans(2) = 3213514.83814162; ans(3) = 3990892.92789074;
    ans(4) = 5394.06894044389; ans(5) = -2365.21290574021; ans(6) = -7061.8448137347;
    
    _assert(m_equals(R, ans, 1e-5));
    
    return 0;
}

int all_tests()
{	
    _verify(m_sum_01);
	_verify(m_sum_02);
    _verify(m_sub_01);
    _verify(m_sub_02);
    _verify(m_zeros_01);
    _verify(m_prod_01);
    _verify(m_prod_02);
    _verify(m_div_01);
    _verify(m_div_02);
    _verify(m_equ_01);
	_verify(m_inv_01);
    _verify(m_trans_01);
    _verify(m_eye_01);
    _verify(m_norm_01);
    _verify(m_dot_01);
    _verify(m_cross_01);
    _verify(m_assign_column_01);
    _verify(m_assign_row_01);
    _verify(m_extract_column_01);
    _verify(m_extract_row_01);
    _verify(m_union_vector_01);
    _verify(m_extract_vector_01); //22 test
	
	_verify(m_r_x_01);
	_verify(m_r_y_01);
	_verify(m_r_z_01);
	_verify(m_AccelPointMass_01);
	_verify(m_Cheb_01);
	_verify(m_EccAnom_01);
	_verify(m_Frac_01);
	_verify(m_MeanObli_01);
	_verify(m_Mjday_01);
	_verify(m_Mjday_02);
	_verify(m_MjdayTDB_01);
	_verify(m_Position_01);
	_verify(m_sign__01);
	_verify(m_sign__02);
	_verify(m_timediff_01);
	_verify(m_AzElPa_01);
	_verify(m_IERS_01);
	_verify(m_Legendre_01);
	_verify(m_NutAngles_01);
	_verify(m_TimeUpdate_01);
	_verify(m_TimeUpdate_02); //43 test

	_verify(m_AccelHarmonic_01);
	_verify(m_EqnEquinox_01);
	_verify(m_JPL_Eph_DE430_01);
	_verify(m_LTC_01);
	_verify(m_NutMatrix_01);
	_verify(m_PoleMatrix_01);
	_verify(m_PrecMatrix_01);
	_verify(m_gmst_01); //51 test
    
	_verify(m_gast_01);
	_verify(m_MeasUpdate_01);
	_verify(m_GAccelHarmonic_01);
	_verify(m_GHAMatrix_01);
	_verify(m_Accel_01);
	_verify(m_VarEqn_01); //57 test

	_verify(m_DEInteg_01); //58 test

    return 0;
}


int main()
{
    eop19620101(21413);
    GGM03S(181);
    DE430Coeff(2285, 1020);
	initAuxParam();

    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
