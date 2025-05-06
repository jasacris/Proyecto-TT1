#include "..\include\matrix.hpp"
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
/*
int m_MeanObli_01() {

	double x = 5;

	double ans = 0.4094;

	double R = MeanObliquity(x);
    
    _assert(fabs(R - ans) < 1e-10);
    
    return 0;
}*/

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
/*
int m_Position_01() {
	
	double lon = 0.5;
	double lat = 0;
	double h = 1;

	Matrix C(3);
	C(1) = -1438078.785611559; C(2) =  -2239675.009373783; C(3) = 5776810.445003163;
	
	Matrix R = Position(lon, lat, h);
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}*/

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

    double ans_Az = 1.32857846;
    double ans_El = 14;
    Matrix ans_dAds = 17;
    Matrix ans_dEds = 37.184;

    _assert(fabs(Az - ans_Az)< 1e-10);
    _assert(fabs(El - ans_El)< 1e-10);
    _assert(m_equals(dAds, ans_dAds, 1e-10));
    _assert(m_equals(dEds, ans_dEds, 1e-10));
	
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
    _verify(m_extract_vector_01);
	
	_verify(m_r_x_01);
	_verify(m_r_y_01);
	_verify(m_r_z_01);
	_verify(m_AccelPointMass_01);
	_verify(m_Cheb_01);
	_verify(m_EccAnom_01);
	_verify(m_Frac_01);
	//_verify(m_MeanObli_01);
	_verify(m_Mjday_01);
	_verify(m_Mjday_02);
	_verify(m_MjdayTDB_01);
	//_verify(m_Position_01);
	_verify(m_sign__01);
	_verify(m_sign__02);
	_verify(m_timediff_01);
	_verify(m_AzElPa_01);

    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
