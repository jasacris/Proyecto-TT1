#include "..\include\matrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
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
	A(1) = 1; A(2) = 2; A(3) = 3;
	
	double ans = 14;
	
	double R = norm(A);
    
    _assert(ans == R);
    
    return 0;
}

int m_dot_01() {
	
	Matrix A(3);
	A(1) = 1; A(2) = 2; A(3) = 3;
	
	Matrix B(3);
	B(1) = 3; B(2) = 2; B(3) = 1;
	
	double ans = 10;
	
	double R = dot(A, B);
    
    _assert(ans == R);
    
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

int all_tests()
{
    _verify(m_sum_01);
	_verify(m_sum_02);
    _verify(m_sub_01);
    _verify(m_sub_02);
    _verify(m_zeros_01);
    _verify(m_prod_01);
    _verify(m_prod_02);
    //_verify(m_div_01);
    _verify(m_div_02);
    _verify(m_equ_01);
	//_verify(m_inv_01);
    _verify(m_trans_01);
    _verify(m_eye_01);
    _verify(m_norm_01);
    _verify(m_dot_01);
    _verify(m_cross_01);
	
	_verify(m_r_x_01);
	_verify(m_r_y_01);
	_verify(m_r_z_01);

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
