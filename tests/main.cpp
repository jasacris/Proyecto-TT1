#include "..\include\matrix.h"
#include <iostream>

using namespace std;

int main() {
    atrix A(2, 3);
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

	cout << A << endl;
	cout << B << endl;
	cout << C << endl;
	cout << R << endl;

    return 0;
}