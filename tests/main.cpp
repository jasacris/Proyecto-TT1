#include "..\include\matrix.h"
#include <iostream>

using namespace std;

int main() {
    Matrix v(3);
	v(2) = 5;
	cout << v;

    return 0;
}