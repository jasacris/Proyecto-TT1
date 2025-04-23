#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix {
public:
    int n_row, n_column;
	double **data;

    // Parameterized constructor
	Matrix();
	Matrix(const int v_size);
    Matrix(const int n_row, const int n_column);
	
	// Member operators
	double& operator () (const int n);
	double& operator () (const int row, const int column);
	Matrix& operator + (Matrix &m);
	Matrix& operator - (Matrix &m);
	Matrix& operator * (Matrix &m);
	Matrix& operator / (Matrix &m);
	Matrix& operator = (Matrix &m);
	Matrix& operator + (const double n);
	Matrix& operator - (const double n);
	Matrix& operator * (const double n);
	Matrix& operator / (const double n);
	
	// Non-member operators
	friend ostream& operator << (ostream &o, Matrix &m);
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int v_size);
Matrix& zeros(const int n_row, const int n_column);
Matrix& inv (Matrix &m);
Matrix& transponse (Matrix &m);
Matrix& eye(const int n);
double norm (Matrix &m);
double dot (Matrix &m1, Matrix &m2);
Matrix& cross (Matrix &m1, Matrix &m2);

#endif