#include "..\include\matrix.hpp"

Matrix::Matrix() {
	this->n_row = 0;
	this->n_column = 0;
	this->data = nullptr;
}

Matrix::Matrix(const int v_size) {
    if (v_size <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = 1;
	this->n_column = v_size;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	this->data[0] = (double *) calloc(v_size,sizeof(double));
}

Matrix::Matrix(const int n_row, const int n_column) {
    if (n_row <= 0 || n_column <= 0) {
		cout << "Matrix create: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	this->n_row = n_row;
	this->n_column = n_column;
	this->data = (double **) malloc(n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(n_column*sizeof(double));
	}
}

double& Matrix::operator () (const int n) {
	if (n <= 0 || n > this->n_row*this->n_column) {
		cout << "Matrix get1: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[(n - 1)/this->n_column][(n - 1)%this->n_column];
}

double& Matrix::operator () (const int row, const int column) {
	if (row <= 0 || row > this->n_row || column <= 0 || column > this->n_column) {
		cout << "Matrix get2: error in row/column\n";
        exit(EXIT_FAILURE);
	}
	
	return this->data[row - 1][column - 1];
}

Matrix& Matrix::operator + (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sum1: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + m(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator + (const double n) {
	if (this->n_row < 0|| this->n_column < 0) {
		cout << "Matrix sum2: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + n;
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator - (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub1: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator - (const double n) {
	if (this->n_row < 0|| this->n_column < 0) {
		cout << "Matrix sub2: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - n;
		}
	}
	
	return *m_aux;
}

ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

Matrix& zeros(const int v_size) {
	Matrix *m_aux = new Matrix(1, v_size);
	
	for(int i = 1; i <= v_size; i++) {
		(*m_aux)(1, i) = 0;
	}
	
	return (*m_aux);
}

Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

Matrix& Matrix::operator * (Matrix &m){
	if (this->n_row != m.n_column) {
		cout << "Matrix prod1: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, m.n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= m.n_column; j++) {
			double suma = 0;
			for(int k = 1; k <= this->n_column; k++){
				suma += (*this)(i, k)*m(k, j);
			}
			(*m_aux)(i,j) = suma;
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator * (const double n){
	if (this->n_row < 0|| this->n_column < 0) {
		cout << "Matrix prod2: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) * n;
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator / (Matrix &m){
	if (this->n_column != m.n_column) {
		cout << "Matrix div1: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix((*this) * inv(m));
	
	return *m_aux;
}

Matrix& Matrix::operator / (const double n){
	if (this->n_column < 0) {
		cout << "Matrix div2: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) / n;
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator = (Matrix &m){	
	this->n_row = m.n_row;
	this->n_column = m.n_column;
	this->data = (double **) malloc(m.n_row*sizeof(double *));
	
    if (this->data == NULL) {
		cout << "Matrix create: error in data\n";
        exit(EXIT_FAILURE);
	}
	
	for(int i = 0; i < n_row; i++) {
		this->data[i] = (double *) malloc(m.n_column*sizeof(double));
	}
	
	for(int i = 1; i <= m.n_row; i++){
		for(int j = 1; j <= m.n_column; j++){
			(*this)(i,j) = m(i,j);
		}
	}
	
	return (*this);
}

Matrix& inv (Matrix &m){	
	if (m.n_row != m.n_column) {
		cout << "Matrix inv: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	int n = m.n_row; 
	
	Matrix *m_aux = new Matrix(n, n);
    Matrix *identity = new Matrix(eye(n));

    // Copia m en m_aux
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*m_aux)(i, j) = m(i, j);
        }
    }

    // Realiza Gauss-Jordan
    for (int i = 1; i <= n; i++) {
		int piv = i;
        double diag = (*m_aux)(i, i);

		for (int j = i + 1; j <= n; j++) {
            if (abs((*m_aux)(j, i)) > abs(diag)) {
                diag = (*m_aux)(j, i);
                piv = j;
            }
        }

        if (diag == 0) {
            cout << "Matrix inv: error, la matriz es singular\n";
            exit(EXIT_FAILURE);
        }

		if(piv != i){
			Matrix& mAux = (*m_aux).extract_row(piv);
			(*m_aux).assign_row((*m_aux).extract_row(i), piv);
			(*m_aux).assign_row(mAux, i);
			
			mAux = (*identity).extract_row(piv);
			(*identity).assign_row((*identity).extract_row(i), piv);
			(*identity).assign_row(mAux, i);
		}

        for (int j = 1; j <= n; j++) {
            (*m_aux)(i, j) /= diag;
            (*identity)(i, j) /= diag;
        }

        for (int k = 1; k <= n; k++) {
            if (k != i) {
                double factor = (*m_aux)(k, i);
                for (int j = 1; j <= n; j++) {
                    (*m_aux)(k, j) -= factor * (*m_aux)(i, j);
                    (*identity)(k, j) -= factor * (*identity)(i, j);
                }
            }
        }
    }

    return *identity;
}


Matrix& transponse (Matrix &m){	
	Matrix *m_aux = new Matrix(m.n_column, m.n_row);
	
    for(int i = 1; i <= m.n_row; i++) {
        for(int j = 1; j <= m.n_column; j++) {
			(*m_aux)(j,i) = m(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& eye (const int n){	
	Matrix *m_aux = new Matrix(n, n);
	
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
			if(i != j){
				(*m_aux)(i, j) = 0;
			}else{
				(*m_aux)(i, j) = 1;
			}
		}
	}
	
	return *m_aux;
}

double norm (Matrix &m){	
	if (m.n_row != 1) {
		cout << "Matrix norm: error in n_row\n";
        exit(EXIT_FAILURE);
	}
	double ans = 0;
	
    for(int i = 1; i <= m.n_column; i++){
		ans += pow(m(i),2);
	}
	
	return sqrt(ans);
}

double dot (Matrix &m1, Matrix &m2){	
	if (m1.n_row != 1 || m2.n_row != 1 || m1.n_column != m2.n_column) {
		cout << "Matrix dot: error in n_row\n";
        exit(EXIT_FAILURE);
	}
	double ans = 0;
	
    for(int i=1; i <= m1.n_column; i++){
		ans += m1(i) * m2(i);
	}
	
	return ans;
}

Matrix& cross (Matrix &m1, Matrix &m2){	
	if (m1.n_row != 1 || m2.n_row != 1 || m1.n_column != 3 || m2.n_column != 3) {
		cout << "Matrix cross: error in n_row\n";
        exit(EXIT_FAILURE);
	}
	Matrix *ans = new Matrix(3);
	
    (*ans)(1) = m1(2) * m2(3) - m1(3) * m2(2);
	(*ans)(2) = m1(3) * m2(1) - m1(1) * m2(3);
	(*ans)(3) = m1(1) * m2(2) - m1(2) * m2(1);
	
	return *ans;
}

Matrix& Matrix::assign_column(Matrix &m, const int col){
	if (col < 1 || col > this->n_column || m.n_row != 1){
		cout << "Matrix assign_column: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}		

	for (int i = 1; i <= this->n_row; i++){
		(*this)(i, col) = m(i);
	}

	return *this;
}

Matrix& Matrix::assign_row(Matrix &m, const int row){
	if (row < 1 || row > this->n_row || m.n_row != 1){
		cout << "Matrix assign_row: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	for (int j = 1; j <= this->n_column; j++){
		(*this)(row, j) = m(j);
	}

	return *this;
}

Matrix& Matrix::extract_column(const int col){
	if (col < 1 || col > this->n_column){
		cout << "Matrix extract_column: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m = new Matrix(this->n_row);

	for (int i = 1; i <= this->n_row; i++){
		(*m)(i) = (*this)(i, col);
	}

	return *m;
}

Matrix& Matrix::extract_row(const int row){
	if (row < 1 || row > this->n_row){
		cout << "Matrix extract_row: error in n_row/n_column\n";
		exit(EXIT_FAILURE);
	}

	Matrix *m = new Matrix(this->n_column);

	for (int j = 1; j <= this->n_column; j++){
		(*m)(1, j) = (*this)(row, j);
	}

	return *m;
}

Matrix& Matrix::union_vector(Matrix &m){
	if (this->n_row != 1 || m.n_row != 1){
		cout << "Matrix union_vector: error in n_row\n";
		exit(EXIT_FAILURE);
	}

	int tamano = this->n_column + m.n_column;
	Matrix *m_aux = new Matrix(tamano);

	int aux = 1;
	for (int i = 1; i <= this->n_column; i++){
		(*m_aux)(1, aux) = (*this)(i);
		aux++;
	}
	for (int i = 1; i <= m.n_column; i++){
		(*m_aux)(1, aux) = m(i);
		aux++;
	}

	return *m_aux;
}

Matrix& Matrix::extract_vector(const int from, const int to){
	if (this->n_row != 1 || from <= 0 || to > this->n_column || from > to){
		cout << "Matrix extract_vector: error in n_row\n";
		exit(EXIT_FAILURE);
	}

	int sub_size = to - from + 1;
	Matrix *m = new Matrix(sub_size);

	int index = 1;
	for (int i = from; i <= to; i++){
		(*m)(index) = (*this)(i);
		index++;
	}

	return *m;
}