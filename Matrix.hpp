// Authors: Nicholas Barton

#ifndef MATRIXHEADERDEF__
#define MATRIXHEADERDEF__

#include <string>
#include "Vector.hpp"

using namespace std;

class Matrix
{
	friend class ::Vector; // make Vector a friend class of Matrix
private:
	// attributes
	int mRows, mCols;
	double** mData;

	// functions only used internally, not for user purposes
	friend int lu_det(Matrix& A, Matrix& U); // LU w/ pivotting but returns the number of permutes
	friend Vector arnoldi(const Matrix& A, Matrix& Q, int k); // does kth iteration of arnoldi
public:
	// constructors
	// no ``default'' constructor! must specify something about the matrix
	// construct matrix of given size
	Matrix(int numentries);
	Matrix(int numrows, int numcols);

	// construct matrix of given size initialised to a certain value
	Matrix(int numentries, double initVal);
	Matrix(int numrows, int numcols, double initVal);

	// special matrix constructors
	Matrix(int numentries, string type);
	Matrix(int numrows, int numcols, string type);
	Matrix(int numentries, string type, double initVal);
	Matrix(int numrows, int numcols, string type, double initVal);

	// overridden copy constructor
	Matrix(const Matrix& A);

	// destructor
	~Matrix();

	// indexing
	double& operator()(int i, int j); // ij index style
	double& operator()(int i); // pretend that the matrix is a really long vector

	// various slices of the Matrix for utility
	friend Vector col_slice(const Matrix& A, int col);
	friend Vector row_slice(const Matrix& A, int row);
	friend Matrix slice(const Matrix& A, int row_start, int row_end, int col_start, int col_end);

	// eye and diag as functions and not constructors
	friend Matrix eye(int numentries); // square identity
	friend Matrix eye(int numrows, int numcols); // non-square
	friend Matrix diag(int numentries, double initVal); // square
	friend Matrix diag(int numrows, int numcols, double initVal); // non-square
	friend Matrix diag(const Vector& v); // square matrix with given vector on diagonals

	// overload unary minus operator
	friend Matrix operator-(const Matrix& A);

	// overload binary arithmetic operators
	friend Matrix operator+(const Matrix& A, const Matrix& B);
	friend Matrix operator-(const Matrix& A, const Matrix& B);
	friend Matrix operator*(const Matrix& A, const Matrix& B); // mat-mat mult
	friend Matrix operator*(const Matrix& A, const double& c); // mat-sca mult
	friend Matrix operator*(const double& c, const Matrix& A); // sca-mat mult
	friend Vector operator*(const Matrix& A, const Vector& x); // mat-vec mult
	friend Vector operator*(const Vector& x, const Matrix& A); // vec-mat mult
	friend Matrix operator/(const Matrix& A, const double& c);

	// overload = (assignment) operator
	Matrix& operator=(const Matrix& A);
	

	// to get the properties of a Matrix
	friend int length(const Matrix& A); // gives the total number of entries
	int length() const;
	friend int width(const Matrix& A); // gives the number of cols
	int width() const;
	friend int height(const Matrix& A); // gives the number of rows
	int height() const;
	friend int* size(const Matrix& A); // gives an array of the size of the Matrix (rows, cols)
	int* size() const;
	friend bool is_symmetric(const Matrix& A); // returns true is A is symmetric
	bool is_symmetric() const;
	friend bool is_singular(const Matrix& A); // checks if the determinant is 0
	bool is_singular() const;

	// output functions
	friend ostream& operator<<(ostream& output, const Matrix& A);
	friend string to_string(const Matrix& A);

	// outer product of two Vectors is a Matrix
	friend Matrix outer(const Vector& v1, const Vector& v2);

	// various properties of the Matrix
	friend double norm(const Matrix& A, int p); 
	friend double norm(const Matrix& A, string type); // for infinity & Frobenius norms
	double norm(int p = 1) const; // only 1 norm is implemented, so we make that the default
	double norm(string) const;
	friend double det(const Matrix& A);
	double det() const;

	// friendly transpose function returns a new Matrix
	friend Matrix transpose(const Matrix& A);

	// in-place transpose
	Matrix& transpose();

	// decompositions and factorisations
	// output L, input A unchanged, U to be filled in -- NO PIVOTING
	friend Matrix lu(const Matrix& A, Matrix& U);
	// output L, input A unchanged, P, U to be filled in
	friend Matrix lu(Matrix& A, Matrix& U, Matrix& P);
	friend Matrix qr(Matrix& A, bool explicit_Q); // A overwritten to R, output is W or Q

	// linear solvers return the solution of Ax = b
	// triangular linear solvers
	friend Vector backsub(const Matrix& A, const Vector& b);
	friend Vector forwardsub(const Matrix& A, const Vector& b);

	// iterative linear solvers
	friend Vector cg(const Matrix& A, const Vector& b, bool verbose, double TOL, int maxit);
	friend Vector cg(const Matrix& A, const Vector& b, const Vector& x0, bool verbose, double TOL, int maxit);
	friend Vector gmres(const Matrix& A, const Vector& b, bool verbose, double TOL, int maxit);
	friend Vector gmres(const Matrix& A, const Vector& b, const Vector& x0, bool verbose, double TOL, int maxit);

};

// all the function prototypes are given below, no need to include member methods
int length(const Matrix& A);
int width(const Matrix& A);
int height(const Matrix& A);
int* size(const Matrix& A);
Vector col_slice(const Matrix& A, int col); // pick out columns of A
Vector row_slice(const Matrix& A, int row); // pick out rows of A
Matrix slice(int row_start, int row_end, int col_start, int col_end); // to slice

string to_string(const Matrix& A);

Matrix operator+(const Matrix& A);
Matrix operator-(const Matrix& B);
Matrix operator-(const Matrix& A, const Matrix& B); 
Matrix operator/(const Matrix& A, const double& c); 
Matrix operator*(const Matrix& A, const Matrix& B); // mat-mat
Matrix operator*(const Matrix& A, const double& c); // mat-sca
Matrix operator*(const double& c, const Matrix& A); // sca-mat
Vector operator*(const Matrix& A, const Vector& x); // mat-vec
Vector operator*(const Vector& x, const Matrix& A); // vec-mat
Matrix outer(const Vector& v1, const Vector& v2); // vec-vec outer prod

int* size(const Matrix& A);
int length(const Matrix& A);
int width(const Matrix& A);
int height(const Matrix& A);
bool is_symmetric(const Matrix& A);
bool is_singular(const Matrix& A);

Matrix eye(int numentries); // square identity
Matrix eye(int numrows, int numcols); // non-square
Matrix diag(int numentries, double initVal); // square
Matrix diag(int numrows, int numcols, double initVal); // non-square
Matrix diag(const Vector& v); // square matrix w vector v on diagonal

double norm(const Matrix& A, int p = 1);
double norm(const Matrix& A, string type); // for infinity norm
double det(const Matrix& A);
Matrix transpose(const Matrix& A); // returns a new Matrix A^T

// factorisations and decompositions
Matrix lu(const Matrix& A, Matrix& U); // LU w/o pivotting
int lu_det(Matrix& A, Matrix& U); // LU w/ pivotting but returns the number of permutes
Matrix lu(Matrix& A, Matrix& U, Matrix& P); // output L, input A unchanged, P, U to be filled in
Matrix qr(Matrix& A, bool explicit_Q = 0); // A overwritten to R, output is Q

// functions for solution of linear systems
Vector backsub(const Matrix& A, const Vector& b);
Vector forwardsub(const Matrix& A, const Vector& b);

// iterative methods
Vector arnoldi(const Matrix& A, Matrix& Q, int k);
Vector cg(const Matrix& A, const Vector& b, bool verbose = 0, double TOL = 1E-15, int maxit = 256);
Vector cg(const Matrix& A, const Vector& b, const Vector& x0, bool verbose = 0, double TOL = 1E-15, int maxit = 256);
Vector gmres(const Matrix& A, const Vector& b, bool verbose = 0, double TOL = 1E-15, int maxit = 256);
Vector gmres(const Matrix& A, const Vector& b, const Vector& x0, bool verbose = 0, double TOL = 1E-15, int maxit = 256);

#endif