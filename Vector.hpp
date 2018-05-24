// Authors: Joe Pitt-Francis, Nicholas Barton

#ifndef VECTORDEF
#define VECTORDEF



//  **********************
//  *  Class of vectors  *
//  **********************


//  Class written in such a way that code similar to Matlab
//  code may be written


#include <cmath>
#include <string>
#include "Exception.hpp"//  This class throws errors using the class "error"

using namespace std;

class Matrix; // forward-declare Matrix

class Vector
{
   friend class ::Matrix; // make Matrix a friend class of Vector
private:
   // member variables
   double* mData;   // data stored in vector
   int mSize;      // size of vector

   // functions only used internally, not for user purposes
   friend Vector arnoldi(const Matrix& A, Matrix& Q, int k); // does kth iteration of arnoldi
public:
   // constructors
   // No default constructor
   // overridden copy constructor
   Vector(const Vector& v1);
   // construct vector of given size
   Vector(int sizeVal);

   // destructor
   ~Vector();

   friend Matrix diag(const Vector& v); // diag constructor in Matrix


   // All "friend" external operators and functions are declared as friend inside the class (here)
   // but their actual prototype definitions occur outside the class.
   // Binary operators
   friend Vector operator+(const Vector& v1, const Vector& v2);
   friend Vector operator-(const Vector& v1, const Vector& v2);
   friend double operator*(const Vector& v1, const Vector& v2);
   friend Vector operator*(const Vector& v, const double& a);
   friend Vector operator*(const double& a, const Vector& v);
   friend Vector operator/(const Vector& v, const double& a);
   friend Vector operator*(const Matrix& A, const Vector& x); // friend func with Matrix
   friend Vector operator*(const Vector& x, const Matrix& A); //friend func with Matrix
   // Unary operator
   friend Vector operator-(const Vector& v);

   //other operators
   //assignment
   Vector& operator=(const Vector& v);
   //indexing
   double& operator()(int i);
   friend Vector slice(const Vector& v, int start, int end);
   friend Vector col_slice(const Matrix& A, int col);
   friend Vector row_slice(const Matrix& A, int row);
   //output
   friend ostream& operator<<(ostream& output, const Vector& v);
   // overload to_string for Vector
   friend string to_string(const Vector& v);

   // outer product e.g. column vec * row vec = matrix
   friend Matrix outer(const Vector& v1, const Vector& v2);

   //norm (as a member method)
   double norm(int p = 2) const; // const means the instance is unchanged
   // functions that are friends
   friend double norm(Vector& v, int p);
   friend int length(const Vector& v);

   // forward declare some functions that will be used in Matrix but need access to Vector
   // functions useful for solving linear systems
   friend Vector backsub(const Matrix& A, const Vector& b);
   friend Vector forwardsub(const Matrix& A, const Vector& b);
   friend Matrix qr(Matrix& A, bool explicit_Q);

   // iterative methods
   friend Vector cg(const Matrix& A, const Vector& b, bool verbose, double TOL, int maxit);
   friend Vector cg(const Matrix& A, const Vector& b, const Vector& x0, bool verbose, double TOL, int maxit);
   friend Vector gmres(const Matrix& A, const Vector& b, bool verbose, double TOL, int maxit);
   friend Vector gmres(const Matrix& A, const Vector& b, const Vector& x0, bool verbose, double TOL, int maxit);
};


// All "friend" external operators and functions are declared as friend inside the class
// but their actual prototype definitions occur outside the class (here).
// Binary operators
Vector operator+(const Vector& v1, const Vector& v2);
Vector operator-(const Vector& v1, const Vector& v2);
double operator*(const Vector& v1, const Vector& v2);
Vector operator*(const Vector& v, const double& a);
Vector operator*(const double& a, const Vector& v);
Vector operator/(const Vector& v, const double& a);
// Unary operator
Vector operator-(const Vector& v);

// function prototypes
double norm(Vector& v, int p = 2);
// Prototype signature of length() friend function
int length(const Vector& v);
string to_string(const Vector& v); // overload to_string
Vector col_slice(const Matrix& A, int col);
Vector row_slice(const Matrix& A, int row);
Vector slice(const Vector& v, int start, int end);

// friend functions with Matrix
Vector operator*(const Matrix& A, const Vector& x); // mat-vec
Vector operator*(const Vector& x, const Matrix& A); // vec-mat
Matrix diag(const Vector& v); // diag construction from vector
Matrix outer(const Vector& v1, const Vector& v2); // vec-vec outer prod

// functions for solving linear systems
Vector backsub(const Matrix& A, const Vector& b);
Vector forwardsub(const Matrix& A, const Vector& b);
Matrix qr(Matrix& A, bool explicit_Q);

// iterative methods
Vector arnoldi(const Matrix& A, Matrix& Q, int k);
Vector cg(const Matrix& A, const Vector& b, bool verbose, double TOL, int maxit);
Vector cg(const Matrix& A, const Vector& b, const Vector& x0, bool verbose, double TOL, int maxit);
Vector gmres(const Matrix& A, const Vector& b, bool verbose, double TOL, int maxit);
Vector gmres(const Matrix& A, const Vector& b, const Vector& x0, bool verbose, double TOL, int maxit);

#endif
