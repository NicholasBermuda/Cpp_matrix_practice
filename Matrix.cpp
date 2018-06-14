// Authors: Nicholas Barton

#include <iostream>
#include <cmath>
#include <string>
#include "Matrix.hpp"
#include "Vector.hpp"
#include "Exception.hpp"

using namespace std;


// constructor allocates memory for Matrix of given size
// initialises the matrix to have zeroes for all entries
Matrix::Matrix(int size)
{
	// set the attributes
	mRows = size;
	mCols = size;

	// allocate memory and fill it in
	mData = new double*[size];

	for (int i = 0; i < size; i++)
	{
		mData[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			mData[i][j] = 0.0;
		}
	}
}


// allocates memory for Matrix of given size
// initialise data to zero for all entries
Matrix::Matrix(int rows, int cols)
{
	// set the attributes
	mRows = rows;
	mCols = cols;

	// allocate memory and fill it in
	mData = new double*[rows];

	for (int r = 0; r < rows; r++)
	{
		mData[r] = new double[cols];
		for (int c = 0; c < cols; c++)
		{
			mData[r][c] = 0.0;
		}
	}
}


// constructor allocates memory for matrix of given size
// initialises the matrix to have given initVal for all entries
Matrix::Matrix(int rows, int cols, double initVal)
{
	mRows = rows;
	mCols = cols;

	// allocate memory and fill it in
	mData = new double*[rows];
	
	for (int r = 0; r < rows; r++)
	{
		mData[r] = new double[cols];
		for (int c = 0; c < cols; c++)
		{
			mData[r][c] = initVal;
		}
	}
}

Matrix::Matrix(int size, double initVal)
{
	mRows = size;
	mCols = size;

	// allocate memory and fill it in
	mData = new double*[size];
	for (int i = 0; i < size; i++)
	{
		mData[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			mData[i][j] = initVal;
		}
	}
}


// extensible constructor for special Matrix types
Matrix::Matrix(int size, string type)
{
	
	// allocate memory
	mRows = size;
	mCols = size;
	mData = new double*[size];
	for (int i = 0; i < size; i++)
	{
		mData[i] = new double[size];
	}
	
	// parameter checks
	if (type != "eye")
	{
		throw Exception("invalid parameter",
			"Attempt to crate invalid special matrix. Currently supported: eye, diag.");
	}

	if (type == "eye")
	{
		for (int i = 0; i < size; i++)
			{
				for (int j = 0; j < size; j++)
				{
					if (i==j)
					{
						mData[i][j] = 1.0;
					}
					else
					{
						mData[i][j] = 0.0;
					}
				}
			}
	}
}

Matrix::Matrix(int rows, int cols, string type)
{
	// allocate memory
	mRows = rows;
	mCols = cols;
	mData = new double*[rows];

	for (int r = 0; r < rows; r++)
	{
		mData[r] = new double[cols];
	}
	
	// parameter checks
	if (type != "eye")
	{
		throw Exception("invalid parameter",
			"Attempt to crate invalid special matrix. Currently supported: eye, diag.");
	}

	if (type == "eye")
	{
		for (int r = 0; r < rows; r++)
		{
			for (int c = 0; c < cols; c++)
			{
				if (r==c)
				{
					mData[r][c] = 1.0;
				} // end if
			} // end for
		} // end for
	} // end if
} // end constructor


Matrix::Matrix(int rows, int cols, string type, double initVal)
{
	// allocate memory
	mRows = rows;
	mCols = cols;
	mData = new double*[rows];

	for (int r = 0; r < rows; r++)
	{
		mData[r] = new double[cols];
	}
	
	// parameter checks
	if (type != "diag")
	{
		throw Exception("invalid parameter",
			"Attempt to crate invalid special matrix. Currently supported: eye, diag.");
	}
	
	if (type == "diag")
	{
		for (int r = 0; r < rows; r++)
		{
			for (int c = 0; c < cols; c++)
			{
				if (r==c)
				{
					mData[r][c] = initVal;
				}
				else
				{
					mData[r][c] = 0.0;
				} // end if
			} // end for
		} // end for
	} // end if
} // end constructor


// 
Matrix::Matrix(int size, string type, double initVal)
{
	// allocate memory
	mRows = size;
	mCols = size;
	mData = new double*[size];
	for (int i = 0; i < size; i++)
	{
		mData[i] = new double[size];
	}
	
	// parameter checks
	if (type != "diag")
	{
		throw Exception("invalid parameter",
			"Attempt to crate invalid special matrix. Currently supported: eye, diag.");
	}

	if (type == "diag")
	{
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (i==j)
				{
					mData[i][j] = initVal;
				}
				else
				{
					mData[i][j] = 0.0;
				} // end if
			} // end for
		} // end for
	} // end if
} // end constructor


// copy constructor, create a new Matrix with the same vals as given Matrix
Matrix::Matrix(const Matrix& A)
{
	// allocate memory
	mRows = A.mRows;
	mCols = A.mCols;
	mData = new double*[mRows];
	for (int r = 0; r < mRows; r++)
	{
		mData[r] = new double[mCols];
	}

	// copy entries
	for (int row = 0; row < mRows; row++)
	{
		for (int col = 0; col < mCols; col++)
		{
			mData[row][col] = A.mData[row][col];
		}
	}
}


// destructor deletes dynamically allocated memory
Matrix::~Matrix()
{
	for (int r = 0; r < mRows; r++)
	{
		delete[] mData[r];
	}
	delete[] mData;
}


// indexing with (), given square array abstraction
double& Matrix::operator()(int i, int j)
{
	// dimension checking
	if (i <= 0 || j <= 0 || i > mRows || j > mCols)
	{
		throw Exception("index range","Attempt to index Matrix out of range.");
	}
	return mData[i-1][j-1];
}


// indexing with (), given long list of stacked columns like MATLAB
double& Matrix::operator()(int i)
{
	// start at the beginning
	i -= 1;

	// dimension checking
	if (i < 0 || i > mRows*mCols)
	{
		throw Exception("index range","Attempt to index Matrix out of range.");
	}

	// find the right entry
	int row, col;
	row = i % mCols;
	col = i / mCols;
	return mData[row][col];
}


// get rows and columns of A as Vectors
Vector col_slice(const Matrix& A, int col)
{
	Vector output(A.mRows);

	for (int i = 0; i < A.mRows; ++i)
	{
		output.mData[i] = A.mData[i][col - 1];
	}

	return output;
}


Vector row_slice(const Matrix& A, int row)
{
	Vector output(A.mCols);

	for (int i = 0; i < A.mCols; ++i)
	{
		output.mData[i] = A.mData[row - 1][i];
	}

	return output;
}


// get sub-matrices: inclusive of row_/col_ start/end
Matrix slice(const Matrix& A, int row_start, int row_end, int col_start, int col_end)
{
	// allocate the right size
	Matrix output(row_end - row_start + 1, col_end - col_start + 1);

	for (int i = 0; i < row_end-row_start+1; i++)
	{
		for (int j = 0; j < col_end-col_start+1; j++)
		{
			output.mData[i][j] = A.mData[row_start+i-1][col_start+j-1];
		}
	}

	return output;
}


// assigment operator
Matrix& Matrix::operator=(const Matrix& A)
{

	if (mRows != A.mRows || mCols != A.mCols)
	{
		throw Exception("dimension mismatch",
			"Attempt to assign to a Matrix with different dimensions");
	}

	for (int row = 0; row < mRows; row++)
	{
		for (int col = 0; col < mCols; col++)
		{
			mData[row][col] = A.mData[row][col];
		}
	}

	return *this;
}


// output for printing with <<
ostream& operator<<(ostream& output, const Matrix& A)
{
	output << "(";
	for (int row = 0; row < A.mRows; row++)
	{
		output << "(";
		for (int col = 0; col < A.mCols; col++)
		{
			output << A.mData[row][col];
			if (col < A.mCols - 1)
			{
				output << ", ";
			}
			else
			{
				output << ")";
			}
		}
		if (row < A.mRows - 1)
		{
			output << ", ";
		}
		else
		{
			output << ")";
		}
	}
	output << ")";

	return output;
}


// overload to_string for Matrix --  good for error checks!
string to_string(const Matrix& A)
{
	string output = "(";

	for (int row = 0; row < A.mRows; row++)
	{
		output += "(";
		for (int col = 0; col < A.mCols; col++)
		{
			output += to_string(A.mData[row][col]);
			if (col < A.mCols - 1)
			{
				output += ", ";
			}
			else
			{
				output += ")";
			}
		}
		if (row < A.mRows - 1)
		{
			output += ", ";
		}
		else
		{
			output += ")";
		}
	}

	return output;
}


// special matrix constructions as functions to be used in assignment
// don't duplicate code though
Matrix eye(int numrows, int numcols)
{
	Matrix output(numrows,numcols,"eye");
	return output;
}


Matrix eye(int size)
{
	Matrix output(size,"eye");
	return output;
}


Matrix diag(int numrows, int numcols, double initVal)
{
	Matrix output(numrows,numcols,"diag",initVal);
	return output;
}


Matrix diag(int size, double initVal)
{
	Matrix output(size,"diag",initVal);
	return output;
}


// square matrix of size v.mSize w/ entries of v on diagonal
Matrix diag(const Vector& v)
{
	int n = v.mSize;
	Matrix output(n);
	for (int i = 0; i < n; i++)
	{
		output.mData[i][i] = v.mData[i];
	}
	return output;
}

// matrix addition
Matrix operator+(const Matrix& A, const Matrix& B)
{
	// dimension checking
	if (A.mRows != B.mRows || A.mCols != B.mCols)
	{
		throw Exception("dimension mismatch", "Attempt to add Matrices of different sizes");
	}

	// create output Matrix and fill it in
	Matrix result(A.mRows,A.mCols);

	for (int row = 0; row < A.mRows; row++)
	{
		for (int col = 0; col < A.mCols; col++)
		{
			result.mData[row][col] = A.mData[row][col] + B.mData[row][col];
		}
	}

	return result;
}


// unary minus
Matrix operator-(const Matrix& A)
{
	// create the matix and fill it in
	Matrix result(A.mRows,A.mCols);

	for (int row = 0; row < A.mRows; row++)
	{
		for (int col = 0; col < A.mCols; col++)
		{
			result.mData[row][col] = -A.mData[row][col];
		}
	}

	return result;
}


// binary minus
// is unary minus of the second argument then addition
Matrix operator-(const Matrix& A, const Matrix& B)
{
	// don't duplicate code
	return A + (-B);
}


// matrix-matrix multiplication
Matrix operator*(const Matrix& A, const Matrix& B)
{
	// dimension checking
	if (A.mCols != B.mRows)
	{
		throw Exception("dimension mismatch",
		"Attempt to multiply incompatible quantities: inner Matrix dimensions must agree!");
	}

	Matrix output(A.mRows,B.mCols);

	for (int row = 0; row < A.mRows; row++)
	{
		for (int col = 0; col < B.mCols; col++)
		{	
			// do the dot product
			for (int j = 0; j < B.mRows; j++)
			{
				output.mData[row][col] += A.mData[row][j]*B.mData[j][col];
			}
		}
	}

	return output;
}


// matrix-scalar multiplication
Matrix operator*(const Matrix& A, const double& c)
{
	Matrix output(A);

	for (int row = 0; row < A.mRows; row++)
	{
		for (int col = 0; col < A.mCols; col++)
		{
			output.mData[row][col] *= c;
		}
	}
	return output;
}


// scalar-matrix multiplication
Matrix operator*(const double& c, const Matrix& A)
{
	// scalar multiplication commutes: don't duplicate code
	return A*c;
}


// mat-vec multiplication
Vector operator*(const Matrix& A, const Vector& x)
{
	// dimension checking
	if (A.mCols != x.mSize)
	{
		throw Exception("dimension mismatch", 
	"Attempt to multiply incompatible quantities: Matrix must have same no. cols as length of Vector");
	}

	// make the vector the right dimension :)
	Vector output(A.mRows);

	// for every row
	for (int i = 0; i < A.mRows; i++)
	{
		// do the dot product
		for (int j = 0; j < x.mSize; j++)
		{
			output.mData[i] += A.mData[i][j]*x.mData[j];
		}
	}

	return output;
}


// vector-matrix multiplication
Vector operator*(const Vector& x, const Matrix& A)
{
	// dimension checking
	if (x.mSize != A.mRows)
	{
		throw Exception("dimension mismatch",
		 "Attempt to multiply incompatible quantities: Vector must have length equal to no. rows of Matrix");
	}

	// make sure the output is the right size :)
	Vector output(A.mCols);

	// for every column
	for (int i = 0; i < A.mCols; i++)
	{
		// do the dot product
		for (int j = 0; j < A.mRows; j++)
		{
			output.mData[i] += x.mData[j]*A.mData[j][i];
		}
	}

	return output;
}


// only implement scalar division of matrices
Matrix operator/(const Matrix& A, const double& c)
{
	// parameter checks
	if (c == 0.0)
	{
		throw Exception("div 0", "Attempt to divide by 0");
	}

	// do this as a scalar multiplication: don't duplicate code
	double d = 1.0/c;

	return A*d;
}


// functions for finding the properties of the matrix
int* size(const Matrix& A)
{
	return A.size();
}


int length(const Matrix& A)
{
	return A.length();
}


int width(const Matrix& A)
{
	return A.width();
}


int height(const Matrix& A)
{
	return A.height();
}


bool is_symmetric(const Matrix& A)
{
	return A.is_symmetric();
}


bool is_singular(const Matrix& A)
{
	return A.is_singular();
}


int* Matrix::size() const
{
	int* a;
	a = new int[2];
	a[0] = mRows;
	a[1] = mCols;
	return a;
}


int Matrix::length() const
{
	return mRows*mCols;
}


int Matrix::width() const
{
	return mCols;
}


int Matrix::height() const
{
	return mRows;
}


bool Matrix::is_symmetric() const
{
	// assume symmetric
	bool symmetric = 1;
	
	// compare entries with the transpose
	Matrix AT(*this);
	AT.transpose();

	for (int row = 0; row < mRows; row++)
	{
		for (int col = 0; col < mCols; col++)
		{
			// if any entry not symmetric, break
			if (mData[row][col] != AT.mData[row][col])
			{
				symmetric = 0;
				break;
			} // end if
		} // end for
	} // end for

	return symmetric;
}


// uses LU to find determinant
bool Matrix::is_singular() const
{
	// assume not singular
	bool singular = 0;

	// check determinant of matrix
	double this_det;
	this_det = this->det();

	// if det is O(eps), assume it's zero
	if (fabs(this_det) <= 1E-15)
	{
		singular = 1;
	}

	return singular;
}


// norm member method
double Matrix::norm(int p) const
{
	switch(p)
	{
		case 1: // matrix 1-norm is max absolute column sum
		{
			int n = mCols, m = mRows;
			double running_sum, running_max = 0.0;
			
			for (int j = 0; j < n; j++)
			{
				// find each column sum
				running_sum = 0.0;
				for (int i = 0; i < m; i++)
				{
					running_sum += fabs(mData[i][j]);
				}

				// keep track of the max; no need to remember indices
				running_max = fmax(running_max, running_sum);
			}
			return running_max;
		}
		case 2: // matrix 2-norm is max spectral norm
		{
			throw Exception("not implemented","Matrix 2-norm not yet implemented");
		}
		default: // default behaviour if no case is met i.e. if we don't get a valid parameter
		{
			throw Exception("not implemented",
			"Invalid Matrix norm: only 1, infinity, Frobenius matrix norms currently supported!");
		}
	}
}


// norms with non-numeric names
double Matrix::norm(string type) const
{
	// parameter checking
	if (type != "infinity" && type != "Frobenius")
	{
		throw Exception("not implemented",
		"Invalid Matrix norm: only 1, infinity, Frobenius matrix norms currently supported!");
	}

	int m = mRows, n = mCols;

	// why have we changed design patterns from above??
	if (type == "infinity") // infinity norm is max absolute row sum
	{
		double running_sum, running_max = 0.0;
		for (int i = 0; i < m; i++)
		{
			running_sum = 0.0;
			for (int j = 0; j < n; j++)
			{
				running_sum += fabs(mData[i][j]);
			}
			// don't need to keep track of index
			running_max = fmax(running_max, running_sum);
		}

		return running_max;
	}
	else // Frobenius norm is basically vector 2 norm on all entries of Matrix
	{
		double running_sum = 0.0;
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < n; j++)
			{
				running_sum += mData[i][j]*mData[i][j];
			}
		}

		return sqrt(running_sum);
	}

}

// norm friend functions
double norm(const Matrix& A, int p)
{
	return A.norm(p);
}


double norm(const Matrix& A, string type)
{
	return A.norm(type);
}


// det method
double Matrix::det() const
{
	// dimension checks --> only square matrices have dets
	if (mRows != mCols)
	{
		throw Exception("dimension mismatch",
			"Attempt to find det of non-square Matrix A");
	}

	double det = 1.0;

	if (mRows == 2)
	{
		det = (mData[0][0] * mData[1][1]) - (mData[0][1] * mData[1][0]);
	}
	else
	{
		// we find det with LU, so initialise the necessary matrices
		Matrix A(*this);
		Matrix U(mRows);
		int perm_count;

		// gets the two important quantities from LUP
		perm_count = lu_det(A, U);

		// det(A) = det(L)det(U)(-1)^s
		// det(L) = 1
		// for s = number of row swaps

		// find the product of U diagonals
		for (int i = 0; i < mRows; i++)
		{
			det *= U.mData[i][i];
		}

		// apply the correct sign
		if (perm_count%2 != 0)
		{
			det *= -1;
		}
	}

	return det;
}


// det friend function
double det(const Matrix& A)
{
	return A.det();
}


// friend function transpose returns a new Matrix
Matrix transpose(const Matrix& A)
{
	Matrix AT(A.mCols,A.mRows);

	for (int i = 0; i < A.mRows; i++)
	{
		for (int j = 0; j < A.mCols; j++)
		{
			AT.mData[j][i] = A.mData[i][j];
		}
	}

	return AT;
}


// member method transpose does it in-place
Matrix& Matrix::transpose()
{
	// member method has to have square matrices!
	if (mRows != mCols)
	{
		throw Exception("illegal action",
			"Attempt to call in-place transpose on non-square Matrix");
	}

	double temp;

	for (int i = 0; i < mRows; i++)
	{
		for (int j = i; j < mCols; j++)
		{
			temp = mData[i][j];
			mData[i][j] = mData[j][i];
			mData[j][i] = temp;
		}
	}

	return *this;
}


// outer product of two vectors
Matrix outer(const Vector& v1, const Vector& v2)
{
    
    // dimension check
    if (v1.mSize != v2.mSize)
    {
      throw Exception("dimension mismatch",
        "Attempt to take outer product of different length Vectors.");
    }
    
    Matrix output(v1.mSize);

    // do the multiplication
    for (int i = 0; i < v1.mSize; i++)
    {
        for (int j = 0; j < v2.mSize; j++)
        {
            output.mData[i][j] = v1.mData[i]*v2.mData[j];
        }
    }

    return output;
}


// linear solver functions
Vector backsub(const Matrix& A, const Vector& b)
{
	// A must be upper triangular

	// A must be square and same size as b
	if (A.mRows != A.mCols)
	{
		throw Exception("dimension mismatch",
			"Attempt to back sub with non-square matrix");
	}
	if (A.mRows != b.mSize)
	{
		throw Exception("dimension mismatch",
			"Attempt to back sub incompatible Matrix and Vector");
	}

	// initialise result
	Vector x(b.mSize);

	// do the back sub
	for (int i = b.mSize-1; i >= 0; i--)
	{	
		x.mData[i] = b.mData[i];
		for (int j = i+1; j < b.mSize; j++)
		{
			x.mData[i] -= x.mData[j]*A.mData[i][j];
		}
		x.mData[i] /= A.mData[i][i];
	}

	return x;
}


Vector forwardsub(const Matrix& A, const Vector& b)
{
	// A must be lower triangular

	// A must be square and same size as b
	if (A.mRows != A.mCols)
	{
		throw Exception("dimension mismatch",
			"Attempt to forward sub with non-square matrix");
	}
	if (A.mRows != b.mSize)
	{
		throw Exception("dimension mismatch",
			"Attempt to forward sub incompatible Matrix and Vector");
	}

	// initialise result
	Vector x(b.mSize);

	// do the forward sub
	for (int i = 0; i < b.mSize; i++)
	{
		x.mData[i] = b.mData[i];
		for (int j = 0; j < i; j++)
		{
			x.mData[i] -= x.mData[j]*A.mData[i][j];
		}
		x.mData[i] /= A.mData[i][i];
	}

	return x;
}


// QR factorisation using Householder reflectors, overwrites A --> R
// returns Matrix W containing reflectors by default
// if explicit_Q = 1, then explicitly form Q
// WARNING: explicitly forming Q takes a lot longer!
Matrix qr(Matrix& A, bool explicit_Q)
{
	// currently only for square matrices
	if (A.mRows != A.mCols)
	{
		throw Exception("not implemented",
			"QR only implemented for square matrices at this time.");
	}
	
	// initialise some storage
	int counter, counter2; // for indexing
	double normV;

	// Initialise matrices W, Q
	Matrix W(A.mRows);
	Matrix Q(A.mRows, "eye");

	// allocate an extra storage column
	Vector b(A.mRows);

	// big loop
	for (int k = 0; k < A.mCols; k++)
	{
		// memory for the kth reflector
		Vector V(A.mRows - k);

		// assign some entries of V
		counter = 0;
		for (int i = k; i < A.mCols; i++)
		{
			V.mData[counter] = A.mData[i][k];
			counter++;
		}
		
		// norm of V
		normV = V.norm(2);
		
		// modify the first entry if needed
		if (V.mData[0] >= 0.0)
		{
			V.mData[0] += normV;
		}
		else
		{
			V.mData[0] -= normV;
		}

		normV = V.norm(2); // new norm of V
		V = V/normV; // normalise V

		// do the matrix multiplication on a slice of A
		Matrix Aslice(V.mSize);
		Aslice = slice(A, k + 1, A.mRows, k + 1, A.mCols);
		Aslice = 2 * outer(V, (V * Aslice));

		// do the reflection
		counter = 0;
		for (int i = k; i < A.mRows; i++)
		{
			counter2 = 0;
			for (int j = k; j < A.mCols; j++)
			{
				A.mData[i][j] -=  Aslice.mData[counter][counter2];
				counter2++;
			}
			W.mData[i][k] = V.mData[counter]; // V has the kth reflector
			counter++;
		}
	}

	// zero out the entries that are definitely 0
	for (int col = 0; col < A.mCols; col++)
	{
		for (int row = col + 1; row < A.mRows; row++)
		{
			A.mData[row][col] = 0.0;
		}
	}

	if (explicit_Q)
	{
		for (int j = 0; j < A.mRows; j++)
		{
			// pick the relevant column of Q
			b.mData[j] = 1.0;

			// subtract off previous reflectors
			for (int l = A.mRows - 1; l >= 0; l--)
			{
				// initialise the right size vector to pick reflectors
				Vector V(A.mRows - l);
				
				// pick the lth reflector from W
				counter = 0;
				for (int ii = l; ii < A.mRows; ii++)
				{
					V.mData[counter] = W.mData[ii][l];
					counter++;
				}
				
				// subtract off this reflector
				// get the relevant matrix product
				Matrix VdotVT(V.mSize);
				VdotVT = outer(V, V);

				// and operate on the correct slice of b
				Vector bslice(A.mRows - l);
				bslice = slice(b, l + 1, A.mRows);
				bslice = 2.0 * V * (V * bslice);

				// subtract off the correct value
				counter = 0;
				for (int ii = l; ii < A.mRows; ii++)
				{
					b.mData[ii] -= bslice.mData[counter];
					counter++;
				}
			}

			// assign the new column of Q and reset b
			for (int i = 0; i < A.mRows; i++)
			{
				Q.mData[i][j] = b.mData[i];
				b.mData[i] = 0.0;
			}
		}

		return Q;
	}
	else
	{
		return W;
	}
}


// LU decomp w/o pivotting, A unchanged, U filled in
Matrix lu(const Matrix& A, Matrix& U)
{
	// dimension checks
	if (U.mRows != A.mRows || U.mCols != A.mCols)
	{
		throw Exception("dimension mismatch",
			"Please provide U with same dimensions as A");
	}
	int m = A.mRows, n = A.mCols;

	// initialise L
	Matrix L(m);

	// big loop
	for (int i = 0; i < m; i++)
	{
		// do the row subtraction and record multipliers
		// use the doolittle algorithm
		
		// fill in row of U
		for (int k = i; k < n; k++)
		{ // use the doolittle formulas
			U.mData[i][k] = A.mData[i][k];
			for (int j = 0; j < m; j++)
			{
				U.mData[i][k] -= L.mData[i][j] * U.mData[j][k];

			} // end for j
		} // end for k

		// fill in col of L
		for (int k = i; k < m; k++)
		{
			if (k == i) // diagonals are 1
			{
				L.mData[i][i] = 1.0;
			}
			else // otherwise use the doolittle formulas
			{
				L.mData[k][i] = A.mData[k][i];
				for (int j = 0; j < i; j++)
				{
					L.mData[k][i] -= L.mData[k][j] * U.mData[j][i];

				} // end for j
				L.mData[k][i] /= U.mData[i][i];
			} // end if
		} // end for k

	} // end for i

	return L;
}


// LU  decomp with partial pivoting, overwrites A with U, returns num permutes
int lu_det(Matrix& A, Matrix& U)
{

	// dimension checks fine since we implement this ourselves
	int m = A.mRows, n = A.mCols; // shorthand

	int p , perm_count = 0; // track pivot locations, counts
	double running_val, tmp; // temporary variables

	// initialise L
	Matrix L(m);

	Matrix P(m);

	// big loop
	for (int i = 0; i < m; i++)
	{
		// step 1: find pivot as greatest entry in the column
		running_val = 0.0;
		for (int j = i; j < m; j++)
		{
			// record the index of row with largest entry in i'th column
			if (fabs(A.mData[j][i]) > running_val)
			{
				running_val = fabs(A.mData[j][i]);
				p = j;
			}
		} // end for j

		// step 2: perform the row swap in A, L, and P
		if (p != i)
		{
			perm_count++;
			for (int j = 0; j < n; j++)
			{
				// swap A
				tmp = A.mData[i][j];
				A.mData[i][j] = A.mData[p][j];
				A.mData[p][j] = tmp;

				// swap P
				tmp = P.mData[i][j];
				P.mData[i][j] = P.mData[p][j];
				P.mData[p][j] = tmp;

				// swap L	
				tmp = L.mData[i][j];
				L.mData[i][j] = L.mData[p][j];
				L.mData[p][j] = tmp;
			} // end for j
		}

		// step 3: do the row subtraction and record multipliers		
		// fill in row of U
		for (int k = i; k < n; k++)
		{ // use the doolittle formulas
			U.mData[i][k] = A.mData[i][k];
			for (int j = 0; j < i; j++)
			{
				U.mData[i][k] -= L.mData[i][j] * U.mData[j][k];

			}
		}

		// fill in col of L
		for (int k = i; k < m; k++)
		{
			if (k == i) // diagonals are 1
			{
				L.mData[i][i] = 1.0;
			}
			else // use the doolittle formulas
			{
				L.mData[k][i] = A.mData[k][i];
				for (int j = 0; j < i; j++)
				{
					L.mData[k][i] -= L.mData[k][j] * U.mData[j][i];

				}
				L.mData[k][i] /= U.mData[i][i];
			}
		}

	}

	// un-permute A so it remains as the original
	A = P*A;

	return perm_count;

}


// LU  decomp with partial pivoting, overwrites A with U, outputs L
Matrix lu(Matrix& A, Matrix& U, Matrix& P)
{

	// dimension checks
	if (P.mRows != P.mCols)
	{
		throw Exception("dimension mismatch",
			"Please provide square identity matrix for P");
	}
	if (P.mCols != A.mRows)
	{
		throw Exception("dimension mismatch",
			"Please provide (mxm) identity P for (mxn) A");
	}
	if (U.mRows != A.mRows || U.mCols != A.mCols)
	{
		throw Exception("dimension mismatch",
			"Please provide U with same dimensions as A");
	}
	if (A.mRows > A.mCols)
	{
		throw Exception("not implemented",
			"LU for (mxn) A, m > n, not currently implemented.");
	}

	int m = A.mRows, n = A.mCols; // shorthand

	int p; // track pivot locations
	double running_val, tmp; // temporary variables

	// initialise L
	Matrix L(m);

	// big loop
	for (int i = 0; i < m; i++)
	{
		// step 1: find pivot as greatest entry in the column
		running_val = 0.0;
		for (int j = i; j < m; j++)
		{
			// record the index of row with largest entry in i'th column
			if (fabs(A.mData[j][i]) > running_val)
			{
				running_val = fabs(A.mData[j][i]);
				p = j;
			}
		} // end for j

		// step 2: perform the row swap in A, L, and P
		for (int j = 0; j < n; j++)
		{
			// swap A
			tmp = A.mData[i][j];
			A.mData[i][j] = A.mData[p][j];
			A.mData[p][j] = tmp;

			// swap P
			tmp = P.mData[i][j];
			P.mData[i][j] = P.mData[p][j];
			P.mData[p][j] = tmp;

			// swap L	
			tmp = L.mData[i][j];
			L.mData[i][j] = L.mData[p][j];
			L.mData[p][j] = tmp;
		} // end for j

		// step 3: do the row subtraction and record multipliers		
		// fill in row of U
		for (int k = i; k < n; k++)
		{ // use the doolittle formulas
			U.mData[i][k] = A.mData[i][k];
			for (int j = 0; j < i; j++)
			{
				U.mData[i][k] -= L.mData[i][j] * U.mData[j][k];

			}
		}

		// fill in col of L
		for (int k = i; k < m; k++)
		{
			if (k == i) // diagonals are 1
			{
				L.mData[i][i] = 1.0;
			}
			else // use the doolittle formulas
			{
				L.mData[k][i] = A.mData[k][i];
				for (int j = 0; j < i; j++)
				{
					L.mData[k][i] -= L.mData[k][j] * U.mData[j][i];

				}
				L.mData[k][i] /= U.mData[i][i];
			}
		}

	}

	// un-permute A so it remains as the original
	A = P*A;

	return L;

}

// iterative methods
Vector cg(const Matrix& A, const Vector& b, bool verbose, double TOL, int maxit)
{
	// if no initial guess is given, just use the empty vector
	Vector x0(b.mSize);
	return cg(A, b, x0, verbose, TOL, maxit);
}

Vector cg(const Matrix& A, const Vector& b, const Vector& x0, bool verbose, double TOL, int maxit)
{

	// dimension checks
	if (A.mRows != A.mCols)
	{
		throw Exception("dimension mismatch",
			"Attempt to use Arnoldi iteration on non-square matrix.");
	}
	if (A.mCols != b.mSize)
	{
		throw Exception("dimension mismatch",
			"Attempt to use Arnoldi on incompatible A, b.");
	}
	// parameter checks
	if (maxit < 1)
	{
		throw Exception("invalid parameter",
			"Must choose max no. iterations >= 1.");
	}
	if (TOL <= 0)
	{
		throw Exception("invalid parameter",
			"Error tolerance must be > 0.");
	}
	// disable for speed comparisons
	// if (A.is_singular())
	// {
	// 	throw Exception("invalid parameter",
	// 		"Attempt to use cg with singular A");
	// }
	// if (!A.is_symmetric())
	// {
	// 	throw Exception("invalid parameter",
	// 		"Attempt to use cg with non-symmetric A");
	// }

	// useful storage
	int k = 0;
	double alpha, beta, residual, rk2;
	Vector xk(x0), xkp1(b.mSize); // xk will start with x0
	Vector rk(b.mSize), rkp1(b.mSize), pkp1(b.mSize); // useful storage
	Vector Apk(b.mSize); // to store the mat-vec multiplication

	// starting points
	rk = b - A*xk;
	Vector pk(rk); // p_0 == r_0

	// big loop
	while (k < maxit)
	{
		// find useful products
		Apk = A*pk;
		rk2 = rk*rk;

		// do the CG calculation
		alpha = (rk2) / (pk * Apk);
		xkp1 = xk + (alpha * pk);
		rkp1 = rk - (alpha * Apk);
		residual = rkp1.norm();

		if (residual < TOL) // check if we converge
		{
			break;
		}

		beta = (rkp1 * rkp1)/(rk2);
		pkp1 = rkp1 + (beta * pk);
		
		// update the iterates
		k++;
		xk = xkp1;
		rk = rkp1;
		pk = pkp1;
	}

	// informative output
	if (verbose)
	{
		if (k > maxit)
		{
			cout << "Conjugate Gradients terminated after reaching maximum iteration count,\n";
			cout << "Residual of solution was = " << residual << "\n";
		}
		else
		{
			cout << "Conjugate Gradients converged to tolerance after " << k << " iteration(s) ";
			cout << "with relative residual " << residual << "\n";
		}
	}
	
	return xkp1;
}

// does the kth iteration of arnoldi and returns the vector hk, updates Q in place
Vector arnoldi(const Matrix& A, Matrix& Q, int k)
{
	// dimension checks
	if (A.mRows != A.mCols)
	{
		throw Exception("dimension mismatch",
			"Attempt to use Arnoldi iteration on non-square matrix.");
	}


	// the output is the next column of the matrix H
	Vector hk(k+1);

	// Vector qk is column k of matrix Q
	Vector qk(Q.mRows);
	qk = col_slice(Q, k);

	qk = A*qk; // form the next vector in the Krylov subspace

	// orthogonalise
	for (int j = 0; j < k; j++)
	{
		// get jth col of Q = qj
		Vector qj(Q.mRows);
		qj = col_slice(Q, j + 1);

		// remove the component in qj direction from qk
		hk.mData[j] = qj * qk;
		qk = qk - (hk.mData[j] * qj);
	}

	// normalise
	hk.mData[k] = qk.norm();
	if (hk.mData[k] < 1E-12 || k == A.mRows) // the convergence critera; we have reproduced the Krylov space
	{
		hk.mData[k] = 0;
	}
	else
	{
		qk = qk / hk.mData[k];
	}

	// fill in Q with the new value
	for (int i = 0; i < Q.mRows; i++)
	{
		Q.mData[i][k] = qk.mData[i];
	}

	return hk;
}


Vector gmres(const Matrix& A, const Vector& b, bool verbose, double TOL, int maxit)
{
	Vector x0(b.mSize); // if no initial guess, start at 0
	return gmres(A, b, x0, verbose, TOL, maxit);
}


Vector gmres(const Matrix& A, const Vector& b, const Vector& x0, bool verbose, double TOL, int maxit)
{
	// dimension checks
	if (A.mRows != A.mCols)
	{
		throw Exception("dimension mismatch",
			"Attempt to use Arnoldi iteration on non-square matrix.");
	}
	if (A.mCols != b.mSize)
	{
		throw Exception("dimension mismatch",
			"Attempt to use Arnoldi on incompatible A, b.");
	}
	// parameter checks
	if (maxit < 1)
	{
		throw Exception("invalid parameter",
			"Must choose max no. iterations >= 1.");
	}
	if (TOL <= 0)
	{
		throw Exception("invalid parameter",
			"Error tolerance must be > 0.");
	}
	// disabled to test convergence speed
	// if (A.is_singular())
	// {
	// 	throw Exception("invalid parameter",
	// 		"Attempt to use gmres with singular A");
	// }
	
	bool stop = 0; // stop criteria flag
	int k = 0; // iteration number
	double tmp1, tmp2, rkp2_p, normb = b.norm(); // declare useful variables
	Matrix H(maxit + 1, maxit + 1); // allocate all the memory needed for H
	Vector hk(b.mSize), normbe1(maxit + 1); // storage for result of arnoldi and "RHS" of LLS
	normbe1.mData[0] = normb;
	
	// storage for the Givens refections
	double* sink;
	sink = new double[maxit];
	sink[0] = 1.0;
	double* cosk;
	cosk = new double[maxit];
	cosk[0] = 1.0;

	// allocate memory to hold the basis vectors of the Krylov space
	Matrix Q(A.mRows, maxit + 1);
	
	double residual = (b - A*x0).norm(); // initialise residual

	// fill in the first columns
	H.mData[0][0] = normb;
	for (int i = 0; i < A.mRows; i++)
	{
		Q.mData[i][0] = b.mData[i] / normb;
	}

	while (k < maxit)
	{
		k++;
		// get qk via Arnoldi

		// get the next column of H from Arnoldi
		Vector hk(k + 1);
		hk = arnoldi(A, Q, k);

		if (hk.mData[k] < 1E-15) // redundant check but don't compare to zero!
		{
			// this is a stop criteria
			stop = 1;
		}


		// apply the old Givens rotation to the new column of H
		if (k > 1)
		{
			for (int j = 0; j < k-1; j++)
			{
				tmp1 = (cosk[j] * hk.mData[j]) - (sink[j] * hk.mData[j+1]);
				tmp2 = (sink[j] * hk.mData[j]) + (cosk[j] * hk.mData[j+1]);
				hk.mData[j] = tmp1;
				hk.mData[j+1] = tmp2;

			}
		}

		
		// find Givens rotation
		rkp2_p = sqrt((hk.mData[k-1] * hk.mData[k-1]) + (hk.mData[k] * hk.mData[k]));
		cosk[k-1] = hk.mData[k-1] / rkp2_p;
		sink[k-1] = -hk.mData[k] / rkp2_p;

		// now do Givens rotation on the new column
		hk.mData[k-1] = (cosk[k-1] * hk.mData[k-1]) - (sink[k-1] * hk.mData[k]);
		hk.mData[k] = 0;

		// apply the Givens rotation to norm(b)*e1:
		tmp1 = (cosk[k-1] * normbe1.mData[k-1]) - (sink[k-1] * normbe1.mData[k]);
		tmp2 = (sink[k-1] * normbe1.mData[k-1]) + (cosk[k-1] * normbe1.mData[k]);
		normbe1.mData[k-1] = tmp1;
		normbe1.mData[k] = tmp2;
		residual = fabs(tmp2);

		// store the new column in H
		for (int i = 0; i < k; i++)
		{
			H.mData[i][k-1] = hk.mData[i];
		}

		if (stop)
		{
			// break if we reach the 
			break;
		}

		// check the residual
		if (residual < TOL)
		{
			// only explicitly calculate the solution if the residual is small enough
		  	break;
		}
	}

	// informative output
	if (verbose)
	{
		if (k >= maxit)
		{
			cout << "GMRES terminated after reaching maximum iteration count.\n";
			cout << "Residual of solution was = " << residual << "\n";
		}
		else if (stop)
		{
			cout << "GMRES terminanted after reaching the dimension of the matrix.\n";
			cout << "Residual of solution was = " << residual << "\n";
		}
		else
		{
			cout << "GMRES converged to tolerance after " << k << " iteration(s) ";
			cout << "with relative residual " << residual << "\n";
		}
	}

	// trim extra values to get right dimensions
	Matrix Qk(Q.mRows, k);
	Qk = slice(Q, 1, Q.mRows, 1, k);

	Matrix Hk(k, k);
	Hk = slice(H, 1, k, 1, k);

	Vector normbe1k(k);
	normbe1k = slice(normbe1, 1, k);
	
	// and solve the problem
	Vector yk(k);
	yk = backsub(Hk, normbe1k);

	return x0 + (Qk * yk);
}
