// Author: Nicholas Barton

#include <iostream>
#include <string>
#include <random>
#include <chrono>
#include "Matrix.hpp"
#include "Vector.hpp"
#include "Exception.hpp"

using namespace std;
using namespace std::chrono;

bool True = 1;

// function prototypes for the timing tests to solve Ax = b
Matrix build_A_matrix(int size); // builds a big SPD matrix A
Vector build_RHS(int size); // builds a RHS b
double lu_time_function(Matrix& A, const Vector& b); // does a timing test for LU
double qr_time_function(Matrix& A, const Vector& b); // does a timing test for QR
double cg_time_function(const Matrix& A, const Vector& b); // does a timing test for CG
double gmres_time_function(const Matrix& A, const Vector& b); // does a timing test for GMRES

// driver function
int main()
{
	try
	{
	int size = 256;
	assert(size%8 == 0); // make sure we have a size that can be used

	Matrix Ahuge(size), Abig(size / 2), Asmall(size / 4), Atiny(size / 8);

	Ahuge = build_A_matrix(size);
	Abig = slice(Ahuge, 1, size / 2, 1, size / 2);
	Asmall = slice(Ahuge, 1, size / 4, 1, size / 4);
	Atiny = slice(Ahuge, 1, size / 8, 1, size / 8);

	// we want to test how long it takes to solve the system using our different solvers
	// we have 4 solvers and 4 different matrices -- 16 total tests
	Matrix times(4);
	
	// build the RHS
	Vector bhuge(size), bbig(size / 2), bsmall(size / 4), btiny(size / 8);
	bhuge = build_RHS(size);
	bbig = slice(bhuge, 1, size / 2);
	bsmall = slice(bhuge, 1, size / 4);
	btiny = slice(bhuge, 1, size / 8);
	
	// run the tests and store them in times Matrix
	times(1, 1) = lu_time_function(Atiny, btiny);
	times(1, 2) = lu_time_function(Asmall, bsmall);
	times(1, 3) = lu_time_function(Abig, bbig);
	times(1, 4) = lu_time_function(Ahuge, bhuge);

	times(2, 1) = qr_time_function(Atiny, btiny);
	times(2, 2) = qr_time_function(Asmall, bsmall);
	times(2, 3) = qr_time_function(Abig, bbig);
	times(2, 4) = qr_time_function(Ahuge, bhuge);

	times(3, 1) = cg_time_function(Atiny, btiny);
	times(3, 2) = cg_time_function(Asmall, bsmall);
	times(3, 3) = cg_time_function(Abig, bbig);
	times(3, 4) = cg_time_function(Ahuge, bhuge);

	times(4, 1) = gmres_time_function(Atiny, btiny);
	times(4, 2) = gmres_time_function(Asmall, bsmall);
	times(4, 3) = gmres_time_function(Abig, bbig);
	times(4, 4) = gmres_time_function(Ahuge, bhuge);

	// just look: no further explanation given here; discussed elsewhere
	cout << times << "\n";

	}
	catch(Exception &ex)
	{
		ex.DebugPrint();
	}

	return 0;
}


// builds a random SPD Matrix of given size
Matrix build_A_matrix(int size)
{
	// a uniform random number generator 
	random_device rd; // gets a seed for the RNG
	mt19937 gen(rd()); // RNG engine
	uniform_real_distribution<> distr(0.0, 1.0); // U([0, 1))

	Matrix A(size);
	Matrix pos(size,"diag",(double) size); // to make the matrix > 0

	for (int row = 1; row <= size; row++)
	{
		for (int col = 1; col <= size; col++)
		{
			// random entries in [0, 1)
			A(row,col) = distr(gen);
		}
	}

	// return a Matrix that is SPD
	return (A * transpose(A)) + pos;
}


Vector build_RHS(int size)
{
	// a uniform random number generator 
	random_device rd; // gets a seed for the RNG
	mt19937 gen(rd()); // RNG engine
	uniform_real_distribution<> distr(0.0, 1.0); // U([0, 1))

	Vector b(size);

	for (int i = 1; i <= size; i++)
	{
		b(i) = 2.0 * size * distr(gen);
	}

	return b;
}

// time the methods! not concerned with actual numerical result, just how long it takes
// since we test the numerical accuracy in test_matrix
// we run each method 30 times and find the average time
double lu_time_function(Matrix& A, const Vector& b)
{
	high_resolution_clock::time_point t1, t2;

	// track progress
	cout << "Running an LU time test...\n";
	flush(cout);

	Matrix L(height(A));
	Matrix P(height(A), "eye");
	Matrix U(height(A));
	Matrix zero(height(A)); // to reset U
	Matrix I(height(A), "eye"); // to reset P
	Vector x(height(A)); // for the solution
	Vector y(height(A)); // intermediary solution

	// start the timer
	t1 = high_resolution_clock::now();
	
	// solve the system with LUP
	for (int test = 0; test < 30; test++)
	{
		L = lu(A, U, P);
		y = backsub(U, b);
		x = forwardsub(L, y);
		x = P.transpose() * x;

		// reset matrices
		P = I;
		U = zero;
	}
	
	// stop the timer and determine average time
	t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( (t2 - t1) / 30 ).count();
	return duration;
}


double qr_time_function(Matrix& A, const Vector& b)
{
	high_resolution_clock::time_point t1, t2;

	Matrix W(height(A));
	Matrix A_copy(A);
	Vector y(length(b)), x(length(b));

	// track progress
	cout << "Running a QR time test...\n";
	flush(cout);
	
	t1 = high_resolution_clock::now();

	// solve the system with QR
	// since we don't care about values, just speed we don't explicitly form Q
	// could apply the Householder reflectors in W to y "properly" but it's still O(mat-vec mult)
	for (int test = 0; test < 30; test++)
	{
		W = qr(A_copy);
		y = backsub(A_copy, b);
		x = W.transpose() * y;

		// reset matrices
		A_copy = A;
	}

	// stop the timer and determine average time
	t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( (t2 - t1) /30 ).count();

	return duration;
}


double cg_time_function(const Matrix& A, const Vector& b)
{
	high_resolution_clock::time_point t1, t2;

	// track progress
	cout << "Running a CG time test...\n";
	flush(cout);

	Vector x(length(b));

	t1 = high_resolution_clock::now();

	for (int test = 0; test < 30; test++)
	{
		if (test == 29) // see how many iterations it needs
		{
			x = cg(A, b, True);
		}
		else
		{
			x = cg(A, b);
		}
		
	}

	// stop the timer and determine average time
	t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( (t2 - t1) /30 ).count();

	return duration;
}


double gmres_time_function(const Matrix& A, const Vector& b)
{
	high_resolution_clock::time_point t1, t2;

	// track progress
	cout << "Running a GMRES time test...\n";
	flush(cout);

	Vector x(height(A));

	// start the timer
	t1 = high_resolution_clock::now();

	// run the test
	for (int test = 0; test < 30; test++)
	{
		if (test == 29) // see how many iterations it needs
		{
			x = gmres(A, b, True);
		}
		else
		{
			x = gmres(A, b);
		}
	}

	// stop the timer and determine average time
	t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( (t2 - t1) /30 ).count();

	return duration;
}
