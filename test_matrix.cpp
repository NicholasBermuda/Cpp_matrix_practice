// Author: Nicholas Barton

#include <iostream>
#include <string>
#include <cassert>
#include "Matrix.hpp"
#include "Vector.hpp"
#include "Exception.hpp"

using namespace std; 

// function prototypes return number of tests failed
void test_constructors_and_index();
void test_slice();
void test_size_properties();
void test_to_string();
void test_binary_operators();
void test_norm();
void test_transpose();
void test_factorisations();
void test_det();
void test_triangular_solvers();
void test_iterative_solvers();
Vector build_test_vec();
Matrix build_test_mat();
Vector build_test_soln();


// tests the main functionality of the different features of the code
int main()
{
	// do all the tests -- we catch all exceptions and errors in the functions
	// we order the tests so we can use the results from previous tests
	test_constructors_and_index();
	test_slice();
	test_size_properties();
	test_to_string();
	test_binary_operators();
	test_norm();
	test_transpose();
	test_factorisations();
	test_det();
	test_triangular_solvers();
	test_iterative_solvers(); 
	
	cout << "All tests finished!\n";

	return 0;
}

void test_constructors_and_index()
{
	cout << "Testing constructors and indexing... ";
	flush(cout);

	try
	{
		// zero and init value constructors
		Matrix A(3);
		Matrix B(3,3);
		Matrix C(3,1.0);
		Matrix D(3,3,1.0);

		// indexing and testing values of above
		assert(A(3) == 0.0);
		assert(A(2,2) == 0.0);
		assert(C(1) == 1.0);
		assert(D(3,1) == 1.0);

		// assignment and special constuctor functions (and assignment operator shape checks!)
		// these "friend" constructors call the "member method" constructors under the hood so this is a 2-4-1
		A = eye(3,3);
		B = eye(3);
		C = diag(3,3,4.3);
		D = diag(3,5.2);

		// test that they have the right form
		assert(A(1,1) == 1.0 && A(2,1) == 0.0 && A(3,2) == 0.0 && A(3,3) == 1.0);
		assert(B(1,1) == 1.0 && B(1,2) == 0.0 && B(2,3) == 0.0 && B(2,2) == 1.0);
		assert(C(1,1) == 4.3 && C(2,1) == 0.0 && C(3,2) == 0.0 && C(3,3) == 4.3);
		assert(D(1,1) == 5.2 && D(1,2) == 0.0 && D(2,3) == 0.0 && D(2,2) == 5.2);

		// assigning a vector to diag
		Vector v(3);
		v(1) = 1.;
		v(2) = 2.;
		v(3) = 3.;
		D = diag(v);
		assert(D(1,1) == 1.0 && D(1,2) == 0.0 && D(2,3) == 0.0 && D(2,2) == 2. && D(3,3) == 3.);

		// test copy constructor
		Matrix C_2(C);
		for(int i = 1; i <= 3; i++)
		{
			for(int j = 1; j <= 3; j++)
			{
				assert(C(i,j) == C_2(i,j));
			}

		}

	}
	catch(Exception &ex)
	{
		ex.DebugPrint();
		cout << "Constructor and indexing test FAIL\n";
		flush(cout);
		return;
	}
	cout << "Success!\n";
	flush(cout);
	
	return;
}


void test_slice()
{
	cout << "Testing slicing... ";
	flush(cout);

	// test the various slicing functions
	try
	{
		// sub-matrix slice
		Matrix A2(5, 4, 1.2);
		A2(1, 1) = 2;
		A2(1, 2) = 3;
		A2(2, 1) = 3;
		A2(2, 2) = 4;

		Matrix B(2, 2);
		B = slice(A2, 1, 2, 1, 2); // tests assignment size checks
		
		for(int i = 1; i <= 2; i++)
		{
			for(int j = 1; j <= 2; j++)
			{ // makes sure we get the right values
				assert(B(i, j) == i+j);
			}
		}

		// search for this value in slices
		A2(3,3) = 5.0;

		// column slice
		Vector Acol(5);
		Acol = col_slice(A2, 3);

		// find value
		for(int i = 1; i <= 5; i++)
		{
			if (i == 3)
			{
				assert(Acol(i) == 5.0);
			}
			else
			{
				assert(Acol(i) == 1.2);
			}
		}

		// row slice
		Vector Arow(4);
		Arow = row_slice(A2, 3);

		// find value
		for(int i = 1; i <= 4; i++)
		{
			if (i == 3)
			{
				assert(Arow(i) == 5.0);
			}
			else
			{
				assert(Arow(i) == 1.2);
			}
		}

	}
	catch (Exception &ex)
	{
		cout << "FAIL!\n";
		flush(cout);
		ex.DebugPrint();
		return;
	}
	cout << "Success!\n";
	flush(cout);

	return;
}


void test_size_properties()
{
	cout << "Testing size functions... ";
	flush(cout);
	
	// test the functions holding the size of a Matrix
	// the friend functions call member methods so it's a 2-4-1
	try
	{
		Matrix A(4, 5);

		assert(length(A) == 20);
		assert(width(A) == 5);
		assert(height(A) == 4);

		int* s;
		s = new int[2];
		s = size(A);
		assert(s[0] == 4);
		assert(s[1] == 5);
	}
	catch (Exception &ex)
	{
		cout << "FAIL!\n";
		flush(cout);
		ex.DebugPrint();
		return;
	}
	cout << "Success!\n";
	flush(cout);

	return;
}


void test_to_string()
{
	cout << "Testing conversion to string... ";
	flush(cout);

	try
	{
		Matrix A(3,"eye");
		string a_string = to_string(A);
		string results = "((1.000000, 0.000000, 0.000000), (0.000000, 1.000000, 0.000000), (0.000000, 0.000000, 1.000000))";
		assert(a_string == results);
	}
	catch (Exception &ex)
	{
		cout << "FAIL!\n";
		flush(cout);
		ex.DebugPrint();
		return;
	}
	cout << "Success!\n";
	flush(cout);

	return;
}


void test_binary_operators()
{
	cout << "Testing Matrix binary operators... ";
	flush(cout);
	
	try
	{
		Matrix A(3,"eye");
		Matrix B(3, 2.5);
		Matrix C(3); C = -A;
		Matrix D(3); D = A + B;
		Matrix E(3); E = B - A;
		Matrix F(3); F = 2.0 * B;
		Matrix G(3); G = B * 2.0;
		Matrix H(3); H = B / 2.5;


		for(int i = 1; i <= 3; i++)
		{
			for(int j = 1; j <= 3; j++)
			{
				if(i == j)
				{
					assert(C(i, j) == -1.0);
					assert(D(i, j) == 3.5);
					assert(E(i, j) == 1.5);
					assert(F(i, j) == 5.0);
					assert(G(i, j) == 5.0);
					assert(H(i, j) == 1.0);
				}
				else
				{
					assert(C(i,j) == 0.0);
					assert(D(i,j) == 2.5);
					assert(E(i,j) == 2.5);
					assert(F(i,j) == 5.0);
					assert(G(i,j) == 5.0);
					assert(H(i,j) == 1.0);
				}
			}
		}

		Matrix J(3,4);
		J(1, 1) = 1;
		J(1, 2) = 2;
		J(1, 3) = 3;
		J(1, 4) = 4;
		J(2, 1) = 4;
		J(2, 2) = 5;
		J(2, 3) = 6;
		J(2, 4) = 7;
		J(3, 1) = 7;
		J(3, 2) = 8;
		J(3, 3) = 9;
		J(3, 4) = 10;

		Matrix EJ(3,4);
		EJ = E*J;
		string ej_string = to_string(EJ);

		// result calculated on MATLAB
		string test_ej = "((29.000000, 35.500000, 42.000000, 48.500000)";
		test_ej += ", (26.000000, 32.500000, 39.000000, 45.500000)";
		test_ej += ", (23.000000, 29.500000, 36.000000, 42.500000))";

		assert(ej_string == test_ej);

		Vector x(3);
		x(1) = 1;
		x(2) = 2;
		x(3) = 3;
		x = E*x;
		assert(x(1) == 14.0);
		assert(x(2) == 13.0);
		assert(x(3) == 12.0);

		x = x*E;
		assert(x(1) == 83.5);
		assert(x(2) == 84.5);
		assert(x(3) == 85.5);
	}
	catch (Exception &ex)
	{
		cout << "FAIL!\n";
		flush(cout);
		ex.DebugPrint();
		return;
	}
	cout << "Success!\n";
	flush(cout);

	return;
}


void test_norm()
{
	cout << "Testing Matrix norm functions... ";
	flush(cout);

	try
	{
		Matrix J(3,4);
		J(1, 1) = 1;
		J(1, 2) = 2;
		J(1, 3) = 3;
		J(1, 4) = 4;
		J(2, 1) = 4;
		J(2, 2) = 5;
		J(2, 3) = 6;
		J(2, 4) = 7;
		J(3, 1) = 7;
		J(3, 2) = 8;
		J(3, 3) = 9;
		J(3, 4) = 10;
		assert(norm(J) == 21);
		assert(norm(J,1) == 21);
		assert(norm(J,"infinity") == 34);
		assert(norm(J,"Frobenius") - 21.2132034 < 1E-6);
	}
	catch (Exception &ex)
	{
		cout << "FAIL!\n";
		flush(cout);
		ex.DebugPrint();
		return;
	}
	cout << "Success!\n";
	flush(cout);

	return;
}


void test_transpose()
{
	cout << "Testing Matrix transpose... ";
	flush(cout);

	try
	{
		Matrix J(2, 3);
		J(1, 1) = 1;
		J(1, 2) = 2;
		J(1, 3) = 3;
		J(2, 1) = 4;
		J(2, 2) = 5;
		J(2, 3) = 6;
		
		Matrix JT(3,2);
		JT = transpose(J);

		Matrix J_3(2);
		J_3 = slice(J, 1, 2, 1, 2);
		J_3.transpose();

		for(int i = 1; i <= 2; i++)
		{
			for(int j = 1; j <= 3; j++)
			{
				assert(J(i, j) == JT(j, i));
			}
		}

		assert(J_3(1,1) == 1);
		assert(J_3(1,2) == 4);
		assert(J_3(2,1) == 2);
		assert(J_3(2,2) == 5);
	}
	catch (Exception &ex)
	{
		cout << "FAIL!\n";
		flush(cout);
		ex.DebugPrint();
		return;
	}
	cout << "Success!\n";
	flush(cout);

	return;
}


// we don't test big factorisations since they don't have a convergence criteria
// but rather are exact
void test_factorisations()
{
	cout << "Testing factorisations and decompositions... ";
	flush(cout);

	try
	{
		Matrix J(3);
		J(1, 1) = 1;
		J(1, 2) = 2;
		J(1, 3) = 3;
		J(2, 1) = 6;
		J(2, 2) = 5;
		J(2, 3) = 4;
		J(3, 1) = 19;
		J(3, 2) = 8;
		J(3, 3) = 17;
		
		Matrix L(3);
		Matrix U(3);
		Matrix P(3,"eye");
		L = lu(J,U);

		Matrix LU(3);
		LU = L*U;

		for(int i = 1; i <= 3; i++)
		{
			for(int j = 1; j <= 3; j++)
			{
				assert(fabs(J(i, j) - LU(i, j)) < 1E-10);
				if(i == j)
				{
					assert(L(i, j) == 1.0);
				}
				if(i < j)
				{
					assert(L(i, j) == 0.0);
				}
				if(i > j)
				{
					assert(U(i, j) == 0.0);
				}
			}
		}

		Matrix J2(3,4);
		J2(1, 1) = 1;
		J2(1, 2) = 2;
		J2(1, 3) = 3;
		J2(1, 4) = 4;
		J2(2, 1) = 4;
		J2(2, 2) = 5;
		J2(2, 3) = 6;
		J2(2, 4) = 7;
		J2(3, 1) = 7;
		J2(3, 2) = 8;
		J2(3, 3) = 9;
		J2(3, 4) = 10;

		Matrix U2(3,4);

		L = lu(J2,U2,P);
		Matrix PJ2(3,4);
		Matrix LU2(3,4);
		PJ2 = P*J2;
		LU2 = L*U2;

		for(int i = 1; i <= 3; i++)
		{
			for(int j = 1; j <= 4; j++)
			{
				assert(fabs(PJ2(i, j) - LU2(i, j)) < 1E-10);
				if(i == j && j != 4)
				{
					assert(L(i, j) == 1.0);
				}
				if(i < j && j != 4)
				{
					assert(L(i, j) == 0.0);
				}
				if(i > j)
				{
					assert(U2(i, j) == 0.0);
				}
			}
		}

		// we actually want to solve correctly here,
		// so explicitly form Q
		bool explicit_Q = 1;
		Matrix Q(3);
		Matrix R(J);
		Q = qr(R, explicit_Q); // qr overwrites J
		Matrix QR(3);
		QR = Q*R;


		for(int i = 1; i <= 3; i++)
		{
			for(int j = 1; j <= 3; j++)
			{
				assert(fabs(QR(i, j) - J(i, j)) < 1E-10);
				if(i > j)
				{
					assert(R(i, j) == 0.0);
				}
			}
		}


	}
	catch (Exception &ex)
	{
		cout << "FAIL!\n";
		flush(cout);
		ex.DebugPrint();
		return;
	}
	cout << "Success!\n";
	flush(cout);

	return;
}


void test_det()
{
	cout << "Testing determinant calculation... ";
	flush(cout);

	try
	{
		Matrix J(3);
		J(1, 1) = 1;
		J(1, 2) = 2;
		J(1, 3) = 3;
		J(2, 1) = 6;
		J(2, 2) = 5;
		J(2, 3) = 4;
		J(3, 1) = 19;
		J(3, 2) = 8;
		J(3, 3) = 17;
		assert(det(J) == -140);

		J(2, 1) = 1;
		J(2, 2) = 2;
		J(2, 3) = 3;
		J(1, 1) = 6;
		J(1, 2) = 5;
		J(1, 3) = 4;
		J(3, 1) = 19;
		J(3, 2) = 8;
		J(3, 3) = 17;
		assert(det(J) == 140);

	}
	catch (Exception &ex)
	{
		cout << "FAIL!\n";
		flush(cout);
		ex.DebugPrint();
		return;
	}
	cout << "Success!\n";
	flush(cout);

	return;
}


void test_triangular_solvers()
{
	cout << "Testing triangular linear solvers... ";
	flush(cout);

	try
	{
		Vector b(3);
		b(1) = 1;
		b(2) = 2;
		b(3) = 3;

		Matrix J(3);
		J(1, 1) = 1;
		J(1, 2) = 2;
		J(1, 3) = 3;
		J(2, 1) = 0;
		J(2, 2) = 4;
		J(2, 3) = 5;
		J(3, 1) = 0;
		J(3, 2) = 0;
		J(3, 3) = 6;

		Vector x(3);

		// compare backsub with MATLAB
		x = backsub(J, b);
		assert(x(1) == -0.25);
		assert(x(2) == -0.125);
		assert(x(3) == 0.5);

		J(1, 1) = 1;
		J(1, 2) = 0;
		J(1, 3) = 0;
		J(2, 1) = 2;
		J(2, 2) = 3;
		J(2, 3) = 0;
		J(3, 1) = 4;
		J(3, 2) = 5;
		J(3, 3) = 6;

		// compare forward with MATLAB
		x = forwardsub(J, b);
		assert(x(1) == 1.0);
		assert(x(2) == 0.0);
		assert(x(3) == -1.0/6.0);
	}
	catch (Exception &ex)
	{
		cout << "FAIL!\n";
		flush(cout);
		ex.DebugPrint();
		return;
	}
	cout << "Success!\n";
	flush(cout);

	return;
}


void test_iterative_solvers()
{
	cout << "Testing iterative linear solvers... \n";
	flush(cout);

	try
	{
		bool verbose = 1;

		// first test that the functions can be called
		Matrix A(3);
		
		A(1,1) = 4;
		A(1,2) = 1;
		A(1,3) = 2;
		A(2,1) = 1;
		A(2,2) = 3;
		A(2,3) = 2;	
		A(3,1) = 2;
		A(3,2) = 2;
		A(3,3) = 4;

		Vector b(3);
		b(1) = 1;
		b(2) = 3;
		b(3) = 5;

		Vector x(3); // to hold solution

		// test iterations without an initial point
		// since they just pass a zero initial vector
		x = gmres(A, b, verbose);
		assert(x(1) + 0.5 < 1E-6);
		assert(x(2) - 0.25 < 1E-6);
		assert(x(3) - 1.375 < 1E-6);

		x = cg(A, b, verbose);
		assert(x(1) + 0.5 < 1E-6);
		assert(x(2) - 0.25 < 1E-6);
		assert(x(3) - 1.375 < 1E-6);

		// now test that they work iteratively
		// the test Matrix and Vectors were built in MATLAB
		// and code was generated for initialisation in Python
		Matrix A2(32);
		A2 = build_test_mat();

		Vector b2(32);
		b2 = build_test_vec();

		Vector xtest(32);
		xtest = build_test_soln();

		Vector x2(32);
		x2 = gmres(A2, b2, verbose);
		assert(fabs((A2 * x2).norm() - b2.norm()) < 1E-12);
		// less precision below because we are limited by the initialisation values
		// of xtest and how precise we write the csvs from MATLAB :)
		assert(fabs((xtest).norm() - x2.norm()) < 1E-5);

		x2 = cg(A2, b2, verbose);
		assert(fabs((A2 * x2).norm() - b2.norm()) < 1E-12);
		assert(fabs(xtest.norm() - x2.norm()) < 1E-5);


	}
	catch (Exception &ex)
	{
		cout << "FAIL!\n";
		flush(cout);
		ex.DebugPrint();
		return;
	}
	cout << "Success!\n";
	flush(cout);

	return;
}


// BEWARE! yucky hard-coded numbers
Vector build_test_vec()
{
	Vector b(32);

	b(1) = 17.61525;
	b(2) = 61.35886;
	b(3) = 26.52781;
	b(4) = 47.42464;
	b(5) = 47.80742;
	b(6) = 32.00974;
	b(7) = 51.47232;
	b(8) = 40.85589;
	b(9) = 57.66453;
	b(10) = 5.88511;
	b(11) = 14.06886;
	b(12) = 10.71291;
	b(13) = 28.58424;
	b(14) = 10.35378;
	b(15) = 42.84752;
	b(16) = 57.79759;
	b(17) = 8.699403;
	b(18) = 6.844801;
	b(19) = 3.314149;
	b(20) = 57.77131;
	b(21) = 34.23105;
	b(22) = 63.02474;
	b(23) = 49.00163;
	b(24) = 56.72047;
	b(25) = 27.35786;
	b(26) = 61.27361;
	b(27) = 21.35585;
	b(28) = 50.42265;
	b(29) = 48.12957;
	b(30) = 58.55449;
	b(31) = 20.74123;
	b(32) = 54.10169;

	return b;
}


Vector build_test_soln()
{
	Vector xtest(32);

	xtest(1) = -0.543337;
	xtest(2) = 0.8642605;
	xtest(3) = 0.07356812;
	xtest(4) = 0.3923438;
	xtest(5) = 0.2729096;
	xtest(6) = -0.05085879;
	xtest(7) = 0.3975401;
	xtest(8) = 0.393314;
	xtest(9) = 0.7880766;
	xtest(10) = -0.6238471;
	xtest(11) = -0.5006098;
	xtest(12) = -0.6255704;
	xtest(13) = -0.1483315;
	xtest(14) = -0.6179134;
	xtest(15) = 0.3458447;
	xtest(16) = 0.5470878;
	xtest(17) = -0.7256269;
	xtest(18) = -0.7679861;
	xtest(19) = -0.7271971;
	xtest(20) = 0.5856957;
	xtest(21) = 0.2039201;
	xtest(22) = 0.681425;
	xtest(23) = 0.3980699;
	xtest(24) = 0.7242048;
	xtest(25) = -0.141964;
	xtest(26) = 0.8220481;
	xtest(27) = -0.4375674;
	xtest(28) = 0.581416;
	xtest(29) = 0.443743;
	xtest(30) = 0.820181;
	xtest(31) = -0.3733699;
	xtest(32) = 0.6320287;

	return xtest;
}


Matrix build_test_mat()
{
	Matrix A(32);

	A(1, 1) = 44.23582;
	A(1, 2) = 9.245092;
	A(1, 3) = 6.649603;
	A(1, 4) = 8.51674;
	A(1, 5) = 10.11802;
	A(1, 6) = 8.956281;
	A(1, 7) = 10.39055;
	A(1, 8) = 7.293228;
	A(1, 9) = 8.911456;
	A(1, 10) = 7.382321;
	A(1, 11) = 9.091552;
	A(1, 12) = 8.534986;
	A(1, 13) = 9.129124;
	A(1, 14) = 8.539641;
	A(1, 15) = 8.374846;
	A(1, 16) = 11.05993;
	A(1, 17) = 9.113927;
	A(1, 18) = 9.069862;
	A(1, 19) = 8.171354;
	A(1, 20) = 10.4536;
	A(1, 21) = 7.026091;
	A(1, 22) = 11.42303;
	A(1, 23) = 9.473822;
	A(1, 24) = 7.802992;
	A(1, 25) = 8.421722;
	A(1, 26) = 9.0935;
	A(1, 27) = 9.107638;
	A(1, 28) = 9.442419;
	A(1, 29) = 9.134119;
	A(1, 30) = 8.025165;
	A(1, 31) = 9.197884;
	A(1, 32) = 8.631121;
	A(2, 1) = 9.245092;
	A(2, 2) = 42.78328;
	A(2, 3) = 5.982738;
	A(2, 4) = 8.28829;
	A(2, 5) = 9.416464;
	A(2, 6) = 7.374202;
	A(2, 7) = 8.982329;
	A(2, 8) = 6.435612;
	A(2, 9) = 9.084778;
	A(2, 10) = 7.619578;
	A(2, 11) = 7.929278;
	A(2, 12) = 8.648113;
	A(2, 13) = 8.741815;
	A(2, 14) = 8.082172;
	A(2, 15) = 7.808906;
	A(2, 16) = 10.43379;
	A(2, 17) = 8.698342;
	A(2, 18) = 8.172884;
	A(2, 19) = 8.154764;
	A(2, 20) = 9.102691;
	A(2, 21) = 7.011811;
	A(2, 22) = 10.44021;
	A(2, 23) = 8.734475;
	A(2, 24) = 7.359481;
	A(2, 25) = 7.521159;
	A(2, 26) = 8.532229;
	A(2, 27) = 9.692117;
	A(2, 28) = 7.942288;
	A(2, 29) = 8.72569;
	A(2, 30) = 6.940404;
	A(2, 31) = 8.204024;
	A(2, 32) = 8.361619;
	A(3, 1) = 6.649603;
	A(3, 2) = 5.982738;
	A(3, 3) = 38.45913;
	A(3, 4) = 6.064802;
	A(3, 5) = 6.677219;
	A(3, 6) = 5.484416;
	A(3, 7) = 6.857371;
	A(3, 8) = 5.114658;
	A(3, 9) = 6.276779;
	A(3, 10) = 6.021932;
	A(3, 11) = 5.728537;
	A(3, 12) = 6.173596;
	A(3, 13) = 5.591511;
	A(3, 14) = 5.924653;
	A(3, 15) = 6.281531;
	A(3, 16) = 7.358663;
	A(3, 17) = 5.875128;
	A(3, 18) = 5.075829;
	A(3, 19) = 5.822643;
	A(3, 20) = 5.956427;
	A(3, 21) = 4.780523;
	A(3, 22) = 7.531752;
	A(3, 23) = 7.020379;
	A(3, 24) = 5.543248;
	A(3, 25) = 5.718482;
	A(3, 26) = 5.395215;
	A(3, 27) = 6.688468;
	A(3, 28) = 6.211776;
	A(3, 29) = 6.30242;
	A(3, 30) = 5.29562;
	A(3, 31) = 5.810936;
	A(3, 32) = 6.367104;
	A(4, 1) = 8.51674;
	A(4, 2) = 8.28829;
	A(4, 3) = 6.064802;
	A(4, 4) = 42.41075;
	A(4, 5) = 9.415558;
	A(4, 6) = 8.163077;
	A(4, 7) = 9.041188;
	A(4, 8) = 7.355028;
	A(4, 9) = 8.564954;
	A(4, 10) = 7.791804;
	A(4, 11) = 8.243576;
	A(4, 12) = 8.208992;
	A(4, 13) = 7.725002;
	A(4, 14) = 8.105379;
	A(4, 15) = 7.031554;
	A(4, 16) = 9.963557;
	A(4, 17) = 7.029001;
	A(4, 18) = 7.571491;
	A(4, 19) = 6.893831;
	A(4, 20) = 9.295841;
	A(4, 21) = 7.409582;
	A(4, 22) = 9.322583;
	A(4, 23) = 8.535307;
	A(4, 24) = 7.756538;
	A(4, 25) = 6.917626;
	A(4, 26) = 7.476358;
	A(4, 27) = 9.164117;
	A(4, 28) = 7.526169;
	A(4, 29) = 9.105786;
	A(4, 30) = 7.531226;
	A(4, 31) = 8.335414;
	A(4, 32) = 9.13942;
	A(5, 1) = 10.11802;
	A(5, 2) = 9.416464;
	A(5, 3) = 6.677219;
	A(5, 4) = 9.415558;
	A(5, 5) = 46.04899;
	A(5, 6) = 9.110519;
	A(5, 7) = 11.0016;
	A(5, 8) = 7.663898;
	A(5, 9) = 9.848471;
	A(5, 10) = 8.460583;
	A(5, 11) = 9.120053;
	A(5, 12) = 8.56692;
	A(5, 13) = 10.50126;
	A(5, 14) = 9.337581;
	A(5, 15) = 9.106943;
	A(5, 16) = 11.58421;
	A(5, 17) = 9.466877;
	A(5, 18) = 9.597993;
	A(5, 19) = 9.374222;
	A(5, 20) = 10.20244;
	A(5, 21) = 8.308143;
	A(5, 22) = 10.78901;
	A(5, 23) = 10.74663;
	A(5, 24) = 9.37942;
	A(5, 25) = 9.058475;
	A(5, 26) = 10.3263;
	A(5, 27) = 10.59962;
	A(5, 28) = 9.52149;
	A(5, 29) = 9.969017;
	A(5, 30) = 8.811292;
	A(5, 31) = 9.604378;
	A(5, 32) = 8.815245;
	A(6, 1) = 8.956281;
	A(6, 2) = 7.374202;
	A(6, 3) = 5.484416;
	A(6, 4) = 8.163077;
	A(6, 5) = 9.110519;
	A(6, 6) = 41.8224;
	A(6, 7) = 8.687003;
	A(6, 8) = 7.292682;
	A(6, 9) = 7.812928;
	A(6, 10) = 6.066884;
	A(6, 11) = 7.722062;
	A(6, 12) = 7.021276;
	A(6, 13) = 8.375001;
	A(6, 14) = 7.128219;
	A(6, 15) = 7.294345;
	A(6, 16) = 9.86902;
	A(6, 17) = 7.554961;
	A(6, 18) = 7.249692;
	A(6, 19) = 6.955856;
	A(6, 20) = 9.467061;
	A(6, 21) = 6.778485;
	A(6, 22) = 9.55513;
	A(6, 23) = 8.789953;
	A(6, 24) = 7.713491;
	A(6, 25) = 7.553156;
	A(6, 26) = 8.095305;
	A(6, 27) = 7.494553;
	A(6, 28) = 7.376231;
	A(6, 29) = 8.644343;
	A(6, 30) = 6.765675;
	A(6, 31) = 8.372642;
	A(6, 32) = 7.215221;
	A(7, 1) = 10.39055;
	A(7, 2) = 8.982329;
	A(7, 3) = 6.857371;
	A(7, 4) = 9.041188;
	A(7, 5) = 11.0016;
	A(7, 6) = 8.687003;
	A(7, 7) = 45.85355;
	A(7, 8) = 7.329371;
	A(7, 9) = 10.00191;
	A(7, 10) = 8.602264;
	A(7, 11) = 9.452686;
	A(7, 12) = 8.927584;
	A(7, 13) = 8.970857;
	A(7, 14) = 9.667432;
	A(7, 15) = 8.645294;
	A(7, 16) = 11.87633;
	A(7, 17) = 8.951361;
	A(7, 18) = 8.583562;
	A(7, 19) = 8.772589;
	A(7, 20) = 10.43226;
	A(7, 21) = 7.149815;
	A(7, 22) = 10.89035;
	A(7, 23) = 10.49166;
	A(7, 24) = 8.705566;
	A(7, 25) = 9.262994;
	A(7, 26) = 9.333155;
	A(7, 27) = 8.90611;
	A(7, 28) = 10.33834;
	A(7, 29) = 8.989779;
	A(7, 30) = 8.447423;
	A(7, 31) = 9.762229;
	A(7, 32) = 8.705503;
	A(8, 1) = 7.293228;
	A(8, 2) = 6.435612;
	A(8, 3) = 5.114658;
	A(8, 4) = 7.355028;
	A(8, 5) = 7.663898;
	A(8, 6) = 7.292682;
	A(8, 7) = 7.329371;
	A(8, 8) = 40.05005;
	A(8, 9) = 7.687752;
	A(8, 10) = 6.006975;
	A(8, 11) = 7.035912;
	A(8, 12) = 7.036352;
	A(8, 13) = 6.767463;
	A(8, 14) = 6.35503;
	A(8, 15) = 6.962544;
	A(8, 16) = 7.935115;
	A(8, 17) = 6.986355;
	A(8, 18) = 6.204019;
	A(8, 19) = 6.200157;
	A(8, 20) = 7.893106;
	A(8, 21) = 5.937622;
	A(8, 22) = 8.228448;
	A(8, 23) = 7.095266;
	A(8, 24) = 7.108769;
	A(8, 25) = 6.249789;
	A(8, 26) = 6.456123;
	A(8, 27) = 6.674888;
	A(8, 28) = 6.155236;
	A(8, 29) = 7.478709;
	A(8, 30) = 5.490128;
	A(8, 31) = 7.22629;
	A(8, 32) = 6.534867;
	A(9, 1) = 8.911456;
	A(9, 2) = 9.084778;
	A(9, 3) = 6.276779;
	A(9, 4) = 8.564954;
	A(9, 5) = 9.848471;
	A(9, 6) = 7.812928;
	A(9, 7) = 10.00191;
	A(9, 8) = 7.687752;
	A(9, 9) = 43.34357;
	A(9, 10) = 8.400121;
	A(9, 11) = 8.369882;
	A(9, 12) = 9.018649;
	A(9, 13) = 8.921093;
	A(9, 14) = 8.646107;
	A(9, 15) = 7.955535;
	A(9, 16) = 10.42008;
	A(9, 17) = 8.728746;
	A(9, 18) = 7.982992;
	A(9, 19) = 8.219064;
	A(9, 20) = 9.144565;
	A(9, 21) = 6.645207;
	A(9, 22) = 9.178903;
	A(9, 23) = 8.999753;
	A(9, 24) = 7.558445;
	A(9, 25) = 7.555898;
	A(9, 26) = 7.915104;
	A(9, 27) = 8.719343;
	A(9, 28) = 8.830105;
	A(9, 29) = 9.526395;
	A(9, 30) = 5.566744;
	A(9, 31) = 8.854994;
	A(9, 32) = 7.674169;
	A(10, 1) = 7.382321;
	A(10, 2) = 7.619578;
	A(10, 3) = 6.021932;
	A(10, 4) = 7.791804;
	A(10, 5) = 8.460583;
	A(10, 6) = 6.066884;
	A(10, 7) = 8.602264;
	A(10, 8) = 6.006975;
	A(10, 9) = 8.400121;
	A(10, 10) = 41.79683;
	A(10, 11) = 8.039888;
	A(10, 12) = 7.202514;
	A(10, 13) = 6.965702;
	A(10, 14) = 7.725535;
	A(10, 15) = 6.45108;
	A(10, 16) = 9.372269;
	A(10, 17) = 6.125163;
	A(10, 18) = 6.904225;
	A(10, 19) = 7.622837;
	A(10, 20) = 7.341396;
	A(10, 21) = 5.895427;
	A(10, 22) = 8.227939;
	A(10, 23) = 8.095408;
	A(10, 24) = 5.737969;
	A(10, 25) = 6.739049;
	A(10, 26) = 6.105246;
	A(10, 27) = 8.051716;
	A(10, 28) = 7.516197;
	A(10, 29) = 8.150392;
	A(10, 30) = 5.744286;
	A(10, 31) = 7.489577;
	A(10, 32) = 7.6337;
	A(11, 1) = 9.091552;
	A(11, 2) = 7.929278;
	A(11, 3) = 5.728537;
	A(11, 4) = 8.243576;
	A(11, 5) = 9.120053;
	A(11, 6) = 7.722062;
	A(11, 7) = 9.452686;
	A(11, 8) = 7.035912;
	A(11, 9) = 8.369882;
	A(11, 10) = 8.039888;
	A(11, 11) = 41.60612;
	A(11, 12) = 7.499724;
	A(11, 13) = 7.729023;
	A(11, 14) = 7.438088;
	A(11, 15) = 7.019857;
	A(11, 16) = 9.989879;
	A(11, 17) = 7.560997;
	A(11, 18) = 7.641013;
	A(11, 19) = 7.181202;
	A(11, 20) = 8.430634;
	A(11, 21) = 6.656947;
	A(11, 22) = 9.157767;
	A(11, 23) = 8.828074;
	A(11, 24) = 6.928482;
	A(11, 25) = 7.316482;
	A(11, 26) = 7.170356;
	A(11, 27) = 7.771206;
	A(11, 28) = 7.801772;
	A(11, 29) = 8.555453;
	A(11, 30) = 6.612849;
	A(11, 31) = 7.937045;
	A(11, 32) = 7.798531;
	A(12, 1) = 8.534986;
	A(12, 2) = 8.648113;
	A(12, 3) = 6.173596;
	A(12, 4) = 8.208992;
	A(12, 5) = 8.56692;
	A(12, 6) = 7.021276;
	A(12, 7) = 8.927584;
	A(12, 8) = 7.036352;
	A(12, 9) = 9.018649;
	A(12, 10) = 7.202514;
	A(12, 11) = 7.499724;
	A(12, 12) = 42.6957;
	A(12, 13) = 8.403827;
	A(12, 14) = 7.874285;
	A(12, 15) = 8.031177;
	A(12, 16) = 9.520031;
	A(12, 17) = 8.204215;
	A(12, 18) = 7.657;
	A(12, 19) = 7.034878;
	A(12, 20) = 8.198113;
	A(12, 21) = 6.78881;
	A(12, 22) = 10.71875;
	A(12, 23) = 8.28483;
	A(12, 24) = 7.845776;
	A(12, 25) = 7.022617;
	A(12, 26) = 6.983741;
	A(12, 27) = 8.927778;
	A(12, 28) = 8.285803;
	A(12, 29) = 7.814824;
	A(12, 30) = 6.1109;
	A(12, 31) = 8.875031;
	A(12, 32) = 8.762041;
	A(13, 1) = 9.129124;
	A(13, 2) = 8.741815;
	A(13, 3) = 5.591511;
	A(13, 4) = 7.725002;
	A(13, 5) = 10.50126;
	A(13, 6) = 8.375001;
	A(13, 7) = 8.970857;
	A(13, 8) = 6.767463;
	A(13, 9) = 8.921093;
	A(13, 10) = 6.965702;
	A(13, 11) = 7.729023;
	A(13, 12) = 8.403827;
	A(13, 13) = 43.18258;
	A(13, 14) = 7.612771;
	A(13, 15) = 8.091345;
	A(13, 16) = 9.841947;
	A(13, 17) = 8.880984;
	A(13, 18) = 8.549732;
	A(13, 19) = 7.945046;
	A(13, 20) = 8.95281;
	A(13, 21) = 7.368425;
	A(13, 22) = 10.17831;
	A(13, 23) = 9.349998;
	A(13, 24) = 8.328472;
	A(13, 25) = 7.320691;
	A(13, 26) = 8.708432;
	A(13, 27) = 9.071458;
	A(13, 28) = 8.533298;
	A(13, 29) = 9.012588;
	A(13, 30) = 6.767958;
	A(13, 31) = 8.550621;
	A(13, 32) = 8.117458;
	A(14, 1) = 8.539641;
	A(14, 2) = 8.082172;
	A(14, 3) = 5.924653;
	A(14, 4) = 8.105379;
	A(14, 5) = 9.337581;
	A(14, 6) = 7.128219;
	A(14, 7) = 9.667432;
	A(14, 8) = 6.35503;
	A(14, 9) = 8.646107;
	A(14, 10) = 7.725535;
	A(14, 11) = 7.438088;
	A(14, 12) = 7.874285;
	A(14, 13) = 7.612771;
	A(14, 14) = 41.66961;
	A(14, 15) = 7.062624;
	A(14, 16) = 9.434391;
	A(14, 17) = 7.677468;
	A(14, 18) = 7.173582;
	A(14, 19) = 7.997362;
	A(14, 20) = 8.257515;
	A(14, 21) = 6.592223;
	A(14, 22) = 8.789155;
	A(14, 23) = 8.49129;
	A(14, 24) = 6.906992;
	A(14, 25) = 6.826281;
	A(14, 26) = 8.074467;
	A(14, 27) = 8.183154;
	A(14, 28) = 7.966928;
	A(14, 29) = 8.061718;
	A(14, 30) = 6.672546;
	A(14, 31) = 7.085254;
	A(14, 32) = 7.930122;
	A(15, 1) = 8.374846;
	A(15, 2) = 7.808906;
	A(15, 3) = 6.281531;
	A(15, 4) = 7.031554;
	A(15, 5) = 9.106943;
	A(15, 6) = 7.294345;
	A(15, 7) = 8.645294;
	A(15, 8) = 6.962544;
	A(15, 9) = 7.955535;
	A(15, 10) = 6.45108;
	A(15, 11) = 7.019857;
	A(15, 12) = 8.031177;
	A(15, 13) = 8.091345;
	A(15, 14) = 7.062624;
	A(15, 15) = 41.46988;
	A(15, 16) = 9.371882;
	A(15, 17) = 8.087656;
	A(15, 18) = 7.178695;
	A(15, 19) = 7.69663;
	A(15, 20) = 8.207889;
	A(15, 21) = 6.159832;
	A(15, 22) = 10.48912;
	A(15, 23) = 8.478971;
	A(15, 24) = 7.962222;
	A(15, 25) = 7.913155;
	A(15, 26) = 7.705297;
	A(15, 27) = 7.734393;
	A(15, 28) = 7.216666;
	A(15, 29) = 8.057176;
	A(15, 30) = 5.870655;
	A(15, 31) = 7.890541;
	A(15, 32) = 7.414451;
	A(16, 1) = 11.05993;
	A(16, 2) = 10.43379;
	A(16, 3) = 7.358663;
	A(16, 4) = 9.963557;
	A(16, 5) = 11.58421;
	A(16, 6) = 9.86902;
	A(16, 7) = 11.87633;
	A(16, 8) = 7.935115;
	A(16, 9) = 10.42008;
	A(16, 10) = 9.372269;
	A(16, 11) = 9.989879;
	A(16, 12) = 9.520031;
	A(16, 13) = 9.841947;
	A(16, 14) = 9.434391;
	A(16, 15) = 9.371882;
	A(16, 16) = 45.71568;
	A(16, 17) = 9.118902;
	A(16, 18) = 10.05577;
	A(16, 19) = 9.303243;
	A(16, 20) = 10.92504;
	A(16, 21) = 8.601168;
	A(16, 22) = 12.15409;
	A(16, 23) = 10.67833;
	A(16, 24) = 8.315734;
	A(16, 25) = 9.747427;
	A(16, 26) = 9.661396;
	A(16, 27) = 10.50385;
	A(16, 28) = 10.15614;
	A(16, 29) = 10.52856;
	A(16, 30) = 8.399444;
	A(16, 31) = 9.488288;
	A(16, 32) = 9.458236;
	A(17, 1) = 9.113927;
	A(17, 2) = 8.698342;
	A(17, 3) = 5.875128;
	A(17, 4) = 7.029001;
	A(17, 5) = 9.466877;
	A(17, 6) = 7.554961;
	A(17, 7) = 8.951361;
	A(17, 8) = 6.986355;
	A(17, 9) = 8.728746;
	A(17, 10) = 6.125163;
	A(17, 11) = 7.560997;
	A(17, 12) = 8.204215;
	A(17, 13) = 8.880984;
	A(17, 14) = 7.677468;
	A(17, 15) = 8.087656;
	A(17, 16) = 9.118902;
	A(17, 17) = 42.23637;
	A(17, 18) = 7.304122;
	A(17, 19) = 6.995418;
	A(17, 20) = 8.028089;
	A(17, 21) = 6.814038;
	A(17, 22) = 9.781071;
	A(17, 23) = 8.816444;
	A(17, 24) = 8.104349;
	A(17, 25) = 7.359839;
	A(17, 26) = 8.698833;
	A(17, 27) = 8.458376;
	A(17, 28) = 8.151227;
	A(17, 29) = 7.824457;
	A(17, 30) = 6.353978;
	A(17, 31) = 8.320092;
	A(17, 32) = 8.007373;
	A(18, 1) = 9.069862;
	A(18, 2) = 8.172884;
	A(18, 3) = 5.075829;
	A(18, 4) = 7.571491;
	A(18, 5) = 9.597993;
	A(18, 6) = 7.249692;
	A(18, 7) = 8.583562;
	A(18, 8) = 6.204019;
	A(18, 9) = 7.982992;
	A(18, 10) = 6.904225;
	A(18, 11) = 7.641013;
	A(18, 12) = 7.657;
	A(18, 13) = 8.549732;
	A(18, 14) = 7.173582;
	A(18, 15) = 7.178695;
	A(18, 16) = 10.05577;
	A(18, 17) = 7.304122;
	A(18, 18) = 42.64573;
	A(18, 19) = 6.454747;
	A(18, 20) = 8.628462;
	A(18, 21) = 6.199064;
	A(18, 22) = 10.18262;
	A(18, 23) = 8.458736;
	A(18, 24) = 6.945464;
	A(18, 25) = 7.284826;
	A(18, 26) = 8.090492;
	A(18, 27) = 8.94317;
	A(18, 28) = 8.657142;
	A(18, 29) = 8.664519;
	A(18, 30) = 7.552076;
	A(18, 31) = 7.330491;
	A(18, 32) = 7.377906;
	A(19, 1) = 8.171354;
	A(19, 2) = 8.154764;
	A(19, 3) = 5.822643;
	A(19, 4) = 6.893831;
	A(19, 5) = 9.374222;
	A(19, 6) = 6.955856;
	A(19, 7) = 8.772589;
	A(19, 8) = 6.200157;
	A(19, 9) = 8.219064;
	A(19, 10) = 7.622837;
	A(19, 11) = 7.181202;
	A(19, 12) = 7.034878;
	A(19, 13) = 7.945046;
	A(19, 14) = 7.997362;
	A(19, 15) = 7.69663;
	A(19, 16) = 9.303243;
	A(19, 17) = 6.995418;
	A(19, 18) = 6.454747;
	A(19, 19) = 41.94908;
	A(19, 20) = 7.963725;
	A(19, 21) = 6.214614;
	A(19, 22) = 8.458713;
	A(19, 23) = 8.154207;
	A(19, 24) = 6.871617;
	A(19, 25) = 7.357392;
	A(19, 26) = 7.528459;
	A(19, 27) = 7.624521;
	A(19, 28) = 6.722559;
	A(19, 29) = 8.055284;
	A(19, 30) = 5.277095;
	A(19, 31) = 7.429545;
	A(19, 32) = 6.107671;
	A(20, 1) = 10.4536;
	A(20, 2) = 9.102691;
	A(20, 3) = 5.956427;
	A(20, 4) = 9.295841;
	A(20, 5) = 10.20244;
	A(20, 6) = 9.467061;
	A(20, 7) = 10.43226;
	A(20, 8) = 7.893106;
	A(20, 9) = 9.144565;
	A(20, 10) = 7.341396;
	A(20, 11) = 8.430634;
	A(20, 12) = 8.198113;
	A(20, 13) = 8.95281;
	A(20, 14) = 8.257515;
	A(20, 15) = 8.207889;
	A(20, 16) = 10.92504;
	A(20, 17) = 8.028089;
	A(20, 18) = 8.628462;
	A(20, 19) = 7.963725;
	A(20, 20) = 43.88015;
	A(20, 21) = 7.185614;
	A(20, 22) = 10.86741;
	A(20, 23) = 9.56881;
	A(20, 24) = 8.486403;
	A(20, 25) = 8.834362;
	A(20, 26) = 9.238183;
	A(20, 27) = 9.052507;
	A(20, 28) = 8.407161;
	A(20, 29) = 9.306524;
	A(20, 30) = 8.352698;
	A(20, 31) = 9.509018;
	A(20, 32) = 8.387606;
	A(21, 1) = 7.026091;
	A(21, 2) = 7.011811;
	A(21, 3) = 4.780523;
	A(21, 4) = 7.409582;
	A(21, 5) = 8.308143;
	A(21, 6) = 6.778485;
	A(21, 7) = 7.149815;
	A(21, 8) = 5.937622;
	A(21, 9) = 6.645207;
	A(21, 10) = 5.895427;
	A(21, 11) = 6.656947;
	A(21, 12) = 6.78881;
	A(21, 13) = 7.368425;
	A(21, 14) = 6.592223;
	A(21, 15) = 6.159832;
	A(21, 16) = 8.601168;
	A(21, 17) = 6.814038;
	A(21, 18) = 6.199064;
	A(21, 19) = 6.214614;
	A(21, 20) = 7.185614;
	A(21, 21) = 39.65204;
	A(21, 22) = 7.73244;
	A(21, 23) = 7.239654;
	A(21, 24) = 6.103824;
	A(21, 25) = 6.095469;
	A(21, 26) = 6.83822;
	A(21, 27) = 7.735146;
	A(21, 28) = 6.095053;
	A(21, 29) = 6.755796;
	A(21, 30) = 6.11867;
	A(21, 31) = 6.287632;
	A(21, 32) = 7.523602;
	A(22, 1) = 11.42303;
	A(22, 2) = 10.44021;
	A(22, 3) = 7.531752;
	A(22, 4) = 9.322583;
	A(22, 5) = 10.78901;
	A(22, 6) = 9.55513;
	A(22, 7) = 10.89035;
	A(22, 8) = 8.228448;
	A(22, 9) = 9.178903;
	A(22, 10) = 8.227939;
	A(22, 11) = 9.157767;
	A(22, 12) = 10.71875;
	A(22, 13) = 10.17831;
	A(22, 14) = 8.789155;
	A(22, 15) = 10.48912;
	A(22, 16) = 12.15409;
	A(22, 17) = 9.781071;
	A(22, 18) = 10.18262;
	A(22, 19) = 8.458713;
	A(22, 20) = 10.86741;
	A(22, 21) = 7.73244;
	A(22, 22) = 46.60693;
	A(22, 23) = 10.39502;
	A(22, 24) = 9.77928;
	A(22, 25) = 9.771961;
	A(22, 26) = 9.633077;
	A(22, 27) = 10.41915;
	A(22, 28) = 9.951712;
	A(22, 29) = 9.582796;
	A(22, 30) = 8.793287;
	A(22, 31) = 10.35086;
	A(22, 32) = 10.46545;
	A(23, 1) = 9.473822;
	A(23, 2) = 8.734475;
	A(23, 3) = 7.020379;
	A(23, 4) = 8.535307;
	A(23, 5) = 10.74663;
	A(23, 6) = 8.789953;
	A(23, 7) = 10.49166;
	A(23, 8) = 7.095266;
	A(23, 9) = 8.999753;
	A(23, 10) = 8.095408;
	A(23, 11) = 8.828074;
	A(23, 12) = 8.28483;
	A(23, 13) = 9.349998;
	A(23, 14) = 8.49129;
	A(23, 15) = 8.478971;
	A(23, 16) = 10.67833;
	A(23, 17) = 8.816444;
	A(23, 18) = 8.458736;
	A(23, 19) = 8.154207;
	A(23, 20) = 9.56881;
	A(23, 21) = 7.239654;
	A(23, 22) = 10.39502;
	A(23, 23) = 43.62555;
	A(23, 24) = 8.144623;
	A(23, 25) = 8.504797;
	A(23, 26) = 8.80019;
	A(23, 27) = 9.413723;
	A(23, 28) = 9.23968;
	A(23, 29) = 9.801824;
	A(23, 30) = 7.998706;
	A(23, 31) = 9.271368;
	A(23, 32) = 8.502528;
	A(24, 1) = 7.802992;
	A(24, 2) = 7.359481;
	A(24, 3) = 5.543248;
	A(24, 4) = 7.756538;
	A(24, 5) = 9.37942;
	A(24, 6) = 7.713491;
	A(24, 7) = 8.705566;
	A(24, 8) = 7.108769;
	A(24, 9) = 7.558445;
	A(24, 10) = 5.737969;
	A(24, 11) = 6.928482;
	A(24, 12) = 7.845776;
	A(24, 13) = 8.328472;
	A(24, 14) = 6.906992;
	A(24, 15) = 7.962222;
	A(24, 16) = 8.315734;
	A(24, 17) = 8.104349;
	A(24, 18) = 6.945464;
	A(24, 19) = 6.871617;
	A(24, 20) = 8.486403;
	A(24, 21) = 6.103824;
	A(24, 22) = 9.77928;
	A(24, 23) = 8.144623;
	A(24, 24) = 42.32337;
	A(24, 25) = 7.153761;
	A(24, 26) = 7.83137;
	A(24, 27) = 7.94244;
	A(24, 28) = 7.0099;
	A(24, 29) = 7.585384;
	A(24, 30) = 6.553133;
	A(24, 31) = 8.329912;
	A(24, 32) = 7.453521;
	A(25, 1) = 8.421722;
	A(25, 2) = 7.521159;
	A(25, 3) = 5.718482;
	A(25, 4) = 6.917626;
	A(25, 5) = 9.058475;
	A(25, 6) = 7.553156;
	A(25, 7) = 9.262994;
	A(25, 8) = 6.249789;
	A(25, 9) = 7.555898;
	A(25, 10) = 6.739049;
	A(25, 11) = 7.316482;
	A(25, 12) = 7.022617;
	A(25, 13) = 7.320691;
	A(25, 14) = 6.826281;
	A(25, 15) = 7.913155;
	A(25, 16) = 9.747427;
	A(25, 17) = 7.359839;
	A(25, 18) = 7.284826;
	A(25, 19) = 7.357392;
	A(25, 20) = 8.834362;
	A(25, 21) = 6.095469;
	A(25, 22) = 9.771961;
	A(25, 23) = 8.504797;
	A(25, 24) = 7.153761;
	A(25, 25) = 40.94436;
	A(25, 26) = 8.150279;
	A(25, 27) = 7.591934;
	A(25, 28) = 7.281586;
	A(25, 29) = 7.511089;
	A(25, 30) = 6.754559;
	A(25, 31) = 7.781699;
	A(25, 32) = 7.098694;
	A(26, 1) = 9.0935;
	A(26, 2) = 8.532229;
	A(26, 3) = 5.395215;
	A(26, 4) = 7.476358;
	A(26, 5) = 10.3263;
	A(26, 6) = 8.095305;
	A(26, 7) = 9.333155;
	A(26, 8) = 6.456123;
	A(26, 9) = 7.915104;
	A(26, 10) = 6.105246;
	A(26, 11) = 7.170356;
	A(26, 12) = 6.983741;
	A(26, 13) = 8.708432;
	A(26, 14) = 8.074467;
	A(26, 15) = 7.705297;
	A(26, 16) = 9.661396;
	A(26, 17) = 8.698833;
	A(26, 18) = 8.090492;
	A(26, 19) = 7.528459;
	A(26, 20) = 9.238183;
	A(26, 21) = 6.83822;
	A(26, 22) = 9.633077;
	A(26, 23) = 8.80019;
	A(26, 24) = 7.83137;
	A(26, 25) = 8.150279;
	A(26, 26) = 42.18546;
	A(26, 27) = 8.857089;
	A(26, 28) = 7.619456;
	A(26, 29) = 8.058824;
	A(26, 30) = 7.75767;
	A(26, 31) = 7.386223;
	A(26, 32) = 7.773975;
	A(27, 1) = 9.107638;
	A(27, 2) = 9.692117;
	A(27, 3) = 6.688468;
	A(27, 4) = 9.164117;
	A(27, 5) = 10.59962;
	A(27, 6) = 7.494553;
	A(27, 7) = 8.90611;
	A(27, 8) = 6.674888;
	A(27, 9) = 8.719343;
	A(27, 10) = 8.051716;
	A(27, 11) = 7.771206;
	A(27, 12) = 8.927778;
	A(27, 13) = 9.071458;
	A(27, 14) = 8.183154;
	A(27, 15) = 7.734393;
	A(27, 16) = 10.50385;
	A(27, 17) = 8.458376;
	A(27, 18) = 8.94317;
	A(27, 19) = 7.624521;
	A(27, 20) = 9.052507;
	A(27, 21) = 7.735146;
	A(27, 22) = 10.41915;
	A(27, 23) = 9.413723;
	A(27, 24) = 7.94244;
	A(27, 25) = 7.591934;
	A(27, 26) = 8.857089;
	A(27, 27) = 43.7624;
	A(27, 28) = 8.759824;
	A(27, 29) = 9.119857;
	A(27, 30) = 8.281866;
	A(27, 31) = 8.421967;
	A(27, 32) = 9.493143;
	A(28, 1) = 9.442419;
	A(28, 2) = 7.942288;
	A(28, 3) = 6.211776;
	A(28, 4) = 7.526169;
	A(28, 5) = 9.52149;
	A(28, 6) = 7.376231;
	A(28, 7) = 10.33834;
	A(28, 8) = 6.155236;
	A(28, 9) = 8.830105;
	A(28, 10) = 7.516197;
	A(28, 11) = 7.801772;
	A(28, 12) = 8.285803;
	A(28, 13) = 8.533298;
	A(28, 14) = 7.966928;
	A(28, 15) = 7.216666;
	A(28, 16) = 10.15614;
	A(28, 17) = 8.151227;
	A(28, 18) = 8.657142;
	A(28, 19) = 6.722559;
	A(28, 20) = 8.407161;
	A(28, 21) = 6.095053;
	A(28, 22) = 9.951712;
	A(28, 23) = 9.23968;
	A(28, 24) = 7.0099;
	A(28, 25) = 7.281586;
	A(28, 26) = 7.619456;
	A(28, 27) = 8.759824;
	A(28, 28) = 42.07356;
	A(28, 29) = 8.375107;
	A(28, 30) = 7.169741;
	A(28, 31) = 8.428096;
	A(28, 32) = 7.784962;
	A(29, 1) = 9.134119;
	A(29, 2) = 8.72569;
	A(29, 3) = 6.30242;
	A(29, 4) = 9.105786;
	A(29, 5) = 9.969017;
	A(29, 6) = 8.644343;
	A(29, 7) = 8.989779;
	A(29, 8) = 7.478709;
	A(29, 9) = 9.526395;
	A(29, 10) = 8.150392;
	A(29, 11) = 8.555453;
	A(29, 12) = 7.814824;
	A(29, 13) = 9.012588;
	A(29, 14) = 8.061718;
	A(29, 15) = 8.057176;
	A(29, 16) = 10.52856;
	A(29, 17) = 7.824457;
	A(29, 18) = 8.664519;
	A(29, 19) = 8.055284;
	A(29, 20) = 9.306524;
	A(29, 21) = 6.755796;
	A(29, 22) = 9.582796;
	A(29, 23) = 9.801824;
	A(29, 24) = 7.585384;
	A(29, 25) = 7.511089;
	A(29, 26) = 8.058824;
	A(29, 27) = 9.119857;
	A(29, 28) = 8.375107;
	A(29, 29) = 42.57956;
	A(29, 30) = 6.742315;
	A(29, 31) = 8.466095;
	A(29, 32) = 7.870807;
	A(30, 1) = 8.025165;
	A(30, 2) = 6.940404;
	A(30, 3) = 5.29562;
	A(30, 4) = 7.531226;
	A(30, 5) = 8.811292;
	A(30, 6) = 6.765675;
	A(30, 7) = 8.447423;
	A(30, 8) = 5.490128;
	A(30, 9) = 5.566744;
	A(30, 10) = 5.744286;
	A(30, 11) = 6.612849;
	A(30, 12) = 6.1109;
	A(30, 13) = 6.767958;
	A(30, 14) = 6.672546;
	A(30, 15) = 5.870655;
	A(30, 16) = 8.399444;
	A(30, 17) = 6.353978;
	A(30, 18) = 7.552076;
	A(30, 19) = 5.277095;
	A(30, 20) = 8.352698;
	A(30, 21) = 6.11867;
	A(30, 22) = 8.793287;
	A(30, 23) = 7.998706;
	A(30, 24) = 6.553133;
	A(30, 25) = 6.754559;
	A(30, 26) = 7.75767;
	A(30, 27) = 8.281866;
	A(30, 28) = 7.169741;
	A(30, 29) = 6.742315;
	A(30, 30) = 41.11817;
	A(30, 31) = 6.762022;
	A(30, 32) = 7.684819;
	A(31, 1) = 9.197884;
	A(31, 2) = 8.204024;
	A(31, 3) = 5.810936;
	A(31, 4) = 8.335414;
	A(31, 5) = 9.604378;
	A(31, 6) = 8.372642;
	A(31, 7) = 9.762229;
	A(31, 8) = 7.22629;
	A(31, 9) = 8.854994;
	A(31, 10) = 7.489577;
	A(31, 11) = 7.937045;
	A(31, 12) = 8.875031;
	A(31, 13) = 8.550621;
	A(31, 14) = 7.085254;
	A(31, 15) = 7.890541;
	A(31, 16) = 9.488288;
	A(31, 17) = 8.320092;
	A(31, 18) = 7.330491;
	A(31, 19) = 7.429545;
	A(31, 20) = 9.509018;
	A(31, 21) = 6.287632;
	A(31, 22) = 10.35086;
	A(31, 23) = 9.271368;
	A(31, 24) = 8.329912;
	A(31, 25) = 7.781699;
	A(31, 26) = 7.386223;
	A(31, 27) = 8.421967;
	A(31, 28) = 8.428096;
	A(31, 29) = 8.466095;
	A(31, 30) = 6.762022;
	A(31, 31) = 43.54748;
	A(31, 32) = 8.217811;
	A(32, 1) = 8.631121;
	A(32, 2) = 8.361619;
	A(32, 3) = 6.367104;
	A(32, 4) = 9.13942;
	A(32, 5) = 8.815245;
	A(32, 6) = 7.215221;
	A(32, 7) = 8.705503;
	A(32, 8) = 6.534867;
	A(32, 9) = 7.674169;
	A(32, 10) = 7.6337;
	A(32, 11) = 7.798531;
	A(32, 12) = 8.762041;
	A(32, 13) = 8.117458;
	A(32, 14) = 7.930122;
	A(32, 15) = 7.414451;
	A(32, 16) = 9.458236;
	A(32, 17) = 8.007373;
	A(32, 18) = 7.377906;
	A(32, 19) = 6.107671;
	A(32, 20) = 8.387606;
	A(32, 21) = 7.523602;
	A(32, 22) = 10.46545;
	A(32, 23) = 8.502528;
	A(32, 24) = 7.453521;
	A(32, 25) = 7.098694;
	A(32, 26) = 7.773975;
	A(32, 27) = 9.493143;
	A(32, 28) = 7.784962;
	A(32, 29) = 7.870807;
	A(32, 30) = 7.684819;
	A(32, 31) = 8.217811;
	A(32, 32) = 42.68885;


	return A;
}
