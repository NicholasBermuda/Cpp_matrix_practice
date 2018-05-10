# Cpp_matrix_practice
Practicing C++ with some linear algebra. The goal is to write Matrix and Vector classes that have similar syntax to MATLAB with some familiar MATLAB functionality (i.e. linear solvers). This is an extension to practical work in a C++ course.

All code is mine unless otherwise noted!

# Files included
* Exception.h/cpp - a useful class for throwing exceptions. Written entirely by Joe Pitt-Francis
* Vector.h/cpp - a class of vectors with some associated linear algebra functionality. Basic functionality by Joe Pitt-Francis, but I've significantly extended it and reformatted.
* Matrix.h/cpp - a class of matrices with some associated linear algebra functionality, and significant crossover with Vector class.
* test_matrix.cpp - some basic unit tests for the Matrix class. Executed whenever Matrix.o is recompiled - handy!
* use_matrix.cpp - a brief exploration of the rate of growth of the LU, QR, CG, and GMRES linear solvers in my Matrix class.
* matrix_writer.py - a utility program to convert a MATLAB csv file of [a matrix A, RHS b, and solution x of the problem Ax = b] into C++ code compatible with my Matrix class -- used for writing a unit test. Requirements: appropriately named csv files
* Makefile - to compile the above appropriately. USAGE: simply `make all`
