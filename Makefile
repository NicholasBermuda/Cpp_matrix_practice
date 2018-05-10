all: use_matrix test_matrix

Vector.o: Vector.cpp Vector.hpp
	g++ -Wall -O -c Vector.cpp

Matrix.o: Matrix.cpp Matrix.hpp
	g++ -Wall -O -c Matrix.cpp

use_matrix: Matrix.o Vector.o Exception.o use_matrix.cpp
	g++ -std=c++11 -Wall -O -o use_matrix Matrix.o Vector.o Exception.o use_matrix.cpp

test_matrix: Matrix.o Vector.o Exception.o test_matrix.cpp
	g++ -std=c++11 -Wall -O -o test_matrix Matrix.o Vector.o Exception.o test_matrix.cpp
	./test_matrix

clean:
	rm -f *.o *~ test_matrix use_matrix