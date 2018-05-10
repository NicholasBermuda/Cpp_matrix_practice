#!/usr/bin/env python
import numpy as np
import os

def main():
	# remove old versions
	if os.path.isfile('Cpp_mat_code.txt'):
		os.remove('Cpp_mat_code.txt')
	if os.path.isfile('Cpp_vec_code.txt'):
		os.remove('Cpp_vec_code.txt')
	if os.path.isfile('Cpp_soln_code.txt'):
		os.remove('Cpp_soln_code.txt')

	# write the matrix code
	A = np.loadtxt('big_mat.csv',delimiter=',')
	with open('Cpp_mat_code.txt','a') as the_file:
		for i in range(A.shape[0]):
			for j in range(A.shape[1]):
					the_file.write('A(' + str(i+1) + ', ' + str(j+1) + ') = ' + str(A[i,j]) + ';\n')
	the_file.close()

	# write the vector code
	b = np.loadtxt('big_vec.csv',delimiter=',')
	with open('Cpp_vec_code.txt','a') as the_file:
		for i in range(b.size):
			the_file.write('b(' + str(i+1) + ') = ' + str(b[i]) + ';\n')
	the_file.close()

	# write the solution code
	b = np.loadtxt('big_soln.csv',delimiter=',')
	with open('Cpp_soln_code.txt','a') as the_file:
		for i in range(b.size):
			the_file.write('xtest(' + str(i+1) + ') = ' + str(b[i]) + ';\n')
	the_file.close()

	return



if __name__ == "__main__":
	main()