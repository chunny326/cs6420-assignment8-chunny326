#define CATCH_CONFIG_MAIN
#include "catch2.h"
#include "alg_matrices.h"
#include <iostream>

using namespace std;

// 2x2 matrix
TEST_CASE("2x2 matrix", "2x2") {
	// 2x2 matrix3
	int n = 2;

	// multiply these two matrices
	double* a1 = new double[n*n] { 3.0, 4.0, 5.0, 0.0 };
	double* b1 = new double[n*n] { 1.0, 2.0, 7.0, 1.0 };

	std::cout << "2 x 2 matrix test..." << std::endl;

	// print A and B matrices to be multiplied
	std::cout << "Multiply A:" << std::endl;
	print_in_2d(n, a1);
	std::cout << "\nand B:" << std::endl;
	print_in_2d(n, b1);

	std::cout << "\nResult for each algorithm: " << std::endl;
	auto r1 = iterative_mult(n, a1, b1);
	std::cout << "Iterative: " << std::endl;
	print_in_2d(n, r1);
	auto s1 = recursive_mult(n, a1, b1);
	std::cout << "\nRecursive: " << std::endl;
	print_in_2d(n, s1);
	double* t1 = strassen_mult(n, a1, b1);
	std::cout << "\nStrassen: " << std::endl;
	print_in_2d(n, t1);

	CHECK(check_equal(n, r1, s1));
	CHECK(check_equal(n, s1, t1));
	
	delete[] a1; 
	delete[] b1;
	delete[] r1;
	delete[] s1;
	delete[] t1;
}

// 4x4 matrix
TEST_CASE("4x4 matrix", "4x4") {
	// -----------------------------------------------------------------------
	// 4x4 matrices
	int n = 4;
	
	// multiply these two matrices
	double* a1 = new double[n*n] { 1.0, 4.0, 5.0, 20.0, 3.0, 4.0, 5.0, 7.0, 
								   3.0, 14.0, 5.0, 2.0, 3.0, 4.0, 5.0, 1.0 };
	double* b1 = new double[n*n] { 1.0, 4.0, 7.0, 8.0, 1.0, 2.0, 5.0, 1.0,
								   12.0, 2.0, 7.0, 1.0, 1.0, 2.0, 17.0, 1.0 };

	std::cout << "\n4 x 4 matrix test..." << std::endl;

	// print A and B matrices to be multiplied
	std::cout << "Multiply A:" << std::endl;
	print_in_2d(n, a1);
	std::cout << "\nand B:" << std::endl;
	print_in_2d(n, b1);

	std::cout << "\nResult for each algorithm: " << std::endl;
	auto r1 = iterative_mult(n, a1, b1);
	std::cout << "Iterative: " << std::endl;
	print_in_2d(n, r1);
	auto s1 = recursive_mult(n, a1, b1);
	std::cout << "\nRecursive: " << std::endl;
	print_in_2d(n, s1);
	double* t1 = strassen_mult(n, a1, b1);
	std::cout << "\nStrassen: " << std::endl;
	print_in_2d(n, t1);

	CHECK(check_equal(n, r1, s1));
	CHECK(check_equal(n, s1, t1));

	delete[] a1; 
	delete[] b1;
	delete[] r1;
	delete[] s1;
	delete[] t1;
}

// 8x8 matrix
TEST_CASE("8x8 matrix", "8x8") {
	// -----------------------------------------------------------------------
	// 8x8 matrices
	int n = 8;
	
	// multiply these two matrices
	double* a1 = new double[n*n] { 1.0, 7.0, 8.0, 11.0, 23.0, 4.0, 2.0, 5.0, 
								   9.0, 5.0, 4.0, 11.0, 22.0, 8.0, 5.0, 6.0,
								   11.0, 14.0, 3.0, 1.0, 7.0, 4.0, 1.0, 27.0,
								   8.0, 4.0, 1.0, 8.0, 17.0, 14.0, 11.0, 9.0,
								   8.0, 8.0, 7.0, 1.0, 3.0, 14.0, 12.0, 2.0, 
								   6.0, 6.0, 8.0, 1.0, 12.0, 5.0, 1.0, 1.0,
								   7.0, 3.0, 21.0, 11.0, 17.0, 14.0, 11.0, 2.0,
								   18.0, 3.0, 4.0, 9.0, 8.0, 7.0, 1.0, 4.0 };
	double* b1 = new double[n*n] { 8.0, 8.0, 7.0, 1.0, 3.0, 14.0, 12.0, 2.0, 
								   6.0, 6.0, 8.0, 1.0, 12.0, 5.0, 1.0, 1.0,
								   7.0, 3.0, 21.0, 11.0, 17.0, 14.0, 11.0, 2.0,
								   18.0, 3.0, 4.0, 9.0, 8.0, 7.0, 1.0, 4.0,
								   1.0, 7.0, 8.0, 11.0, 23.0, 4.0, 2.0, 5.0, 
								   9.0, 5.0, 4.0, 11.0, 22.0, 8.0, 5.0, 6.0,
								   11.0, 14.0, 3.0, 1.0, 7.0, 4.0, 1.0, 27.0,
								   8.0, 4.0, 1.0, 8.0, 17.0, 14.0, 11.0, 9.0 };

	std::cout << "\n8 x 8 matrix test..." << std::endl;

	// print A and B matrices to be multiplied
	std::cout << "Multiply A:" << std::endl;
	print_in_2d(n, a1);
	std::cout << "\nand B:" << std::endl;
	print_in_2d(n, b1);

	std::cout << "\nResult for each algorithm: " << std::endl;
	auto r1 = iterative_mult(n, a1, b1);
	std::cout << "Iterative: " << std::endl;
	print_in_2d(n, r1);
	auto s1 = recursive_mult(n, a1, b1);
	std::cout << "\nRecursive: " << std::endl;
	print_in_2d(n, s1);
	double* t1 = strassen_mult(n, a1, b1);
	std::cout << "\nStrassen: " << std::endl;
	print_in_2d(n, t1);

	CHECK(check_equal(n, r1, s1));
	CHECK(check_equal(n, s1, t1));

	delete[] a1; 
	delete[] b1;
	delete[] r1;
	delete[] s1;
	delete[] t1;
}

// 16x16 matrix
TEST_CASE("16x16 matrix", "16x16") {
	// -----------------------------------------------------------------------
	// 16x16 matrices
	int n = 16;
	
	// multiply these two matrices
	double* a1 = new double[n*n] { 1.0, 7.0, 8.0, 11.0, 23.0, 4.0, 2.0, 5.0, 9.0, 5.0, 4.0, 11.0, 22.0, 8.0, 5.0, 6.0,
								   11.0, 14.0, 3.0, 1.0, 7.0, 4.0, 1.0, 27.0, 8.0, 4.0, 1.0, 8.0, 17.0, 14.0, 11.0, 9.0,
								   8.0, 8.0, 7.0, 1.0, 3.0, 14.0, 12.0, 2.0, 6.0, 6.0, 8.0, 1.0, 12.0, 5.0, 1.0, 1.0,
								   7.0, 3.0, 21.0, 11.0, 17.0, 14.0, 11.0, 2.0, 18.0, 3.0, 4.0, 9.0, 8.0, 7.0, 1.0, 4.0,
								   1.0, 7.0, 8.0, 11.0, 23.0, 4.0, 2.0, 5.0, 9.0, 5.0, 4.0, 11.0, 22.0, 8.0, 5.0, 6.0,
								   11.0, 14.0, 3.0, 1.0, 7.0, 4.0, 1.0, 27.0, 8.0, 4.0, 1.0, 8.0, 17.0, 14.0, 11.0, 9.0,
								   8.0, 8.0, 7.0, 1.0, 3.0, 14.0, 12.0, 2.0, 6.0, 6.0, 8.0, 1.0, 12.0, 5.0, 1.0, 1.0,
								   7.0, 3.0, 21.0, 11.0, 17.0, 14.0, 11.0, 2.0, 18.0, 3.0, 4.0, 9.0, 8.0, 7.0, 1.0, 4.0,
								   1.0, 7.0, 8.0, 11.0, 23.0, 4.0, 2.0, 5.0, 9.0, 5.0, 4.0, 11.0, 22.0, 8.0, 5.0, 6.0,
								   11.0, 14.0, 3.0, 1.0, 7.0, 4.0, 1.0, 27.0, 8.0, 4.0, 1.0, 8.0, 17.0, 14.0, 11.0, 9.0,
								   8.0, 8.0, 7.0, 1.0, 3.0, 14.0, 12.0, 2.0, 6.0, 6.0, 8.0, 1.0, 12.0, 5.0, 1.0, 1.0,
								   7.0, 3.0, 21.0, 11.0, 17.0, 14.0, 11.0, 2.0, 18.0, 3.0, 4.0, 9.0, 8.0, 7.0, 1.0, 4.0,
								   1.0, 7.0, 8.0, 11.0, 23.0, 4.0, 2.0, 5.0, 9.0, 5.0, 4.0, 11.0, 22.0, 8.0, 5.0, 6.0,
								   11.0, 14.0, 3.0, 1.0, 7.0, 4.0, 1.0, 27.0, 8.0, 4.0, 1.0, 8.0, 17.0, 14.0, 11.0, 9.0,
								   8.0, 8.0, 7.0, 1.0, 3.0, 14.0, 12.0, 2.0, 6.0, 6.0, 8.0, 1.0, 12.0, 5.0, 1.0, 1.0,
								   7.0, 3.0, 21.0, 11.0, 17.0, 14.0, 11.0, 2.0, 18.0, 3.0, 4.0, 9.0, 8.0, 7.0, 1.0, 4.0 };
	double* b1 = new double[n*n] { 8.0, 8.0, 7.0, 1.0, 3.0, 14.0, 12.0, 2.0, 6.0, 6.0, 8.0, 1.0, 12.0, 5.0, 1.0, 1.0,
								   7.0, 3.0, 21.0, 11.0, 17.0, 14.0, 11.0, 2.0, 18.0, 3.0, 4.0, 9.0, 8.0, 7.0, 1.0, 4.0,
								   1.0, 7.0, 8.0, 11.0, 23.0, 4.0, 2.0, 5.0, 9.0, 5.0, 4.0, 11.0, 22.0, 8.0, 5.0, 6.0,
								   11.0, 14.0, 3.0, 1.0, 7.0, 4.0, 1.0, 27.0, 8.0, 4.0, 1.0, 8.0, 17.0, 14.0, 11.0, 9.0,
								   8.0, 8.0, 7.0, 1.0, 3.0, 14.0, 12.0, 2.0, 6.0, 6.0, 8.0, 1.0, 12.0, 5.0, 1.0, 1.0,
								   7.0, 3.0, 21.0, 11.0, 17.0, 14.0, 11.0, 2.0, 18.0, 3.0, 4.0, 9.0, 8.0, 7.0, 1.0, 4.0,
								   1.0, 7.0, 8.0, 11.0, 23.0, 4.0, 2.0, 5.0, 9.0, 5.0, 4.0, 11.0, 22.0, 8.0, 5.0, 6.0,
								   11.0, 14.0, 3.0, 1.0, 7.0, 4.0, 1.0, 27.0, 8.0, 4.0, 1.0, 8.0, 17.0, 14.0, 11.0, 9.0,
								   8.0, 8.0, 7.0, 1.0, 3.0, 14.0, 12.0, 2.0, 6.0, 6.0, 8.0, 1.0, 12.0, 5.0, 1.0, 1.0,
								   7.0, 3.0, 21.0, 11.0, 17.0, 14.0, 11.0, 2.0, 18.0, 3.0, 4.0, 9.0, 8.0, 7.0, 1.0, 4.0,
								   1.0, 7.0, 8.0, 11.0, 23.0, 4.0, 2.0, 5.0, 9.0, 5.0, 4.0, 11.0, 22.0, 8.0, 5.0, 6.0,
								   11.0, 14.0, 3.0, 1.0, 7.0, 4.0, 1.0, 27.0, 8.0, 4.0, 1.0, 8.0, 17.0, 14.0, 11.0, 9.0,
								   8.0, 8.0, 7.0, 1.0, 3.0, 14.0, 12.0, 2.0, 6.0, 6.0, 8.0, 1.0, 12.0, 5.0, 1.0, 1.0,
								   7.0, 3.0, 21.0, 11.0, 17.0, 14.0, 11.0, 2.0, 18.0, 3.0, 4.0, 9.0, 8.0, 7.0, 1.0, 4.0,
								   1.0, 7.0, 8.0, 11.0, 23.0, 4.0, 2.0, 5.0, 9.0, 5.0, 4.0, 11.0, 22.0, 8.0, 5.0, 6.0,
								   11.0, 14.0, 3.0, 1.0, 7.0, 4.0, 1.0, 27.0, 8.0, 4.0, 1.0, 8.0, 17.0, 14.0, 11.0, 9.0 };

	std::cout << "\n16 x 16 matrix test..." << std::endl;

	// print A and B matrices to be multiplied
	std::cout << "Multiply A:" << std::endl;
	print_in_2d(n, a1);
	std::cout << "\nand B:" << std::endl;
	print_in_2d(n, b1);

	std::cout << "\nResult for each algorithm: " << std::endl;
	auto r1 = iterative_mult(n, a1, b1);
	std::cout << "Iterative: " << std::endl;
	print_in_2d(n, r1);
	auto s1 = recursive_mult(n, a1, b1);
	std::cout << "\nRecursive: " << std::endl;
	print_in_2d(n, s1);
	double* t1 = strassen_mult(n, a1, b1);
	std::cout << "\nStrassen: " << std::endl;
	print_in_2d(n, t1);

	CHECK(check_equal(n, r1, s1));
	CHECK(check_equal(n, s1, t1));
	
	delete[] a1; 
	delete[] b1;
	delete[] r1;
	delete[] s1;
	delete[] t1;
}