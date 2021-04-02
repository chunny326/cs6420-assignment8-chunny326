#ifndef _ADV_ALG_MATRICES_H_
#define _ADV_ALG_MATRICES_H_

// Multiplies matrix A by matrix B and returns the result based on the Square-Matrix-Multiply algorithm. 
// Both A and B matrices are of square n×n size and n is a power of 2.
double* iterative_mult(int n, double *A, double *B);

// Multiplies matrix A by matrix B and returns the result based on the Square-Matrix-Multiply-Recursive 
// algorithm. Both A and B matrices are of square n×n size and n is a power of 2.
double* recursive_mult(int n, double* A, double* B);

// Multiplies matrix A by matrix B and returns the result using the Strassen's method. 
// Both A and B matrices are of square n×n size and n is a power of 2.
double* strassen_mult(int n, double* A, double* B);

// takes a column-major array A and prints in its a two-dimensional matrix format
void print_in_2d(int n, double* A);

// add two matrices A and B, return sum matrix C
double* add_matrices(int n, double* A, double* B);

// check if two matrices A and B are equivalent
bool check_equal(int n, double* A, double* B);

#endif