#include <iostream>
#include "alg_matrices.h"

// Multiplies matrix A by matrix B and returns the result based on the Square-Matrix-Multiply algorithm. 
// Both A and B matrices are of square n×n size and n is a power of 2.
double* iterative_mult(int n, double *A, double *B)
{
    double* C = new double[n * n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            C[j * n + i] = 0;
            for (int k = 0; k < n; k++)
            {
                C[j * n + i] += A[k * n + i] * B[j * n + k];
            }
        }
    }
    return C;
}

// return the sum of two column-major matrices that are n x n size
double* add_matrices(int n, double* A, double* B)
{
    double* C = new double[n * n];

    for (int i = 0; i < n; i++) 
    {
      for (int j = 0; j < n; j++)
      {
        C[j * n + i] = A[j * n + i] + B[j * n + i];
      }
   }
    return C;
}

// return the difference of two column-major matrices that are n x n size
double* sub_matrices(int n, double* A, double* B)
{
    double* C = new double[n * n];

    for (int i = 0; i < n; i++) 
    {
      for (int j = 0; j < n; j++)
      {
        C[j * n + i] = A[j * n + i] - B[j * n + i];
      }
   }
    return C;
}

// Multiplies matrix A by matrix B and returns the result based on the Square-Matrix-Multiply-Recursive 
// algorithm. Both A and B matrices are of square n×n size and n is a power of 2.
double* recursive_mult(int n, double* A, double* B)
{
    if (n == 1)
    {
        double* C = new double[1];
        C[0] = A[0] * B[0];
        return C;
    }

    double* C = new double[n * n];
    int k = n / 2;    
    
    double A11[k*k];
    double A12[k*k];
    double A21[k*k];
    double A22[k*k];
    double B11[k*k];
    double B12[k*k];
    double B21[k*k];
    double B22[k*k];

    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            A11[j * k + i] = A[j * n + i];
            A12[j * k + i] = A[(j + k) * n + i];
            A21[j * k + i] = A[j * n + (i + k)];
            A22[j * k + i] = A[(j + k) * n + (i + k)];
            B11[j * k + i] = B[j * n + i];
            B12[j * k + i] = B[(j + k) * n + i];
            B21[j * k + i] = B[j * n + (i + k)];
            B22[j * k + i] = B[(j + k) * n + (i + k)];
        }
    }
    
    auto C11 = add_matrices(k, recursive_mult(k, A11, B11), recursive_mult(k, A12, B21));   
    auto C12 = add_matrices(k, recursive_mult(k, A11, B12), recursive_mult(k, A12, B22));
    auto C21 = add_matrices(k, recursive_mult(k, A21, B11), recursive_mult(k, A22, B21));
    auto C22 = add_matrices(k, recursive_mult(k, A21, B12), recursive_mult(k, A22, B22));
    
    for (int i = 0; i < k; i++) 
    {
        for (int j = 0; j < k; j++)
        {
            // partition A, B, C as in equations 4.9
            C[j * n + i]             = C11[j * k + i];
            C[(j + k) * n + i]       = C12[j * k + i];
            C[j * n + (i + k)]       = C21[j * k + i];
            C[(j + k) * n + (i + k)] = C22[j * k + i];
        }
    }
    return C;
}

// Multiplies matrix A by matrix B and returns the result using the Strassen's method. 
// Both A and B matrices are of square n×n size and n is a power of 2.
double* strassen_mult(int n, double* A, double* B)
{
    if (n == 1)
    {
        double* C = new double[1];
        C[0] = A[0] * B[0];
        return C;
    }

    double* C = new double[n * n];
    int k = n / 2;        
    
    double A11[k*k];
    double A12[k*k];
    double A21[k*k];
    double A22[k*k];
    double B11[k*k];
    double B12[k*k];
    double B21[k*k];
    double B22[k*k];

    for (int i = 0; i < k; i++)
    {
        for (int j = 0; j < k; j++)
        {
            A11[j * k + i] = A[j * n + i];
            A12[j * k + i] = A[(j + k) * n + i];
            A21[j * k + i] = A[j * n + (i + k)];
            A22[j * k + i] = A[(j + k) * n + (i + k)];
            B11[j * k + i] = B[j * n + i];
            B12[j * k + i] = B[(j + k) * n + i];
            B21[j * k + i] = B[j * n + (i + k)];
            B22[j * k + i] = B[(j + k) * n + (i + k)];
        }
    }
    
    auto P1 = strassen_mult(k, A11, sub_matrices(k, B12, B22));   
    auto P2 = strassen_mult(k, add_matrices(k, A11, A12), B22);
    auto P3 = strassen_mult(k, add_matrices(k, A21, A22), B11);
    auto P4 = strassen_mult(k, A22, sub_matrices(k, B21, B11));
    auto P5 = strassen_mult(k, add_matrices(k, A11, A22), add_matrices(k, B11, B22));
    auto P6 = strassen_mult(k, sub_matrices(k, A12, A22), add_matrices(k, B21, B22));
    auto P7 = strassen_mult(k, sub_matrices(k, A11, A21), add_matrices(k, B11, B12));

    auto C11 = sub_matrices(k, add_matrices(k, add_matrices(k, P5, P4), P6), P2);
    auto C12 = add_matrices(k, P1, P2);
    auto C21 = add_matrices(k, P3, P4);
    auto C22 = sub_matrices(k, sub_matrices(k, add_matrices(k, P5, P1), P3), P7);
    
    for (int i = 0; i < k; i++) 
    {
        for (int j = 0; j < k; j++)
        {
            // partition A, B, C as in equations 4.9
            C[j * n + i]             = C11[j * k + i];
            C[(j + k) * n + i]       = C12[j * k + i];
            C[j * n + (i + k)]       = C21[j * k + i];
            C[(j + k) * n + (i + k)] = C22[j * k + i];
        }
    }
    return C;
}

// takes a column-major array A and prints in its a two-dimensional matrix format
void print_in_2d(int n, double* A)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << A[j*n + i] << " ";
        }
        std::cout << std::endl;
    }
}

bool check_equal(int n, double* A, double* B)
{
    for (int i = 0; i < n*n; i++)
    {
        if (A[i] != B[i])
        {
            return false;
        }
    }
    return true;
}