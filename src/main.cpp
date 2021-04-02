#include <iostream>
#include "alg_matrices.h"
#include <vector>
using namespace std;

int main() 
{
  // -----------------------------------------------------------------------
  // 2x2 matrices
  int n = 2;
  // solution should be --> r1:  r2: 
  double* a1 = new double[n*n] { 3.0, 4.0, 5.0, 0.0 };
  double* b1 = new double[n*n] { 1.0, 2.0, 7.0, 1.0 };

  cout << "2 x 2 matrix test..." << endl;
  auto r1 = iterative_mult(n, a1, b1);
  cout << "Iterative: " << endl;
  print_in_2d(n, r1);
  auto s1 = recursive_mult(n, a1, b1);
  cout << "\nRecursive: " << endl;
  print_in_2d(n, s1);
  double* t1 = strassen_mult(n, a1, b1);
  cout << "\nStrassen: " << endl;
  print_in_2d(n, t1);

  delete[] a1; 
  delete[] b1;
  delete[] r1;
  delete[] s1;
  delete[] t1;

  // -----------------------------------------------------------------------
  // 4x4 matrices
  n = 4;
  // solution should be --> r1:  r2:
  double* a2 = new double[n*n] { 3.0, 4.0, 5.0, 0.0, 3.0, 4.0, 5.0, 0.0, 
                                 3.0, 4.0, 5.0, 0.0, 3.0, 4.0, 5.0, 0.0 };
  double* b2 = new double[n*n] { 1.0, 2.0, 7.0, 1.0, 1.0, 2.0, 7.0, 1.0,
                                 1.0, 2.0, 7.0, 1.0, 1.0, 2.0, 7.0, 1.0 };

  cout << "\n4 x 4 matrix test..." << endl;
  auto r2 = iterative_mult(n, a2, b2);
  cout << "Iterative: " << endl;
  print_in_2d(n, r2);
  double* s2 = recursive_mult(n, a2, b2);
  cout << "\nRecursive: " << endl;
  print_in_2d(n, s2);
  double* t2 = strassen_mult(n, a2, b2);
  cout << "\nStrassen: " << endl;
  print_in_2d(n, t2);

  delete[] a2; 
  delete[] b2;
  delete[] r2;
  delete[] s2;
  delete[] t2;

  return 0;
}
