#include <iostream>
#include <iomanip>
#include <random>
#include "alg_stopwatch.h"
#include "alg_matrices.h"

using namespace std;

int main() 
{
  // use StopWatch class to record running times of algorithms
  StopWatch sw;
  sw.reset();

  // get different n x n matrix sizes up to 256, powers of 2
  vector<int> n_arr { 2, 4, 8, 16, 32, 64, 128, 256 };
  int n { 2 };
  int size { 4 };

  // arrays to store times for running each algorithm
  double iter_times[n_arr.size()];
  double rec_times[n_arr.size()];
  double strass_times[n_arr.size()];

  // run for size of n times for 2^1, 2^2, ... , 2^8
  // power of n is doubled at the end of each iteration
  for (int j = 0; j < n_arr.size(); j++)
  {
    n = n_arr[j];
    size = n * n;
    double* A = new double[size];
    double* B = new double[size];

    // generate random doubles from 0.0 - 10.0
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine gen(seed);
    uniform_real_distribution<double> num(0.0, 10.0);

    // fill n x n matrices with random doubles from 0.0 - 10.0, rounded to whole numbers
    for (int i = 0; i < size; i++)
    {
      A[i] = round(num(gen));
      B[i] = round(num(gen));
    }
    
    // compare running times for 3 matrix matrix multiplication algorithms
    sw.reset();
    double* C1 = iterative_mult(n, A, B);
    iter_times[j] = sw.elapsed_time();
    sw.reset();
    double* C2 = recursive_mult(n, A, B);
    rec_times[j] = sw.elapsed_time();
    sw.reset();
    double* C3 = strassen_mult(n, A, B);
    strass_times[j] = sw.elapsed_time();
    sw.reset();

    delete[] A;
    delete[] B;
    delete[] C1;
    delete[] C2;
    delete[] C3;
  }
   
  // display algorithm timing results for each n x n matrix multiplication in a table-like format
  cout << "----------------------------------------------Matrix-Matrix-Multiply Algorithms Timing Results----------------------------------------------" << endl;
  cout << "N:            ";
  for (int i = 0; i < n_arr.size(); i++)
  {
    cout << setw(12) << n_arr[i] << "    ";
  }
  
  cout << "\nIterative:    ";
  for (int i = 0; i < n_arr.size(); i++)
  {
    cout << setw(12) << iter_times[i] << "    ";
  }

  cout << "\nRecursive:    ";
  for (int i = 0; i < n_arr.size(); i++)
  {
    cout << setw(12) << rec_times[i] << "    ";
  }

  cout << "\nStrassen:     ";
  for (int i = 0; i < n_arr.size(); i++)
  {
    cout << setw(12) << strass_times[i] << "    ";
  }
  cout << "\n";

  return 0;
}