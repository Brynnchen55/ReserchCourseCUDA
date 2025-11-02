#include <iostream>
#include <string>
#include <fstream>
#include <bits/stdc++.h>
#include <math.h>

using namespace std;

void luDecomposition(double *A, double *upper, double *lower, double *p, int n)
{
    double mm, l;

    double *q = new double [n];
    double *b = new double [n];
    double *A2 = new double [n * n];;
                
    // Decomposing matrix into Upper and Lower
    // triangular matrix
    /* Initialize L */
    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            if (i == j) 
            {
                lower[i * n + j] = 1.0;
            }
            else 
            {
                lower[i * n + j] = 0.0;
            }
        }
    }

    /* Initialize U */
    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            if (i > j) 
            {
                upper[i * n + j] = 0.0;
            }
        }
    }
    
    /* LU decomposition */
    for(int k = 0; k < n; ++k) 
    {
        upper[k * n + k] = A[k * n + k];
        for(int i = k + 1; i < n; ++i) 
        {
            lower[i * n + k] = A[i * n + k] / upper[k * n + k];
            upper[k * n + i] = A[k * n + i];
        }
        for(int i = k + 1; i < n; ++i) 
        {
            for(int j = k + 1; j < n; ++j) 
            {
                A[i * n + j] = A[i * n + j] - lower[i * n + k] * upper[k * n + j];
            }
        }
    }
        
    delete [] A2;
    delete [] q;
    delete [] b;

    return;
}

