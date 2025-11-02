#include <iostream>
#include <string>
#include <fstream>
#include <bits/stdc++.h>
#include <math.h>
#include <time.h>

using namespace std;

/********************************************************
  * Calculate the eigenvalues of a real square matrix by  *
  * Rutishauser's Method.                                 *
  * ----------------------------------------------------- *
  * INPUTS:                                               *
  *        eps: absolute precision (double)               *
  *        dta: relative precision (double)               *
  *         m : maximum number of iterations (integer)    *
  *         n : size of matrix A (integer)                *
  *         A : input real square matrix (n x n)          *
  * OUTPUTS:                                              *
  *         it: flag, =-1 if convergence is not obtained  *
  *                   =1 if convergence is ok.            *
  *         R : contains in output the n eigenvalues of A *
  *                                                       *         
  ********************************************************/
void RUTIS(double dta, double *A, double *R, double *upper, double *lower, double *p, int n, int m) 
{
    int l;
    double phi, s, t0;
          
    t0 = 0.0;
    l = 0;
    int it = -1;
    while (l < m && it != 1)
    {
        for (int i = 0; i < n; i++)
        {
            R[i] = A[i * n + i];
        }
        //clock_t t;
        //t=clock();
        
        luDecomposition(A, upper, lower, p, n);
        //t = clock() - t;
        //printf ("It took me for LU %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
        
        if (it == 0)
        {
            for (int i = 0; i < n; i++)
            {
                A[i * n + i] = A[i * n + i] + 1.0;
            }
            t0 += 1.0;
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    s = 0.0;
                    for (int k = 0; k < n; k++)
                    {
                        s = s + upper[i * n + k] * lower[k * n + j];   //?
                    }
                    A[i * n + j] = s;
                }
            }
		  
	        phi = 0.0;
            for (int i = 0; i < n; i++)
            {
                s = fabs(A[i * n + i] - R[i]);
                if (s > phi)
                {
                    phi = s;
                }  
		    }
	        if (phi < dta)
            {
                for (int i = 0; i < n; i++)
                {
                    R[i] = A[i * n + i] - t0;
                }
                it = 1;  
            }
	        else
            {
                l++;
                it = -1;
		    }
        }
    } //while
}

  
