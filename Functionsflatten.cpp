#include <iostream>
#include <string>
#include <fstream>
#include <bits/stdc++.h>
#include <math.h>

#include "LUflatten.cpp"
#include "LUPDflatten.cpp"
#include "Rutishauserflatten.cpp"
#include "InverseIterationflatten.cpp"

using namespace std;

  /*****************************************************
  * INPUTS:                                            *
  * EPS : precision (Double)                           *
  * D1  : precision d1 (Double)                        *
  * M   : maximum number of iterations (integer)       *
  * N   : order of matrix A (integer)                  *
  * A   : input matrix to study (of MAT type)          *
  * -------------------------------------------------- *
  * OUTPUTS:                                           *
  * R   : table of eigenvalues (of VEC type)           *
  * VX  : table of eigenvectors (of MAT type)          *
  *****************************************************/

void CAL(double *A, double *R, double *VX, int n, double d1, int m, double tol, double *lambda)                     //InverseIteration.cpp
{   
	double *X = new double [n];
    double *p  = new double [n];
    
    double *A1 = new double [n * n];    
    double *lower = new double [n * n];
    double *upper = new double [n * n];
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A1[i * n + j] = A[i * n + j];
        }        
    }

    clock_t t;
    t=clock();    
    RUTIS(d1, A1, R, upper, lower, p, n, m);
    t = clock() - t;
    printf ("It took me for rutishaser %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
   
    // restore A1 after VAMR
   
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A1[i * n + j] = A[i * n + j];
        }        
    }
    
    clock_t t1;
    t1=clock();
    int l = 0;
    while (l < n)
    {        
        //if (R[l] < 1)        
        //{
            IIM(A1, R[l], l, R, X, upper, lower, p, n, tol, m, lambda);                                //InverseIteration.cpp
        //}                                                  
  
        // restore A1 after IIM
        for (int i = 0; i < n; i++)
        {
            for (int k = 0; k < n; k++)
            {
                A1[i * n + k] = A[i * n + k];
            }
            VX[i * n + l] = X[i];        
        }
        for (int i = 0; i < n; i++)
        {
            for (int k = 0; k < n; k++)
            {
                upper[i * n + k] = {0};
                lower[i * n + k] = {0};
            }                    
        }
        l++;     
    }
    t1 = clock() - t1;
    printf ("It took me for inverse power %d clicks (%f seconds).\n",t1,((float)t1)/CLOCKS_PER_SEC);

    delete [] X;
    delete [] A1;
    delete [] upper;
    delete [] lower;
    delete [] p;
}


