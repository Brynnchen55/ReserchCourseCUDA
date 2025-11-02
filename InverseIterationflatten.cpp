#include <iostream>
#include <string>
#include <fstream>
#include <bits/stdc++.h>
#include <math.h>

using namespace std;

  /***********************************************************
  * Procedure IIM calculates a real eigenvalue and the asso- *
  * ciated eigenvector of a real square matrix the inverse   *
  * iteration method.                                        *
  * -------------------------------------------------------- *
  * INPUTS:                                                  *
  *         eps : absolute precision (double)                *
  *         dta : relative precision (double)                *
  *          m  : maximum number of iterations (integer)     *
  *          n  : size of matrix A                           *
  *          A  : input real square matrix (n x n)           *
  * OUTPUTS:                                                 *
  *        Gamma: starting value for the eigenvalue as input *
  *               approximation of the eigenvalue with preci-*
  *               sion dta in output.                        *
  *          X1 : contains in output the associated eigen-   *
  *               vector.                                    *
  *                                                          *
  ***********************************************************/
void IIM(double *A, double gamma, int ll, double *R, double *X1, double *upper, double *lower, double *p, int n, double tol, int m, double *sigma_all)
{
    int l, l0, number;
    double  p0, phi, s, t0, q0, q00, z0, costheta, qq, temp, zz, relres, sigma;

    double *W = new double [n];
    double *b = new double [n];
    double *X0 = new double [n];
    double *q = new double [n];
    double *q1 = new double [n];
    double *q2 = new double [n];
    double *y = new double [n];
    double *z = new double [n];

    double *r = new double [n * n];
    double *u = new double [n * n];
    double *r1 = new double [n * n];
    double *u1 = new double [n * n];
    double *A3 = new double [n * n];
    double *pp = new double [n * n];
            
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A3[i * n + j] = A[i * n + j]; 
        } 
    }         
    
    for (int i = 0; i < n; i++)
    {
        A3[i * n + i] = A3[i * n + i] - gamma;
    }
    LUP_Descomposition(A3, lower, upper, p, n);
    

    for (int i = 0; i < n; i++)
    {
        number = p[i];
        pp[i * n + number] = 1;
    }

    u = LUP_solve_inverse(upper, p, n);
    r = LUP_solve_inverse(lower, p, n);
   
    transpose(upper, n, n);
    transpose(lower, n, n);

    u1 = LUP_solve_inverse(upper, p, n);
    r1 = LUP_solve_inverse(lower, p, n);
   
    //Norm of X0
    q0 = 0;
    for (int i = 0; i < n; i++)
    {
        X0[i] = 1.0;
        q0 = q0 + X0[i] * X0[i];
    } 
    for (int i = 0; i < n; i++)
    {
        q[i] = X0[i] / sqrt (q0);
        q2[i] = q[i];       
    } 

    relres = tol + 1;
    l = 0;
    sigma = 0;
    while (relres > tol && l < m)
    {
        l++;

        for (int i = 0; i < n; i++)
        {
            b[i] = 0;            
        }
        for (int i = 0; i < n ; i++)
        {
            for (int j = 0; j < n; j++)
            {
                b[i] = b[i] + pp[i * n + j] * q[j];          
            }
        } 

        for (int i = 0; i < n; i++)
        {
            y[i] = 0;            
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                y[i] = y[i] + r[i * n + j] * b[j]; 
            } 
        } 
        for (int i = 0; i < n; i++)
        {
            z[i] = 0;            
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                z[i] = z[i] + u[i * n + j] * y[j];
            } 
        }

        z0 = 0;
        for (int i = 0; i < n; i++)
        {
            z0 = z0 + z[i] * z[i];
        }
        
        for (int i = 0; i < n; i++)
        {
            q[i] = z[i] / sqrt (z0);
        }

        for (int i = 0; i < n; i++)
        {
            z[i] = 0;            
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                z[i] = z[i] + A[i * n + j] * q[j]; 
            }
        } 

        sigma = 0;
        for (int i = 0; i < n; i++)
        {
            sigma = sigma + q[i] * z[i];
        } 
        
        
        for (int i = 0; i < n ; i++)
        {
            b[i] = q2[i];
        }
        
        for (int i = 0; i < n; i++)
        {
            y[i] = 0;            
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                y[i] = y[i] + u1[i * n + j] * b[j];
            } 
        }         

        for (int i = 0; i < n; i++)
        {
            W[i] = 0;            
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                W[i] = W[i] + r1[i * n + j] * y[j];
            } 
        } 

        for (int i = 0; i < n; i++)
        {
            q2[i] = 0;            
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                q2[i] = q2[i] + W[j] * pp[j * n + i];
            }              
        }

        q00 = 0;
        for (int i = 0; i < n; i++)
        {
            q00 = q00 + q2[i] * q2[i];
        } 
        for (int i = 0; i < n; i++)
        {
            q2[i] = q2[i] / sqrt (q00); 
        }
        qq = 0;
        for (int i = 0; i < n; i++)
        {
            qq = qq + q2[i] * q[i];            
        }
        costheta = fabs(qq); 

        
        if (costheta >= 5e-2)
        {
            zz = 0;
            for (int i = 0; i < n; i++)
            {
                zz = zz + ((z[i] - (sigma * q[i])) * (z[i] - (sigma * q[i])));
            }
            temp = sqrt(zz) / costheta; 
            relres = temp; 
        } 
        for (int i = 0; i < n; i++)
        {
            X1[i] = q[i];
        }
    }
    sigma_all[ll] = sigma;
    
    delete [] X0;
    delete [] W;
    delete [] b;
    delete [] q;
    delete [] q2;
    delete [] y;
    delete [] z;
    
    delete [] r;
    delete [] u;
    delete [] r1;
    delete [] u1;
    delete [] A3;
}    	  						    


