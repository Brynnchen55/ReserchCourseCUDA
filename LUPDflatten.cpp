#include <iostream>
#include <string>
#include <fstream>
#include <bits/stdc++.h>
#include <math.h>

using namespace std;

//LUP Decomposition
void LUP_Descomposition(double *A, double *L, double *U, double *P, int n)
{
    int row = 0;
    for (int i = 0; i < n; i++)
    {
        P[i] = i;
    }
    for (int i = 0; i < n - 1; i++)
    {
        double p = 0.0;
        for (int j = i; j < n; j++)
            {
                if (abs(A[j * n + i]) > p)
                {
                    p = abs(A[j * n + i]);
                    row = j;
                }
            }
        if (0 == p)
        {
            cout<< "Can not compute the inverse of the matrix" <<endl;  
            return ;
        }
  
        //Switch P[i] and P[row]
        int tmp = P[i];
        P[i] = P[row];
        P[row] = tmp;

        double tmp2 = 0.0;
       for (int j = 0; j < n; j++)
        {
            //Switch A[i][j] and A[row][j]
            tmp2 = A[i * n + j];
            A[i * n + j] = A[row * n + j];
            A[row * n + j] = tmp2;
        }

        //LU Decomposition
        double u = A[i * n + i], l = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            l = A[j * n + i] / u;
            A[j * n + i] = l;
            for (int k = i + 1; k < n; k++)
            {
                A[j * n + k] = A[j * n + k] - A[i * n + k] * l;
            }
        }
 
    }
 
    //Construct L and U
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            if (i != j)
            {
                L[i * n + j] = A[i * n + j];
            }
            else
            {
                L[i * n + j] = 1;
            }
        }
        for (int k = i; k < n; k++)
        {
            U[i * n + k] = A[i * n + k];
        }
    }
   

}

//LUP solve
double * LUP_Solve(double *L, double *U, double *P, double *b, int n)
{
    double *x = new double[n]();
    double *y = new double[n]();

    //Forward
    for (int i = 0; i < n; i++)
    {
        int k = P[i];
        y[i] = b[k];
        for (int j = 0; j < i; j++)
        {
            y[i] = y[i] - L[i * n + j] * y[j];
        }
    }
    
    
    //Backward
    for (int i = n - 1; i >= 0; i--)
    {
        x[i] = y[i];
        for (int j = n - 1; j > i; j--)
        {
            x[i] = x[i] - U[i * n + j] * x[j];
        }
        x[i] /= U[i * n + i];
    }
    
    return x;
}

/*****************Matrix Tranpose (start)********************/

int getnext(int i, int m, int n)
{
    return (i % n) * m + i / n;
}

int getPre(int i, int m, int n)
{
    return (i % m) * n + i / m;
}

void movedata(double *mtx, int i, int m, int n)
{
    double temp = mtx[i]; 
    int cur = i;    
    int pre = getPre(cur, m, n);
    while (pre != i)
    {
        mtx[cur] = mtx[pre];
        cur = pre;
        pre = getPre(cur, m, n);
    }
    mtx[cur] = temp;
}
 
void transpose(double *mtx, int m, int n)
{
    for (int i = 0; i < m * n; ++i)
    {
        int next = getnext(i, m, n);
        while (next > i) 
            next = getnext(next, m, n);
            if(next == i)  
                movedata(mtx, i, m, n);
    }
}
/*****************Matrix Tranpose (end)********************/
 
//LUP Inverse
double * LUP_solve_inverse(double *A, double *P, int n)
{
    double *A_mirror   = new double[n*n]();
    double *inv_A      = new double[n*n]();
    double *inv_A_each = new double[n](); 
    double *b          = new double[n]();
 
    for (int i = 0; i < n; i++)
    {
        double *L = new double[n*n]();
        double *U = new double[n*n]();
                
        for (int i = 0; i < n; i++)
        {
            b[i] = 0;
        }
        b[i] = 1;
 
        for (int i = 0; i < n * n; i++)
        {
            A_mirror[i] = A[i];
        }

        LUP_Descomposition(A_mirror, L, U, P, n);
        
        inv_A_each = LUP_Solve (L, U, P, b, n);
        
        memcpy(inv_A + i * n, inv_A_each, n*sizeof(double));
    }
    transpose(inv_A,n,n);

    return inv_A;
}
