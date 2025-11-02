#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>     /*  for  fabs, sqrt, pow, exp, sin, cos, log,   */
                      /*       atan, acos                             */
#include <string.h>
#include <algorithm>
#include <time.h>

using namespace std;
const int n = 16;

// CUDA kernel to preform Jacobi
//__global__ void Jacobi(float *A, int N, float *D, float *V, int *NROT, float s, float tau, float *B, float *Z) 
__global__ void Jacobi(int N, float *D, float *B, float *Z)
{
    //int idy = threadIdx.x + blockDim.x * blockIdx.x;
    //int idy = threadIdx.y + blockDim.y * blockIdx.y;
    //int idy = threadIdx.y;
    int idy = threadIdx.x;
    //int idx = threadIdx.x;          
   
    if (idy < N)
    {
        B[idy] += Z[idy];
        D[idy] = B[idy];
        Z[idy] = 0.0;
    }                     
        
    return;  //too many iterations
}    

int cmp(float a,float b)
{
    return a>b;
}

void Sorting(float *A, int N, float *D, float *V) 
{
    float   *V1 = new float [N*N];           // eigenvectors solution matrix
    float   *D1 = new float [N];              // eigenvalues solution vector
    int     *nu = new int [N];
    int     i, j;
     
    for (i = 0; i < N; i++)
    {
        D1[i] = D[i];
    }

    sort(D, D + N, cmp);
    
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (D[i] == D1[j])
            {
                nu[i] = j;
            }
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            V1[j * N + i] = V[j * N + nu[i]];
        }
    }

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            V[i * N + j] = V1[i * N + j];
        }
    }

    delete [] D1;
    delete [] V1;
    return;
}

int main()  {
    
    int i, j, ip, iq;                         // loop variables
    float vmax, sm, tresh, c, g, h, s, t, tau, theta;                     // normalization value

    int nrot;                        // Number of matrix rotations
    int iter;                        // Number of iterations
    
    //allocate vector and matrices on host
    float   *Mat = new float [n*n];                    // given square symmetric matrix
    float   *Mat1 = new float [n*n];
    float   *Eigenvectors = new float [n*n];           // eigenvectors solution matrix
    float   *Eigenvalues = new float [n];              // eigenvalues solution vector
    
    float   *B = new float [n];
    float   *Z = new float [n];
              
    //allocate vectors B, Z
    float   *d_B;
    float   *d_Z;
    //float   *d_Mat;
    float   *d_Eigenvalues;
    //float   *d_Eigenvectors;
    
    // Allocate Unified Memory -- accessible from CPU or GPU
    //cudaMalloc((float**)&d_Mat, n*n*sizeof(float *));
    //cudaMalloc((float**)&d_Eigenvectors, n*n*sizeof(float *));

    cudaMalloc((float**)&d_Eigenvalues, n*sizeof(float *));
    cudaMalloc((float**)&d_B, n*sizeof(float *));    
    cudaMalloc((float**)&d_Z, n*sizeof(float *));
    
    //cudaMalloc((int**)&d_nrot, sizeof(int *));


    /*float Rxx[][4] = {{1.0000,    0.5000,    0.3333,    0.2500},
                     {0.5000,    1.0000,    0.6667,    0.5000},
                     {0.3333,    0.6667,    1.0000,    0.7500},
                     {0.2500,    0.5000,    0.7500,    1.0000}};*/
    float Rxx[][16] = {
                        {3.7450,    0.0839,    0.8675,    0.8287,    2.0374,    1.7158,   -0.8208,    2.3104,         0,    0.6090,    0.8836,   -2.2513,    1.1669,    0.1956,    0.6991,   -1.3373},
                        {0.0839,    3.0862,    0.2318,    0.5619,    0.5974,    1.7734,    1.5466,   -0.7386,   -0.6090,         0,    0.6258,    0.8013,   -1.6870,    0.8037,    0.3392,    0.6059},
                        {0.8675,    0.2318,    3.3715,    0.0973,    0.8360,    0.8869,    1.9295,    1.5585,   -0.8836,   -0.6258,         0,    0.3838,    0.7301,   -1.9622,    1.0107,   -0.0538},
                        {0.8287,    0.5619,    0.0973,    3.5665,   -0.0063,    0.6810,    0.6569,    2.0806,    2.2513,   -0.8013,   -0.3838,         0,    0.8428,    0.8335,    -1.8880,    1.2303},
                        {2.0374,    0.5974,    0.8360,   -0.0063,    3.1479,    0.2583,    0.6119,    0.7910,   -1.1669,    1.6870,   -0.7301,   -0.8428,         0,    0.3315,    0.7888,   -2.0354},
                        {1.7158,    1.7734,    0.8869,    0.6810,    0.2583,    3.2225,    0.1897,    0.8974,   -0.1956,   -0.8037,    1.9622,   -0.8335,   -0.3315,         0,    0.7787,    0.5938},
                       {-0.8208,    1.5466,    1.9295,    0.6569,    0.6119,    0.1897,    3.3317,    0.0532,   -0.6991,   -0.3392,   -1.0107,    1.8880,   -0.7888,   -0.7787,         0,    0.5647},    
                        {2.3104,   -0.7386,    1.5585,    2.0806,    0.7910,    0.8974,    0.0532,    3.5365,    1.3373,   -0.6059,    0.0538,   -1.2303,    2.0354,   -0.5938,    -0.5647,         0},
                        {     0,   -0.6090,   -0.8836,    2.2513,   -1.1669,   -0.1956,   -0.6991,    1.3373,    3.7450,    0.0839,    0.8675,    0.8287,    2.0374,    1.7158,    -0.8208,    2.3104},
                        {0.6090,         0,   -0.6258,   -0.8013,    1.6870,   -0.8037,   -0.3392,   -0.6059,    0.0839,    3.0862,    0.2318,    0.5619,    0.5974,    1.7734,    1.5466,   -0.7386},
                        {0.8836,    0.6258,         0,   -0.3838,   -0.7301,    1.9622,   -1.0107,    0.0538,    0.8675,    0.2318,    3.3715,    0.0973,    0.8360,    0.8869,    1.9295,    1.5585},
                       {-2.2513,    0.8013,    0.3838,         0,   -0.8428,   -0.8335,    1.8880,   -1.2303,    0.8287,    0.5619,    0.0973,    3.5665,   -0.0063,    0.6810,    0.6569,    2.0806},
                        {1.1669,   -1.6870,    0.7301,    0.8428,         0,   -0.3315,   -0.7888,    2.0354,    2.0374,    0.5974,    0.8360,   -0.0063,    3.1479,    0.2583,    0.6119,    0.7910},
                        {0.1956,    0.8037,   -1.9622,    0.8335,    0.3315,         0,   -0.7787,   -0.5938,    1.7158,    1.7734,    0.8869,    0.6810,    0.2583,    3.2225,    0.1897,    0.8974},
                        {0.6991,    0.3392,    1.0107,   -1.8880,    0.7888,    0.7787,         0,   -0.5647,   -0.8208,    1.5466,    1.9295,    0.6569,    0.6119,    0.1897,    3.3317,    0.0532},
                       {-1.3373,    0.6059,   -0.0538,    1.2303,   -2.0354,    0.5938,    0.5647,         0,    2.3104,   -0.7386,    1.5585,    2.0806,    0.7910,    0.8974,    0.0532,    3.5365}
                        };

    // Flattening 2D-data to 1D-data
    for (i = 0; i < n; i++) 
    {
        for (j = 0; j < n; j++)
        {
            Mat[(i * n) + j] = Rxx[i][j];
            Mat1[(i * n) + j] = Mat[(i * n) + j];
        }
    }

    // Initialize Eigenvectors to identity matrix
    for (i = 0; i < n; i++)         
    {
        for (j = 0; j < n; j++)
        {
            Eigenvectors[i * n + j] = 0.0; 
        }
        Eigenvectors[i * n + i] = 1.0;
    } 

    for (i = 0; i < n; i++) 
    {
        B[i] = Mat[i * n + i];
        Eigenvalues[i] = B[i];
        Z[i] = 0.0;
    }  
    
    nrot = 0;
    
    clock_t t1;
    t1=clock();
    for (iter = 0; iter < 50; iter++)
    {
        //iter = 50;
        sm = 0;
            
        for (i = 0; i < n - 1; i++)    //sum off-diagonal elements
        {            
            for (j = i + 1; j < n; j++)
            {
                sm = sm + fabs(Mat1[(i * n) + j]);
            }
            //sm = 0;
        }

        if (sm == 0)
        {
            break;       //normal return
            //exit (0);            
        }

        if (iter < 4)
        {
            tresh = 0.2 * sm * sm;
            //tresh = 0.2 * sm / (n * n);
        }
        else
        {
            tresh = 0.0;
        }
        
        for (ip = 0; ip < n - 1; ip++)        
        {
            for (iq = ip + 1; iq < n; iq++)            
            {
                g = 100 * fabs(Mat1[ip * n + iq]);            // after 4 sweeps, skip the rotation if the off-diagonal element is small
                //cout << g << "\t";
                //cout << endl;
                if ((i > 4) && (fabs(Eigenvalues[ip]) + g == fabs(Eigenvalues[ip])) && (fabs(Eigenvalues[iq]) + g == fabs(Eigenvalues[iq])))
                {
                    Mat1[ip * n + iq] = 0.0;
                }
                else if (fabs(Mat1[ip * n + iq]) > tresh)
                {
                    h = Eigenvalues[iq] - Eigenvalues[ip];
                    //cout << h << "\t";
                    //cout << endl;
                    if (fabs(h) + g == fabs(h))
                    {
                        t = (Mat1[ip * n + iq]) / h;
                    }
                    else
                    {
                        theta = 0.5 * h / Mat1[ip * n + iq];
                        t = 1 / (fabs(theta) + sqrt(1.0 + theta * theta));
                        //cout << t << "\t";
                        //cout << endl;
                        if (theta < 0)
                        {
                            t = -t;
                        }
                    }

                    c = 1.0 / sqrt(1.0 + t * t);
                    s = t * c;
                    tau = s / (1.0 + c);
                    h = t * Mat1[ip * n + iq];
                    
                    Z[ip] -= h;
                    Z[iq] += h;
                    Eigenvalues[ip] -= h;
                    Eigenvalues[iq] += h;
                    Mat1[ip * n + iq] = 0.0;
                    //cout << s << "\t";
                    //cout << endl;
                    
                    for (j = 0; j < ip; j++)
                    {
                        g = Mat1[j * n + ip];
                        h = Mat1[j * n + iq];
                        Mat1[j * n + ip] = g - s * (h + g * tau);
                        Mat1[j * n + iq] = h + s * (g - h * tau);
                    }                    
                    for (j = ip + 1; j < iq; j++)
                    {
                        g = Mat1[ip * n + j];
                        h = Mat1[j * n + iq];
                        Mat1[ip * n + j] = g - s * (h + g * tau);
                        Mat1[j * n + iq] = h + s * (g - h * tau);
                    }
                    for (j = iq + 1; j < n; j++)
                    {
                        g = Mat1[ip * n + j];
                        h = Mat1[iq * n + j];
                        Mat1[ip * n + j] = g - s * (h + g * tau);
                        Mat1[iq * n + j] = h + s * (g - h * tau);
                    }
                    for (j = 0; j < n; j++)
                    {
                        g = Eigenvectors[j * n + ip];
                        h = Eigenvectors[j * n + iq];
                        Eigenvectors[j * n + ip] = g - s * (h + g * tau);
                        Eigenvectors[j * n + iq] = h + s * (g - h * tau);
                    }

                    nrot = nrot + 1;                        
                }
            }
        }

        // Copy the data from host to device
        //cudaMemcpy(d_Mat, Mat1, n*n*sizeof(float), cudaMemcpyHostToDevice);                               
        //cudaMemcpy(d_Eigenvectors, Eigenvectors, n*n*sizeof(float), cudaMemcpyHostToDevice);
        
        cudaMemcpy(d_Eigenvalues, Eigenvalues, n*sizeof(float), cudaMemcpyHostToDevice);        
        cudaMemcpy(d_B, B, n*sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_Z, Z, n*sizeof(float), cudaMemcpyHostToDevice);

        //cudaMemcpy(d_s, s, sizeof(float), cudaMemcpyHostToDevice);
                        
        // Launch kernel on 1M elements on the GPU
        dim3 block(16, 1, 1);
        dim3 grid(1, 1, 1);
        
        //Jacobi<< <grid, block>> >(d_Mat, n, d_Eigenvalues, d_Eigenvectors, nrot, d_s, tau, d_B, d_Z);
        Jacobi<< <grid, block>> >(n, d_Eigenvalues, d_B, d_Z);
 
        // Wait for GPU to finish before accessing on host
        cudaDeviceSynchronize();
        
        // Copy the data from device to host
        cudaMemcpy(Eigenvalues, d_Eigenvalues, n*sizeof(float), cudaMemcpyDeviceToHost);        
        
    }
    t1 = clock() - t1;
    printf ("It took me %d clicks (%f seconds).\n",t1,((double)t1)/CLOCKS_PER_SEC);

    //cudaMemcpy(nrot, d_nrot, sizeof(int), cudaMemcpyDeviceToHost);
    
    //Displaying the orignal matrix
    cout << "---------------------------------------------------\n" << endl;
    cout <<  "Symmetric matrix is:" << endl;
    for (i = 0; i < n; i++) 
    {
        for (j = 0; j < n; j++)
        {
            cout << Mat[(i * n) + j] << "\t";
        }
        cout << "\t" << endl;
    }
    cout << "\t" << endl;
    
    Sorting(Mat,n,Eigenvalues,Eigenvectors);  

    //Displaying the results 
    cout << "Eigenvalues:" << endl;
    for (i = 0; i < n; i++)
    {
        cout << Eigenvalues[i] << "\t";
    }
    cout << "\t" << endl;
    cout << "\t" << endl;

    cout << "Eigenvectors (in lines):" << endl;
    for (j = 0; j < n; j++) 
    {
        vmax = Eigenvectors[j * n];
	    for (i = 0; i < n; i++)
	    if(fabs(Eigenvectors[(i * n) + j]) > fabs(vmax))
        {
            vmax = Eigenvectors[(i * n) + j];
        }  
        for (i = 0; i < n; i++)
        {
            cout << Eigenvectors[(i * n) + j] / vmax << "\t";
        }
	    cout << "\t" << endl;
    }
    cout << "\t" << endl;
  
    cout << "Number of rotations:" << nrot << endl;
    cout << "\t" << endl;

     
    // Free memory
    cudaFree(d_B);
    cudaFree(d_Z);
    //cudaFree(d_Mat);
    cudaFree(d_Eigenvalues);
    //cudaFree(d_Eigenvectors);
    
    delete [] Mat;
    delete [] Eigenvalues;
    delete [] Eigenvectors; 
    delete [] B;
    delete [] Z; 

}
