#include "field_ops.h"
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include "params.h"

// Finds the maximum value of a field
double max(double *f, size_t N){
    if (N == 0) return 0.0; // avoid empty array

    double a = f[0];
    #pragma omp parallel for reduction(max:a)
    for (size_t i = 0; i < N; i++)
        a = (a < f[i]) ? f[i] : a;
    return a;
}

// Find the average value of the field
double mean(double *f, size_t N){
    double C = 0;
    #pragma omp parallel for reduction(+:C)
    for (size_t i = 0; i < N; i++)
        C += f[i];
    return C /= N;
}

double *sqrt_field(double *f, size_t N){
    double *f_sqrt = malloc(sizeof(double) * N);

    #pragma omp parallel for
    for (size_t i = 0; i < N; i++)
        f_sqrt[i] = sqrt(f[i]);
    

    return f_sqrt;
}   

// Input two real fields and compute the dot product 
double *dot_product_r2r(RealField *A, RealField *B){
    size_t Nx = A->Nx;
    size_t Ny = A->Ny;
    size_t Nz = A->Nz;
    double* C = calloc(Nx*Ny*Nz, sizeof(double));
    if (!C) return NULL;

    #pragma omp parallel for collapse(3)
    for (size_t i = 0; i < Nx; i++){
        for (size_t j=0; j < Ny; j++){
            for (size_t k=0; k < Nz; k++){
                int idx = (i*Ny + j) * Nz + k;
                C[idx] += A->x[idx] * B->x[idx];
                C[idx] += A->y[idx] * B->y[idx];
                C[idx] += A->z[idx] * B->z[idx];
            }
        }
    }
    return C;
}

// Input two real fields and compute the dot product, expect a real output
double *dot_product_c2r(ComplexField *A, ComplexField *B){
    size_t Nx = A->Nx;
    size_t Ny = A->Ny;
    size_t Nz = A->Nz;
    
    double* C = calloc(Nx*Ny*Nz, sizeof(double));
    if (!C) return NULL;

    #pragma omp parallel for collapse(3)
    for (size_t i = 0; i < Nx; i++){
        for (size_t j=0; j < Ny; j++){
            for (size_t k=0; k < Nz; k++){
                int idx = (i*Ny + j) * Nz + k;

                C[idx] += (A->x[idx][0] * B->x[idx][0]) + ((A->x[idx][1] * B->x[idx][1]));
                C[idx] += (A->y[idx][0] * B->y[idx][0]) + ((A->y[idx][1] * B->y[idx][1]));
                C[idx] += (A->z[idx][0] * B->z[idx][0]) + ((A->z[idx][1] * B->z[idx][1]));
            }
        }
    }
    return C;
}

// Compute the sum of a 3D scalar field
double sum_3D_field(double *A, size_t Nx, size_t Ny, size_t Nz){
    double C = 0.0;

    #pragma omp parallel for collapse(3) reduction(+:C)
    for (size_t i = 0; i < Nx; i++){
        for (size_t j=0; j < Ny; j++){
            for (size_t k=0; k < Nz; k++){
                int idx = (i*Ny + j) * Nz + k;
                C += A[idx];
            }
        }
    }
    return C;
}

// Compute the dot product and sum it (like to compute the kinetic energy)
double dot_product_2_sum(RealField *A, RealField *B){
    size_t Nx = A->Nx;
    size_t Ny = A->Ny;
    size_t Nz = A->Nz;
    double C = 0.0;

    #pragma omp parallel for collapse(3) reduction(+:C)
    for (size_t i = 0; i < Nx; i++){
        for (size_t j=0; j < Ny; j++){
            for (size_t k=0; k < Nz; k++){
                int idx = (i*Ny + j) * Nz + k;
                C += A->x[idx] * B->x[idx];
                C += A->y[idx] * B->y[idx];
                C += A->z[idx] * B->z[idx];
            }
        }
    }
    return C;
}

// Compute the curl of a vector field in Fourier Space
void compute_curl_fftw(ComplexField *comega, ComplexField *cv, Wavenumbers *kk){

    fftw_complex *vx = cv->x;
    fftw_complex *vy = cv->y;
    fftw_complex *vz = cv->z;

    #pragma omp parallel for collapse(3)
    for (size_t i = 0; i < kk->Nx; i++)
    for (size_t j = 0; j < kk->Ny; j++)
    for (size_t k = 0; k < kk->Nz; k++){
        int idx = (i*kk->Ny + j) * kk->Nz + k;

        double kx = kk->kx[i];
        double ky = kk->ky[j];
        double kz = kk->kz[k];

        comega->x[idx][0] =  ky * vz[idx][1] - kz * vy[idx][1];
        comega->x[idx][1] = -ky * vz[idx][0] + kz * vy[idx][0];

        comega->y[idx][0] =  kz * vx[idx][1] - kx * vz[idx][1];
        comega->y[idx][1] = -kz * vx[idx][0] + kx * vz[idx][0];

        comega->z[idx][0] =  kx * vy[idx][1] - ky * vx[idx][1];
        comega->z[idx][1] = -kx * vy[idx][0] + ky * vx[idx][0];
    }
}

// Compute the gradient of f_x
void Dx(fftw_complex *Ain, fftw_complex *Aout, Wavenumbers *kk){

    #pragma omp parallel for collapse(3)
    for (size_t i = 0; i < kk->Nx; i++)
    for (size_t j = 0; j < kk->Ny; j++)
    for (size_t k = 0; k < kk->Nz; k++){
        int idx = (i*kk->Ny + j) * kk->Nz + k;

        double kx = kk->kx[i];
    
        Aout[idx][0] = -kx * Ain[idx][1];
        Aout[idx][1] =  kx * Ain[idx][0];
    }
}

// Compute the gradient of f_x
void Dy(fftw_complex *Ain, fftw_complex *Aout, Wavenumbers *kk){

    #pragma omp parallel for collapse(3)
    for (size_t i = 0; i < kk->Nx; i++)
    for (size_t j = 0; j < kk->Ny; j++)
    for (size_t k = 0; k < kk->Nz; k++){
        int idx = (i*kk->Ny + j) * kk->Nz + k;

        double ky = kk->ky[j];
    
        Aout[idx][0] = -ky * Ain[idx][1];
        Aout[idx][1] =  ky * Ain[idx][0];
    }
}

// Compute the gradient of f_x
void Dz(fftw_complex *Ain, fftw_complex *Aout, Wavenumbers *kk){

    #pragma omp parallel for collapse(3)
    for (size_t i = 0; i < kk->Nx; i++)
    for (size_t j = 0; j < kk->Ny; j++)
    for (size_t k = 0; k < kk->Nz; k++){
        int idx = (i*kk->Ny + j) * kk->Nz + k;

        double kz = kk->ky[k];
    
        Aout[idx][0] = -kz * Ain[idx][1];
        Aout[idx][1] =  kz * Ain[idx][0];
    }
}
// Use the kinetic energy spectrum to compute the mean energy dissipation rate
double compute_dissipation_spectral(double *spEk, Wavenumbers *kk){
    double eps = 0.0;
    for (size_t ik = 0; ik < kk->Nx; ik++)
        eps += 2.0 * nu * (ik*ik) * spEk[ik];
    
    return eps;
}