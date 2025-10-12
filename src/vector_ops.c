#include <stdlib.h>
#include <omp.h>
#include <complex.h>
#include "Fields.h"
#include "wavenumbers.h"

// Compute the dot product of real fields
double *dot_product_real(RealField *A, RealField *B){
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



