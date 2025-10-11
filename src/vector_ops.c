#include <stdlib.h>
#include "Fields.h"


// Compute the dot product of real fields
double *dot_product_real(RealField *A, RealField *B){
    size_t Nx = A->Nx;
    size_t Ny = A->Ny;
    size_t Nz = A->Nz;
    double* C = calloc(Nx*Ny*Nz, sizeof(double));
    if (!C) return NULL;

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