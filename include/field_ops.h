#ifndef FIELD_OPS_H
#define FIELD_OPS_H

#include <fftw3.h>
#include "fields.h"
#include "wavenumbers.h"

double *dot_product_r2r(RealField *A, RealField *B);
double *dot_product_c2r(ComplexField *A, ComplexField *B);

double sum_3D_field(double *A, size_t Nx, size_t Ny, size_t Nz);
double dot_product_2_sum(RealField *A, RealField *B);

void compute_curl_fftw(ComplexField *comega, ComplexField *cv, Wavenumbers *kk);

double max(double *f, size_t Nx, size_t Ny, size_t Nz);
#endif