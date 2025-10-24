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
double compute_dissipation_spectral(double *spEk, Wavenumbers *kk);

void Dx(fftw_complex *Aout, fftw_complex *Ain, Wavenumbers *kk);
void Dy(fftw_complex *Aout, fftw_complex *Ain, Wavenumbers *kk);
void Dz(fftw_complex *Aout, fftw_complex *Ain, Wavenumbers *kk);

void check_divergence(ComplexField *cv, Wavenumbers *kk, const char* label);

double max(double *f, size_t N);
double mean(double *f, size_t N);
double *sqrt_field(double *f, size_t N);

#endif