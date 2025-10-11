#ifndef FIELDS_H
#define FIELDS_H

#include <fftw3.h>

// Real-valued field
typedef struct {
    size_t Nx, Ny, Nz;
    double *x;
    double *y;
    double *z;
} RealField;

RealField *create_real_field(size_t Nx, size_t Ny, size_t Nz);
void free_real_field(RealField *f);

// Complex-valued field for fft
typedef struct {
    size_t Nx, Ny, Nz;
    fftw_complex *x;
    fftw_complex *y;
    fftw_complex *z;
} ComplexField;

ComplexField *create_complex_field(size_t Nx, size_t Ny, size_t Nz);
void free_complex_field(ComplexField *f);

#endif // !FIELDS_H

#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

#include <fftw3.h>
#include <wavenumbers.h>

double *dot_product_real(RealField *A, RealField *B);

#endif