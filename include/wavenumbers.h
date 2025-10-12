#ifndef WAVENUMBERS
#define WAVENUMBERS

#include "fields.h"
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>

// Struct to hold wavenumber arrays
typedef struct {
    size_t Nx;
    size_t Ny;
    size_t Nz;
    double *kx;
    double *ky;
    double *kz;
} Wavenumbers;

//FFTs for fields
void execute_fftw_PS(fftw_plan plan, RealField *f, ComplexField *cf);
void execute_fftw_SP(fftw_plan plan, ComplexField *cf, RealField *f);

//Wavenumbers
Wavenumbers *create_wavenumbers(size_t Nx, size_t Ny, size_t Nz, double Lx, double Ly, double Lz);
void free_wavenumbers(Wavenumbers *kk);


double *spec1D(double *cf, Wavenumbers *kk);
#endif