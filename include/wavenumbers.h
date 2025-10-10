#ifndef WAVENUMBERS
#define WAVENUMBERS

#include <fftw3.h>

// Struct to hold wavenumber arrays
typedef struct {
    size_t Nx;
    size_t Ny;
    size_t Nz;
    double *kx;
    double *ky;
    double *kz;
} Wavenumbers;

//Wavenumbers
Wavenumbers *create_wavenumbers(size_t Nx, size_t Ny, size_t Nz, double Lx, double Ly, double Lz);
void free_wavenumbers(Wavenumbers *kk);


#endif