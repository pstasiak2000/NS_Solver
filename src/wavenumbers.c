#include "wavenumbers.h"
#include <math.h>
#include <omp.h>


#define TWO_PI 6.283185307179586

void execute_fftw_PS(fftw_plan plan, RealField *f, ComplexField *cf){
    fftw_execute_dft_r2c(plan, f->x, cf->x);
    fftw_execute_dft_r2c(plan, f->y, cf->y);
    fftw_execute_dft_r2c(plan, f->z, cf->z);
}

void execute_fftw_SP(fftw_plan plan, ComplexField *cf, RealField *f){
    size_t Nx = f->Nx;
    size_t Ny = f->Ny;
    size_t Nz = f->Nz;

    fftw_execute_dft_c2r(plan, cf->x, f->x);
    fftw_execute_dft_c2r(plan, cf->y, f->y);
    fftw_execute_dft_c2r(plan, cf->z, f->z);

    #pragma omp parallel for 
    for (size_t i = 0; i < Nx*Ny*Nz; i++){
        f->x[i] /= Nx*Ny*Nz;
        f->y[i] /= Nx*Ny*Nz;
        f->z[i] /= Nx*Ny*Nz;
    }
}

Wavenumbers *create_wavenumbers(size_t Nx, size_t Ny, size_t Nz, double Lx, double Ly, double Lz){
    Wavenumbers *kk = malloc(sizeof(Wavenumbers));
    if (!kk) return NULL; //Checks if there is enough memory to allocate the wavenumber structure, if not return NULL

    kk->Nx = Nx;
    kk->Ny = Ny;
    kk->Nz = Nz/2+1;

    kk->kx = malloc(kk->Nx * sizeof(double));
    kk->ky = malloc(kk->Ny * sizeof(double));
    kk->kz = malloc(kk->Nz * sizeof(double));

    if(!kk->kx || !kk->ky || !kk->kz) {
         free(kk->kx); free(kk->ky); free(kk->kz); free(kk);
         return NULL;  
    } //Check if there is enough memory to allocate the arrays, if not return NULL

    for (size_t i = 0; i < kk->Nx; i++) // Allocate kx
        kk->kx[i] = (i <= (kk->Nx)/2) ? (TWO_PI * i / Lx) : (TWO_PI * ( (int) i - (int) kk->Nx) / Lx);
        

    for (size_t i = 0; i < kk->Ny; i++) // Allocate ky
        kk->ky[i] = (i <= (kk->Ny)/2) ? (TWO_PI * i / Ly) : (TWO_PI * ((int) i - (int) kk->Ny) / Ly);
        
    for (size_t i = 0; i < kk->Nz; i++) // Allocate kz
        kk->kz[i] = TWO_PI * i / Lz;
    

    return kk;
}
// Free wavenumbers from memory
void free_wavenumbers(Wavenumbers *kk){
    if(!kk) return;
    free(kk->kx); free(kk->ky); free(kk->kz);
    free(kk);
}

// Compute the 1D spectrum from 3D field in Fourier space (sum over concentric shells)
double *spec1D(double *cf, Wavenumbers *kk){
    double *spec = calloc(kk->Nx, sizeof(double));
    int Nx = kk->Nx;
    int Ny = kk->Ny;
    int Nz = kk->Nz;

    for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
    for (size_t k = 0; k < Nz; k++){
        int idx = (i*Ny + j) * Nz + k;
        double fac = 2.0;
        if(k==0) fac = 1.0;

        double kx = kk->kx[i];
        double ky = kk->ky[j];
        double kz = kk->kz[k];
        
        double ks = pow((kx*kx) + (ky*ky) + (kz*kz),0.5);
        int iks = fmin(floor(ks)+1,Nz);

        spec[iks] += (double) 0.5 * fac *cf[idx];
    }
    

    return spec;
}