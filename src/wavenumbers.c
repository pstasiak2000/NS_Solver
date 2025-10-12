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

    for (size_t i = 0; i < kk->Nx; i++){ // Allocate kx
        kk->kx[i] = (i <= (kk->Nx)/2) ? (TWO_PI * i / Lx) : (TWO_PI * ( (int) i - (int) kk->Nx) / Lx);
        
    }

    for (size_t i = 0; i < kk->Ny; i++) // Allocate ky
        kk->ky[i] = (i <= (kk->Ny)/2) ? (TWO_PI * i / Ly) : (TWO_PI * ((int) i - (int) kk->Ny) / Ly);
        
    for (size_t i = 0; i < kk->Nz; i++) // Allocate kz
        kk->kz[i] = TWO_PI * i / Lz;
    

    return kk;
}

void free_wavenumbers(Wavenumbers *kk){
    if(!kk) return;
    free(kk->kx); free(kk->ky); free(kk->kz);
    free(kk);
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