#include "wavenumbers.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TWO_PI 6.283185307179586

Wavenumbers *create_wavenumbers(size_t Nx, size_t Ny, size_t Nz, double Lx, double Ly, double Lz){
    Wavenumbers *kk = malloc(sizeof(Wavenumbers));
    if (!kk) return NULL; //Checks if there is enough memory to allocate the wavenumber structure, if not return NULL

    kk->Nx = Nx;
    kk->Ny = Ny;
    kk->Nz = Nz/2+1;

    kk->kx = malloc(kk->Nx * sizeof(double));
    kk->ky = malloc(kk->Ny * sizeof(double));
    kk->kz = malloc(kk->Nz * sizeof(double));

    
    for (size_t i = 0; i < kk->Nx; i++) // Allocate kx
        kk->kx[i] = (i <= (kk->Nx)/2) ? (TWO_PI * i / Lx) : (TWO_PI * (i - Nx) / Lx);

    for (size_t i = 0; i < kk->Nx; i++) // Allocate ky
        kk->ky[i] = (i <= (kk->Ny)/2) ? (TWO_PI * i / Ly) : (TWO_PI * (i - Ny) / Ly);
        
    for (size_t i = 0; i < kk->Nz; i++) // Allocate kz
        kk->kz[i] = TWO_PI * i / Lz;
    
    if(!kk->kx || !kk->ky || !kk->kz) {
         free(kk->kx); free(kk->ky); free(kk->kz); free(kk);
         return NULL;  
    } //Check if there is enough memory to allocate the arrays, if not return NULL

    return kk;
}

void free_wavenumbers(Wavenumbers *kk){
    if(!kk) return;
    free(kk->kx); free(kk->ky); free(kk->kz);
    free(kk);
}


