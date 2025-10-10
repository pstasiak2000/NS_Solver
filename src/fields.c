#include "Fields.h"
#include <stdlib.h>

// Allocate memory for a 3D Real field
RealField *create_real_field(size_t Nx, size_t Ny, size_t Nz){
    RealField *f = malloc(sizeof(RealField));
    if (!f) return NULL;

    f->Nx = Nx; f->Ny = Ny; f->Nz = Nz;
    f->x = calloc(Nx*Ny*Nz, sizeof(double));
    f->y = calloc(Nx*Ny*Nz, sizeof(double));
    f->z = calloc(Nx*Ny*Nz, sizeof(double));

    if (!f->x || !f->y || !f->z) {
        free(f->x); free(f->y); free(f->z); free(f);
        return NULL;
    }
    return f;
}
// Free the memory from a real field
void free_real_field(RealField *f) {
    if (!f) return;
    free(f->x); free(f->y); free(f->z);
    free(f);
}

ComplexField *create_complex_field(size_t Nx, size_t Ny, size_t Nz){
    ComplexField *f = malloc(sizeof(ComplexField));
    if(!f) return NULL;
        
    f->Nx = Nx; f->Ny = Ny; f->Nz = Nz;
    size_t NzC = Nz/2 + 1;  // FFTW R2C length
    f->x = NULL;
    f->y = NULL;
    f->z = NULL;

    f->x = fftw_malloc(sizeof(fftw_complex) * Nx * Ny *NzC);
    if (!f->x) { free(f); return NULL; }

    f->y = fftw_malloc(sizeof(fftw_complex) * Nx * Ny *NzC);
    if (!f->y) { fftw_free(f->x); free(f); return NULL; }

    f->z = fftw_malloc(sizeof(fftw_complex) * Nx * Ny *NzC);
    if (!f->z) { fftw_free(f->x); fftw_free(f->y); free(f); return NULL; }

    return f;
}

// Free the memory from a complex field
void free_complex_field(ComplexField *f) {
    if (!f) return;
    fftw_free(f->x); fftw_free(f->y); fftw_free(f->z);
    free(f);
}
