#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#define Nx 16
#define Ny 16
#define Nz 16

double Lx = 1.0;
double Ly = 1.0;
double Lz = 1.0;

void init_TG(double *vx, double *vy, double *vz);
void create_wavenumbers(double *kk, size_t N);

int main() {
    int N = Nx * Ny * Nz;
    Lx *= 2*M_PI;
    Ly *= 2*M_PI;
    Lz *= 2*M_PI;
    printf("Domain size is (%f,%f,%f)\n", Lx, Ly, Lz);

    // Allocate the memory for the real 3D fields here
    double *vx = fftw_malloc(sizeof(double) * N);
    double *vy = fftw_malloc(sizeof(double) * N);
    double *vz = fftw_malloc(sizeof(double) * N);

    // Wavenumbers
    double dkx = 2.0 * M_PI / Lx;
    double dky = 2.0 * M_PI / Ly;
    double dkz = 2.0 * M_PI / Lz;
    double kx[Nx], ky[Ny], kz[Nz/2+1];

    // FFTW outputs: Complex arrays of size NX * NY * (NZ/2 + 1)
    fftw_complex *cx = fftw_malloc(sizeof(fftw_complex) * Nx * Ny * (Nz/2 + 1));
    fftw_complex *cy = fftw_malloc(sizeof(fftw_complex) * Nx * Ny * (Nz/2 + 1));
    fftw_complex *cz = fftw_malloc(sizeof(fftw_complex) * Nx * Ny * (Nz/2 + 1));

    // Create FFTW plans
    fftw_plan plan_vx = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, vx, cx, FFTW_ESTIMATE);
    fftw_plan plan_vy = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, vy, cx, FFTW_ESTIMATE);
    fftw_plan plan_vz = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, vz, cx, FFTW_ESTIMATE);


    // Creates the wavenumber grid 
    create_wavenumbers(kx, Nx);
    create_wavenumbers(ky, Ny);
    create_wavenumbers(kz, Nz/2+1);
    

    // // Initialising a Taylor-Green flow
    // init_TG(vx,vy,vz);

    // // Executing FFTs
    // fftw_execute(plan_vx);
    // fftw_execute(plan_vy);
    // fftw_execute(plan_vz);

    //

    //  Clean up
    fftw_free(vx); fftw_free(vy); fftw_free(vz);
    fftw_free(cx); fftw_free(cy); fftw_free(cz);

    fftw_destroy_plan(plan_vx);
    fftw_destroy_plan(plan_vy);
    fftw_destroy_plan(plan_vz);
    return 0;
}


void create_wavenumbers(double *kk, size_t N){
    for (size_t i = 0; i < N; i++){
        if (i <= N/2)
            kk[i] = (double) i;
        else
            kk[i] = (double) i - N;
    }
}

void init_TG(double *vx, double *vy, double *vz){
    for (size_t i = 0; i < Nx; i++){
        for (size_t j = 0; j < Ny; j++){
            for (size_t k = 0; k < Nz; k++){
                int idx = (i*Ny + j) * Nz + k;
                double x = (double) i / Nx;
                double y = (double) j / Ny;
                double z = (double) k / Nz;
                vx[idx] = sin(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
                vy[idx] = -cos(2*M_PI*x)*sin(2*M_PI*y)*cos(2*M_PI*z);
                vz[idx] = 0.0;
            }
        }
    }
}