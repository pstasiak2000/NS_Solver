#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "project.h"

#define Nx 16
#define Ny 16
#define Nz 16

// Box size, in units of 2*pi
double Lx = 1.0;
double Ly = 1.0;
double Lz = 1.0;

/*  Set the initial condition of the vortex flow
    [1]  --- Taylor-Green flow
*/
int init_cond = 1; // Set the inital condition of NS here

int main() {
    
    int N = Nx * Ny * Nz;

    Lx *= 2*M_PI;
    Ly *= 2*M_PI;
    Lz *= 2*M_PI;
    printf("Domain size is (%f,%f,%f)\n", Lx, Ly, Lz);

    // Allocate the memory for the real 3D fields here
    RealField *v = create_real_field(Nx,Ny,Nz);

    // FFTW outputs: Complex arrays of size NX * NY * (NZ/2 + 1)
    ComplexField *cv = create_complex_field(Nx,Ny,Nz);    

    // Create FFTW plans
    // fftw_plan plan_vx = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, vx, cx, FFTW_ESTIMATE);
    // fftw_plan plan_vy = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, vy, cx, FFTW_ESTIMATE);
    // fftw_plan plan_vz = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, vz, cx, FFTW_ESTIMATE);

    // Create the wavenumbers
    Wavenumbers *kk = create_wavenumbers(Nx,Ny,Nz,Lx,Ly,Lz);

    // Initialising a Taylor-Green flow
    set_initial_condition(v, init_cond);
    

    // Executing FFTs
    // fftw_execute(plan_vx);
    // fftw_execute(plan_vy);
    // fftw_execute(plan_vz);


    double *Ekin = dot_product_real(v,v);


    
    //  Clean up
    free_real_field(v);
    free_complex_field(cv);

    // fftw_destroy_plan(plan_vx);
    // fftw_destroy_plan(plan_vy);
    // fftw_destroy_plan(plan_vz);

    free_wavenumbers(kk);
    return 0;
}

