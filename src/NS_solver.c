#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "project.h"
#include <omp.h>

int Nx = 16;
int Ny = 16;
int Nz = 16;

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

    printf("\n");
    printf("|----------------------------------------------------|\n");
    printf("|--- Running pseudo-spectral Navier-Stokes solver ---|\n");
    printf("|----------------------------------------------------|\n");
    printf("\n");
    printf("Max threads : %d\n", omp_get_max_threads());
    printf("Available processors: %d\n", omp_get_num_procs());

    printf("Domain size is (%f,%f,%f)\n", Lx, Ly, Lz);

    // Allocate the memory for the real 3D fields here
    RealField *v = create_real_field(Nx,Ny,Nz);
    // RealField *omega = create_real_field(Nx,Ny,Nz);

    // FFTW outputs: Complex arrays of size NX * NY * (NZ/2 + 1)
    ComplexField *cv = create_complex_field(Nx,Ny,Nz);    
    // ComplexField *comega = create_complex_field(Nx,Ny,Nz); 

    // Create FFTW plans
    fftw_plan plan_PS = fftw_plan_dft_r2c_3d(Nx, Ny, Nz, v->x, cv->x, FFTW_ESTIMATE);
    fftw_plan plan_SP = fftw_plan_dft_c2r_3d(Nx, Ny, Nz, cv->x, v->x, FFTW_ESTIMATE);

    // Create the wavenumbers
    Wavenumbers *kk = create_wavenumbers(Nx,Ny,Nz,Lx,Ly,Lz);

    // Initialising a Taylor-Green flow
    set_initial_condition(v, init_cond);
    

    // Executing FFT
    execute_fftw_PS(plan_PS, v, cv);

    double *cv_2 = dot_product_c2r(cv,cv);
    double *Ekin = spec1D(cv_2, kk);


    
    for (size_t i = 0; i < kk->Nx; i++){
        printf("k=%2d, E_k=%10.5e\n", i, Ekin[i]);
    }
    

    // compute_curl_fftw(comega,cv,kk);

    // execute_fftw_SP(plan_SP, comega, omega);
    
    // double enstr = (double) dot_product_2_sum(omega,omega) / (N);
    // printf("Enstrophy = %f\n", enstr);

    //  Clean up
    free_real_field(v);
    free_complex_field(cv);

    // free_real_field(omega);
    // free_complex_field(comega);

    fftw_destroy_plan(plan_PS);
    fftw_destroy_plan(plan_SP);

    free_wavenumbers(kk);
    return 0;
}

