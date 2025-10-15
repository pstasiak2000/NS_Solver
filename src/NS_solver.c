#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <omp.h>

#include "project.h"



int main() {
    
    int N = Nx * Ny * Nz;

    printf("\n");
    printf("|----------------------------------------------------|\n");
    printf("|--- Running pseudo-spectral Navier-Stokes solver ---|\n");
    printf("|----------------------------------------------------|\n");
    printf("\n");
    printf("Max threads : %d\n", omp_get_max_threads());
    printf("Available processors: %d\n", omp_get_num_procs());

    read_params("parameterNS.txt");

    Lx *= 2*M_PI;
    Ly *= 2*M_PI;
    Lz *= 2*M_PI;

    printf("Domain size is (%f,%f,%f)\n", Lx, Ly, Lz);

    // Allocate the memory for the real 3D fields here
    RealField *v = create_real_field(Nx,Ny,Nz);
    RealField *tv = create_real_field(Nx,Ny,Nz);

    // FFTW outputs: Complex arrays of size NX * NY * (NZ/2 + 1)
    ComplexField *cv = create_complex_field(Nx,Ny,Nz);     
    ComplexField *ctv = create_complex_field(Nx,Ny,Nz);
    // ComplexField *ctv_rk1 = create_complex_field(Nx,Ny,Nz);

    // Create the wavenumbers
    Wavenumbers *kk = create_wavenumbers(Nx,Ny,Nz,Lx,Ly,Lz);

    // Initialising a Taylor-Green flow
    set_initial_condition(v, init_cond);
    execute_fftw_PS(kk->plan_PS, v, cv);

    // 
    printf("---------------------------\n");
    int it_shots = 0;
    for (size_t it = 0; it <= steps; it++)
    {
        execute_fftw_PS(kk->plan_PS,v,cv);

        if(it % shots == 0){    
            double max_value = sqrt(max(dot_product_r2r(v,v),kk->Nx,kk->Ny,kk->Nz));
            double Reyn = max_value * Lx / nu;
            printf("| %2d | v_max = %10.5f | Re = %10.2f\n",it_shots, max_value, Reyn);
            save_vecfield_2_bin(v,it_shots);
            ++it_shots;
        }

        TransportVel(ctv,cv,v,kk);
        add_visc(ctv, cv, kk);        

        EulerStepNS(cv->x, cv->x, ctv->x, dt);
        EulerStepNS(cv->y, cv->y, ctv->y, dt);
        EulerStepNS(cv->z, cv->z, ctv->z, dt);


    // // First Runge-Kutta step
    //     TransportVel(ctv,cv,v,kk);
    //     add_visc(ctv, cv, kk);        

    //     EulerStepNS(cv->x, ctv_rk1->x, ctv->x, 0.5*dt);
    //     EulerStepNS(cv->y, ctv_rk1->y, ctv->y, 0.5*dt);
    //     EulerStepNS(cv->z, ctv_rk1->z, ctv->z, 0.5*dt);

    // // Second Runge-Kutta step
    //     execute_fftw_SP(kk->plan_SP, ctv_rk1, v);

    //     TransportVel(ctv_rk1, ctv, v, kk);
    //     add_visc(ctv, ctv_rk1, kk);

    //     EulerStepNS(cv->x, cv->x, ctv->x, dt);
    //     EulerStepNS(cv->y, cv->y, ctv->y, dt);
    //     EulerStepNS(cv->z, cv->z, ctv->z, dt);

        execute_fftw_SP(kk->plan_SP,cv,v);


    }

    

    // printf("Max value = %f\n", max(dot_product_c2r(cv,cv),kk->Nx,kk->Ny,kk->Nz));

    // double *cv_2 = dot_product_c2r(cv,cv);
    // double *Ekin = spec1D(cv_2, kk);

    // printf("---------------------\n");
    // printf("|  k  |     E(k)    |\n");
    // printf("---------------------\n");
    // for (size_t i = 0; i < kk->Nx; i++){
    //     printf("| %3d | %10.5f |\n", i, (double) Ekin[i]);
    // }
    // printf("---------------------\n");

    // save_vecfield_2_bin(v, 0);
    // printf("Saved vector field to binary!\n");

    //  Clean up
    free_real_field(v); free_real_field(tv);
    free_complex_field(cv); free_complex_field(ctv);

    free_wavenumbers(kk);

    printf("Run finished successfully!\n");
    return 0;
}

