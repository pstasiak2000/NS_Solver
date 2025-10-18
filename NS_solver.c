#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <omp.h>

#include "project.h"



int main() {
    double t = 0.0;
    int N = Nx * Ny * Nz;

    printf("\n");
    printf("|----------------------------------------------------|\n");
    printf("|--- Running pseudo-spectral Navier-Stokes solver ---|\n");
    printf("|----------------------------------------------------|\n");
    printf("\n");
    printf("Max threads : %d\n", omp_get_max_threads());
    printf("Available processors: %d\n", omp_get_num_procs());

    read_params("parameterNS.txt");
    create_file_structure();

    Lx *= 2*M_PI;
    Ly *= 2*M_PI;
    Lz *= 2*M_PI;

    double dx = (double) Lx / Nx;
    double dy = (double) Ly / Ny;
    double dz = (double) Lz / Nz;

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
    printf("-----------------------------------------------------\n");
    printf("| it |   t   |   v_max  |   v_avg  |   CFL  |   Re  |\n");
    printf("-----------------------------------------------------\n");
    int it_shots = 0;
    for (size_t it = 0; it <= steps; it++)
    {
        execute_fftw_PS(kk->plan_PS,v,cv);

        if(it % shots == 0){
            double *v_mag = sqrt_field(dot_product_r2r(v,v),N);

            double v_max = max(v_mag,N);
            double v_avg = mean(v_mag,N);
            double CFL_x = compute_CFL(v_max, dx, dt);
            double CFL_y = compute_CFL(v_max, dy, dt);
            double CFL_z = compute_CFL(v_max, dz, dt);

            double CFL = max((double[]){CFL_x, CFL_y, CFL_z},3);

            double Reyn = v_avg * mean((double[]){Lx,Ly,Lz},3) / nu;
            printf("| %2d | %5.3f | %8.5f | %8.5f | %5.4f | %5.0f |\n", it_shots, t, v_max, v_avg, CFL, Reyn);
            free(v_mag);

            double *Ekin3D = dot_product_c2r(cv,cv);

            // Save spectrum data
            double *spEk = spec1D(Ekin3D, kk);
            save_spectrum("spEk.dat", spEk, Nx);
            free(Ekin3D); free(spEk);
            

            // Save outputs to field
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

        t += dt;
    }
    printf("-----------------------------------------------------\n");
    

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

