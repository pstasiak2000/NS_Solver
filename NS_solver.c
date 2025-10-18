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
    printf("--------------------------------------------------------------------------\n");
    printf("|------------- Running pseudo-spectral Navier-Stokes solver -------------|\n");
    printf("--------------------------------------------------------------------------\n");
    printf("\n");
    printf("Max threads : %d\n", omp_get_max_threads());
    printf("Available processors: %d\n", omp_get_num_procs());
    printf("\n");

    read_params("parameterNS.txt");
    create_file_structure();

    double dx = (double) Lx / Nx;
    double dy = (double) Ly / Ny;
    double dz = (double) Lz / Nz;

    double dt_max = (dx*dx) / (2 * nu);
    if(dt>dt_max){
	printf("ERROR: timestep is too large - %f>%f", dt, dt_max);
	abort();
    }

    printf("Timestep passed: dt = %f < %f\n", dt, dt_max);

    // Allocate the memory for the real 3D fields here
    RealField *v = create_real_field(Nx,Ny,Nz);
    RealField *tv = create_real_field(Nx,Ny,Nz);

    // FFTW outputs: Complex arrays of size NX * NY * (NZ/2 + 1)
    ComplexField *cv = create_complex_field(Nx,Ny,Nz);     
    ComplexField *ctv = create_complex_field(Nx,Ny,Nz);
    ComplexField *ctv_rk1 = create_complex_field(Nx,Ny,Nz);

    // Create the wavenumbers
    Wavenumbers *kk = create_wavenumbers(Nx,Ny,Nz,Lx,Ly,Lz);

    // Initialising a Taylor-Green flow
    set_initial_condition(v, init_cond);
    execute_fftw_PS(kk->plan_PS, v, cv);

    printf("Executing loop...\n");
    printf("--------------------------------------------------------------------------------------\n");
    printf("| it |   t   |   v_max  |   v_avg  |   Re   |   eps   |    eta  |   CFL  | dt/dt_max |\n");
    printf("--------------------------------------------------------------------------------------\n");
    int it_shots = 0;
    for (size_t it = 0; it <= steps; it++)
    {
        execute_fftw_PS(kk->plan_PS,v,cv);

        if(it % shots == 0){
            double *v_mag = sqrt_field(dot_product_r2r(v,v),N);
            
            double v_max = max(v_mag,N);
            double v_avg = mean(v_mag,N);
            free(v_mag);

            double CFL_x = compute_CFL(v_max, dx, dt);
            double CFL_y = compute_CFL(v_max, dy, dt);
            double CFL_z = compute_CFL(v_max, dz, dt);

            double CFL = max((double[]){CFL_x, CFL_y, CFL_z},3);

            double Reyn = v_avg * mean((double[]){Lx,Ly,Lz},3) / nu;


            // Save spectrum data
            double *Ekin3D = dot_product_c2r(cv,cv);
            double *spEk = spec1D(Ekin3D, kk);
            save_spectrum("spEk.dat", spEk, Nx);
            double eps = compute_dissipation_spectral(spEk, kk);
            free(Ekin3D); free(spEk);
            
            double eta = pow(pow(nu,3)/eps,0.25); // Computes the Kolmogorov length scale
            dt_max = max((double[]){dx,dy,dz},3) / v_avg; // Computes the max timestep based on avg vel


            printf("| %2d | %5.3f | %8.5f | %8.5f | %6.1f | %6.5f | %6.5f | %5.4f |   %6.4f  | \n", it_shots, t, v_max, v_avg, Reyn, eps, eta, CFL, dt/dt_max);

            // Save global quantities
            save_global("CFL.dat", CFL, t);   // Save the CFL output
            save_global("Reyn.dat", Reyn, t); // Save the Reynolds number out
            save_global("eps.dat", eps, t);   // Save the CFL output
            save_global("eta.dat", eta, t); // Save the Reynolds number out

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
    printf("--------------------------------------------------------------------------\n");
    

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
    free_complex_field(ctv_rk1);

    free_wavenumbers(kk);

    printf("Run finished successfully!\n");
    return 0;
}

