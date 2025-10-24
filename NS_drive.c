#include "NS_drive.h"
#include "field_ops.h"

void RK2Step(ComplexField *cv, ComplexField *ctv, ComplexField *ctv_rk1, RealField *v, Wavenumbers *kk, double dtt){
    size_t N = kk -> Nx * kk->Ny * kk->Nz;

     // --- Stage 1: Compute intermediate state at t + dt/2 ---
    
     // Compute RHS at current state: ctv = RHS(cv)
    TransportVel(ctv, cv, v, kk);
    // add_visc(ctv, cv, kk);

    // Compute intermediate state: ctv_rk1 = cv + 0.5*dt * ctv
    #pragma omp parallel for
    for (size_t i = 0; i < N; i++){
        ctv_rk1->x[i][0] = cv->x[i][0] + 0.5 * dtt * ctv->x[i][0];
        ctv_rk1->x[i][1] = cv->x[i][1] + 0.5 * dtt * ctv->x[i][1];
        
        ctv_rk1->y[i][0] = cv->y[i][0] + 0.5 * dtt * ctv->y[i][0];
        ctv_rk1->y[i][1] = cv->y[i][1] + 0.5 * dtt * ctv->y[i][1];
        
        ctv_rk1->z[i][0] = cv->z[i][0] + 0.5 * dtt * ctv->z[i][0];
        ctv_rk1->z[i][1] = cv->z[i][1] + 0.5 * dtt * ctv->z[i][1];
    }

    // Transform intermediate state to real space for TransportVel
    Proj(ctv_rk1, kk);
    execute_fftw_SP(kk->plan_SP, ctv_rk1, v);


    // --- Stage 2: Compute final state ---

    // Compute RHS at intermediate state: ctv = RHS(ctv_rk1)
    TransportVel(ctv, ctv_rk1, v, kk);
    // add_visc(ctv, ctv_rk1, kk);

    // Final update: cv = cv + dt * ctv
    #pragma omp parallel for
    for (size_t i = 0; i < N; i++) {
        cv->x[i][0] += dtt * ctv->x[i][0];
        cv->x[i][1] += dtt * ctv->x[i][1];
        
        cv->y[i][0] += dtt * ctv->y[i][0];
        cv->y[i][1] += dtt * ctv->y[i][1];
        
        cv->z[i][0] += dtt * ctv->z[i][0];
        cv->z[i][1] += dtt * ctv->z[i][1];
    }
}

void EulerStep(ComplexField *fieldOLD, ComplexField *fieldNew, ComplexField *fieldNL, double dtt){
    int N = Nx * Ny * (Nz/2+1);
    
    #pragma omp parallel for 
    for (size_t i = 0; i < N; i++){
        fieldNew->x[i][0] = fieldOLD->x[i][0] - dtt * fieldNL->x[i][0]; 
        fieldNew->x[i][1] = fieldOLD->x[i][1] - dtt * fieldNL->x[i][1];

        fieldNew->y[i][0] = fieldOLD->y[i][0] - dtt * fieldNL->y[i][0]; 
        fieldNew->y[i][1] = fieldOLD->y[i][1] - dtt * fieldNL->y[i][1]; 

        fieldNew->z[i][0] = fieldOLD->z[i][0] - dtt * fieldNL->z[i][0]; 
        fieldNew->z[i][1] = fieldOLD->z[i][1] - dtt * fieldNL->z[i][1]; 
    }
}


void TransportVel(ComplexField *ctv, ComplexField *cv, RealField *v, Wavenumbers *kk){
        
    // Compute the non-linear advection terms in Fourier space
    Transport(ctv->x, cv->x, v, kk);
    Transport(ctv->y, cv->y, v, kk);
    Transport(ctv->z, cv->z, v, kk);

    // Now do the in-place projection
    // check_divergence(ctv, kk, "Before projection");
    Proj(ctv, kk);
    // check_divergence(ctv, kk, "After projection");
}

// Compute the transport of a single velocity component v_x, v_y or v_z
void Transport(fftw_complex *Advec, fftw_complex *field, RealField *v, Wavenumbers *kk){

    size_t Nx = v->Nx;
    size_t Ny = v->Ny;
    size_t Nz = v->Nz;

    size_t N = Nx*Ny*Nz;

    double *Rtemp = malloc(sizeof(double) * Nx * Ny * Nz);
    double *Atemp = malloc(sizeof(double) * Nx * Ny * Nz);
    fftw_complex *Ctemp = fftw_malloc(sizeof(fftw_complex) * Nx * Ny * (Nz/2+1));


    // Compute the partial derivate wrt x of a velocity component and take the dot product in real space
    Dx(field, Ctemp, kk);
    fftw_execute_dft_c2r(kk->plan_SP, Ctemp, Rtemp);
    #pragma omp parallel for 
    for (size_t i = 0; i < N; i++)
        Atemp[i] = v->x[i] * Rtemp[i];
    
    // Compute the partial derivate wrt y of a velocity component and take the dot product in real space
    Dy(field, Ctemp, kk);
    fftw_execute_dft_c2r(kk->plan_SP, Ctemp, Rtemp);
    #pragma omp parallel for 
    for (size_t i = 0; i < N; i++)
        Atemp[i] += v->y[i] * Rtemp[i];

    // Compute the partial derivate wrt z of a velocity component and take the dot product in real space
    Dz(field, Ctemp, kk);
    fftw_execute_dft_c2r(kk->plan_SP, Ctemp, Rtemp);
    #pragma omp parallel for 
    for (size_t i = 0; i < N; i++)
        Atemp[i] += v->z[i] * Rtemp[i];

    // Take advection term back to Fourier space
    fftw_execute_dft_r2c(kk->plan_PS, Atemp, Advec);

    // Normalize the advected field due to Fourier transform
    for (size_t i = 0; i < (Nx*Ny*kk->Nz); i++){
            Advec[i][0] /= N;
            Advec[i][1] /= N;
    }

    // De-alias the field
    for (size_t i = 0; i < (Nx*Ny*kk->Nz); i++){
        Advec[i][0] *= kk->kill[i];
        Advec[i][1] *= kk->kill[i];
    }


    free(Rtemp); free(Atemp);
    free(Ctemp);
}

// Compute the projection of a complex field
void Proj(ComplexField *cv, Wavenumbers *kk){
    size_t Nx = kk->Nx;
    size_t Ny = kk->Ny;
    size_t Nz = kk->Nz;

    double eps=1e-15;

    #pragma omp parallel for collapse(3)
    for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
    for (size_t k = 0; k < Nz; k++){
        int idx = (i*Ny + j) * Nz + k;
        fftw_complex dive;

        double kx = kk->kx[i];
        double ky = kk->ky[j];
        double kz = kk->kz[k];
        double k2 = kx*kx + ky*ky + kz*kz + eps;

        double Re = (kx * cv->x[idx][0]) + (ky * cv->y[idx][0]) + (kz * cv->z[idx][0]);
        double Im = (kx * cv->x[idx][1]) + (ky * cv->y[idx][1]) + (kz * cv->z[idx][1]);

        dive[0] = Re;
        dive[1] = Im;

        cv->x[idx][0] -= kx * dive[0] / k2;
        cv->x[idx][1] -= kx * dive[1] / k2;

        cv->y[idx][0] -= ky * dive[0] / k2;
        cv->y[idx][1] -= ky * dive[1] / k2;

        cv->z[idx][0] -= kz * dive[0] / k2;
        cv->z[idx][1] -= kz * dive[1] / k2;
    }
    
}

void add_visc(ComplexField *ctv, ComplexField *cv, Wavenumbers *kk){
    #pragma omp parallel for collapse(3)
    for (size_t i = 0; i < kk->Nx; i++)
    for (size_t j = 0; j < kk->Ny; j++)
    for (size_t k = 0; k < kk->Nz; k++){
        int idx = (i*Ny + j) * (kk->Nz) + k;

        double kx = kk->kx[i];
        double ky = kk->ky[j];
        double kz = kk->kz[k];
        double k2 = kx*kx + ky*ky + kz*kz ;

        ctv->x[idx][0] +=  nu * k2 * cv->x[idx][0];
        ctv->x[idx][1] +=  nu * k2 * cv->x[idx][1];

        ctv->y[idx][0] +=  nu * k2 * cv->y[idx][0];
        ctv->y[idx][1] +=  nu * k2 * cv->y[idx][1];

        ctv->z[idx][0] +=  nu * k2 * cv->z[idx][0];
        ctv->z[idx][1] +=  nu * k2 * cv->z[idx][1];
    } 
}