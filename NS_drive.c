#include "NS_drive.h"
#include "field_ops.h"

void EulerStepNS(fftw_complex *fieldOLD, fftw_complex *fieldNew, fftw_complex *fieldNL, double dtt){
    int N = Nx * Ny * (Nz/2+1);

    #pragma omp parallel for 
    for (size_t i = 0; i < N; i++){
        fieldNew[i][0] = fieldOLD[i][0] - dtt * fieldNL[i][0]; 
        fieldNew[i][1] = fieldOLD[i][1] - dtt * fieldNL[i][1]; 
    }

}


void TransportVel(ComplexField *ctv, ComplexField *cv, RealField *v, Wavenumbers *kk){
        
    // Compute the non-linear advection terms in Fourier space
    Transport(ctv->x, cv->x, v, kk);
    Transport(ctv->y, cv->y, v, kk);
    Transport(ctv->z, cv->z, v, kk);

    // Now do the in-place projection
    Proj(ctv, kk);
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

    fftw_complex *dive = fftw_malloc(sizeof(fftw_complex)*Nx*Ny*Nz);

    #pragma omp parallel for collapse(3)
    for (size_t i = 0; i < Nx; i++)
    for (size_t j = 0; j < Ny; j++)
    for (size_t k = 0; k < Nz; k++){
        int idx = (i*Ny + j) * Nz + k;

        double kx = kk->kx[i];
        double ky = kk->ky[j];
        double kz = kk->kz[k];
        double k2 = kx*kx + ky*ky + kz*kz + eps;

        double Re = (kx * cv->x[idx][0]) + (ky * cv->x[idx][0]) + (ky * cv->x[idx][0]);
        double Im = (kx * cv->x[idx][1]) + (ky * cv->x[idx][1]) + (ky * cv->x[idx][1]);

        dive[idx][0] = -Im;
        dive[idx][1] = Re;

        cv->x[idx][0] -= kx * dive[idx][0] / k2;
        cv->x[idx][1] -= kx * dive[idx][1] / k2;

        cv->y[idx][0] -= ky * dive[idx][0] / k2;
        cv->y[idx][1] -= ky * dive[idx][1] / k2;

        cv->z[idx][0] -= kz * dive[idx][0] / k2;
        cv->z[idx][1] -= kz * dive[idx][1] / k2;
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