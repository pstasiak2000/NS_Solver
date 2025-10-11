#include "IC.h"
#include <math.h>
#include <stdio.h>

void set_initial_condition(RealField *f, int init_cond){
    switch(init_cond) {
        case 1:
            printf("Taylor-Green flow initial condition\n");
            init_TG(f);
            break;
        default:
            printf("Unknown initial condition: %d\n", init_cond);
            break;
    }

}

// Creates a Taylor-Green initial vector field
void init_TG(RealField *f){
    int Nx = f->Nx;
    int Ny = f->Ny;
    int Nz = f->Nz;

    for (size_t i = 0; i < Nx; i++){
        for (size_t j = 0; j < Ny; j++){
            for (size_t k = 0; k < Nz; k++){
                int idx = (i*Ny + j) * Nz + k;
                double x = (double) i / Nx;
                double y = (double) j / Ny;
                double z = (double) k / Nz;
                f->x[idx] = sin(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
                f->y[idx] = -cos(2*M_PI*x)*sin(2*M_PI*y)*cos(2*M_PI*z);
                f->z[idx] = 0.0;
            }
        }
    }
}

