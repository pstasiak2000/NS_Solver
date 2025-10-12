#include "outputs.h"



// Save a vector field to binary 
void save_vecfield_2_bin(RealField *f, size_t i){
    char filename[32];
    size_t N = f->Nx * f->Ny * f->Nz;

    snprintf(filename, sizeof(filename), "../vx.%03zu.dat", i);
    FILE *fvx = fopen(filename, "wb");
    if (!fvx) { perror(filename); return; }
    fwrite(f->x, sizeof(double), N, fvx);
    fclose(fvx);

    snprintf(filename, sizeof(filename), "../vy.%03zu.dat", i);
    FILE *fvy = fopen(filename, "wb");
    if (!fvy) { perror(filename); return; }
    fwrite(f->y, sizeof(double), N, fvy);
    fclose(fvy);

    snprintf(filename, sizeof(filename), "../vz.%03zu.dat", i);
    FILE *fvz = fopen(filename, "wb");
    if (!fvz) { perror(filename); return; }
    fwrite(f->z, sizeof(double), N, fvz);
    fclose(fvz);
}
