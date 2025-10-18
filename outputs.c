#include "outputs.h"
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

double compute_CFL(double v_max, double dx, double dt){
    return (v_max * dt) / dx;
}


void create_file_structure(){
    char path[64];

    sprintf(path,"%sOUTPUTS",dire);
    printf("Generating file structure in %s\n",dire);
    mkdir(path, 0700); //Create outputs folder

    sprintf(path,"%sOUTPUTS/fields",dire);
    mkdir(path, 0700); //Create outputs folder

    sprintf(path,"%sOUTPUTS/spectral",dire);
    mkdir(path, 0700); //Create outputs folder

    sprintf(path,"%sOUTPUTS/global",dire);
    mkdir(path, 0700); //Create outputs folder

    sprintf(path,"cp parameterNS.txt %sOUTPUTS/parameterNS.txt_0",dire);
    system(path);
}

// Save a vector field to binary 
void save_vecfield_2_bin(RealField *f, size_t i){
    char filename[64];
    size_t N = f->Nx * f->Ny * f->Nz;

    snprintf(filename, sizeof(filename), "%s/OUTPUTS/fields/vx.%03zu.dat", dire, i);
    FILE *fvx = fopen(filename, "wb");
    if (!fvx) { perror(filename); return; }
    fwrite(f->x, sizeof(double), N, fvx);
    fclose(fvx);

    snprintf(filename, sizeof(filename), "%s/OUTPUTS/fields/vy.%03zu.dat", dire, i);
    FILE *fvy = fopen(filename, "wb");
    if (!fvy) { perror(filename); return; }
    fwrite(f->y, sizeof(double), N, fvy);
    fclose(fvy);

    snprintf(filename, sizeof(filename), "%s/OUTPUTS/fields/vz.%03zu.dat", dire, i);
    FILE *fvz = fopen(filename, "wb");
    if (!fvz) { perror(filename); return; }
    fwrite(f->z, sizeof(double), N, fvz);
    fclose(fvz);
}

void save_spectrum(const char *filename, double *f,  size_t N){
    FILE *fptr;
    char fullpath[64];

    snprintf(fullpath, sizeof(fullpath), "%sOUTPUTS/spectral/%s",dire,filename);
    fptr = fopen(fullpath, "a");
    for (size_t i = 0; i < N; i++)
        fprintf(fptr, "%e", f[i]);
    fprintf(fptr, "\n");
}