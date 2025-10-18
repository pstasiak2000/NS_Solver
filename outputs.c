#include "outputs.h"
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

void create_file_structure(){
    char path[64];

    sprintf(path,"%sOUTPUTS",dire);
    printf("Generating file structure in %s\n",dire);
    mkdir(path, 0700); //Create outputs folder

    sprintf(path,"%sOUTPUTS/fields",dire);
    mkdir(path, 0700); //Create outputs folder

    sprintf(path,"%sOUTPUTS/spectral",dire);
    mkdir(path, 0700); //Create outputs folder

    sprintf(path,"cp parameterNS.txt %sOUTPUTS/parameterNS.txt_0",dire);
    system(path);
}

// Save a vector field to binary 
void save_vecfield_2_bin(RealField *f, size_t i){
    char filename[32];
    size_t N = f->Nx * f->Ny * f->Nz;

    snprintf(filename, sizeof(filename), "./fields/vx.%03zu.dat", i);
    FILE *fvx = fopen(filename, "wb");
    if (!fvx) { perror(filename); return; }
    fwrite(f->x, sizeof(double), N, fvx);
    fclose(fvx);

    snprintf(filename, sizeof(filename), "./fields/vy.%03zu.dat", i);
    FILE *fvy = fopen(filename, "wb");
    if (!fvy) { perror(filename); return; }
    fwrite(f->y, sizeof(double), N, fvy);
    fclose(fvy);

    snprintf(filename, sizeof(filename), "./fields/vz.%03zu.dat", i);
    FILE *fvz = fopen(filename, "wb");
    if (!fvz) { perror(filename); return; }
    fwrite(f->z, sizeof(double), N, fvz);
    fclose(fvz);
}
