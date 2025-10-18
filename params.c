#include "params.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

// Default values
char dire[64] = "./";
int Nx = 32, Ny = 32, Nz = 32;
double Lx = 1.0, Ly = 1.0, Lz = 1.0;
double nu = 0.04;
double dt = 0.0001;
int steps = 1;
int init_cond = 0;
int shots = 1;


void read_params(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Warning: could not open %s, using defaults.\n", filename);
        return;
    }

    char line[256];
    while (fgets(line, sizeof(line), file)) {
            // Remove comments safely (only if "//" follows whitespace)
            char* comment = strstr(line, "//");
            if (comment && (comment == line || isspace((unsigned char)*(comment - 1))))
                *comment = '\0';

        // Skip empty lines
        char* p = line;
        while (*p == ' ' || *p == '\t') p++;
        if (*p == '\0') continue;

        // Parse integer parameters
        if (sscanf(line, "dire=%s", dire) == 1) continue;
        if (sscanf(line, "Nx=%d", &Nx) == 1) continue;
        if (sscanf(line, "Ny=%d", &Ny) == 1) continue;
        if (sscanf(line, "Nz=%d", &Nz) == 1) continue;
        if (sscanf(line, "steps=%d", &steps) == 1) continue;
        if (sscanf(line, "init_cond=%d", &init_cond) == 1) continue;
        if (sscanf(line, "shots=%d", &shots) == 1) continue;

        // Parse double parameters
        if (sscanf(line, "Lx=%lf", &Lx) == 1) continue;
        if (sscanf(line, "Ly=%lf", &Ly) == 1) continue;
        if (sscanf(line, "Lz=%lf", &Lz) == 1) continue;
        if (sscanf(line, "nu=%lf", &nu) == 1) continue;
        if (sscanf(line, "dt=%lf", &dt) == 1) continue;
    }
    fclose(file);

    Lx *= 2*M_PI;
    Ly *= 2*M_PI;
    Lz *= 2*M_PI;
    
    

    printf("=== PARAMETERS ===\n");
    printf("Box resolution is (%d,%d,%d)\n", Nx, Ny, Nz);
    printf("Domain size is (%f,%f,%f)\n", Lx, Ly, Lz);
    printf("Timestep size is %f\n", dt);
    printf("The kinematic viscosity is %f\n", nu);
    printf("Saving every %d timesteps\n", shots);
    printf("\n");
}

