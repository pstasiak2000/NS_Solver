#ifndef PARAMS_H
#define PARAMS_H

extern char dire[];

extern int Nx;
extern int Ny;
extern int Nz;

// Box size, in units of 2*pi
extern double Lx;
extern double Ly;
extern double Lz;

extern int steps; // Total number of steps
extern double dt;
extern int shots; // Save every shots

/*  Set the initial condition of the vortex flow
    [1]  --- Taylor-Green flow
    [2]  --- ABC flow
*/
extern int init_cond; // Set the inital condition of NS here
extern double nu; // Numerical viscosity



// Function to read parameters from file
void read_params(const char* filename);

#endif