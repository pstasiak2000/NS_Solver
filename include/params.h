#ifndef PARAMS_H
#define PARAMS_H

extern int Nx;
extern int Ny;
extern int Nz;

// Box size, in units of 2*pi
extern double Lx;
extern double Ly;
extern double Lz;

extern double dt;
/*  Set the initial condition of the vortex flow
    [1]  --- Taylor-Green flow
    [2]  --- ABC flow
*/
extern int init_cond; // Set the inital condition of NS here
extern double nu; // Numerical viscosity
extern int shots; // Save every shots

// Function to read parameters from file
void read_params(const char* filename);

#endif