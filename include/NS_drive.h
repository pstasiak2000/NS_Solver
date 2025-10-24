#ifndef NS_DRIVE_H
#define NS_DRIVE_H

#include <fftw3.h>
#include "params.h"
#include "fields.h"
#include "field_ops.h"
#include "wavenumbers.h"

// Euler step of field
void EulerStep(ComplexField *fieldOLD, ComplexField *fieldNew, ComplexField *fieldNL, double dtt);

// RK2 step of field
void RK2Step(ComplexField *cv, ComplexField *ctv, ComplexField *ctv_rk1, RealField *v, Wavenumbers *kk, double dtt);

// Compute the non-linear advection in fourier space
void TransportVel(ComplexField *cv, ComplexField *ctv, RealField *v, Wavenumbers *kk);

// Computes the transport of one-component
void Transport(fftw_complex *Advec, fftw_complex *field, RealField *v, Wavenumbers *kk);

// Project a complex field into the the space of divergence free functions
void Proj(ComplexField *cv, Wavenumbers *kk);

void add_visc(ComplexField *ctv, ComplexField *cv, Wavenumbers *kk);

#endif