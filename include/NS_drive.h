#ifndef NS_DRIVE_H
#define NS_DRIVE_H

#include <fftw3.h>
#include "params.h"
#include "fields.h"
#include "field_ops.h"
#include "wavenumbers.h"

// Compute the non-linear advection of a single component of velocity


// Compute the non-linear advection in fourier space
void TransportVel(ComplexField *cv, ComplexField *ctv, RealField *v, Wavenumbers *kk);

// Computes the transport of one-component
void Transport(fftw_complex *Advec, fftw_complex *field, RealField *v, Wavenumbers *kk);

// Project a complex field into the the space of divergence free functions
void Proj(ComplexField *cv, Wavenumbers *kk);

void add_visc(ComplexField *ctv, ComplexField *cv, Wavenumbers *kk);

#endif