#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>
#include "params.h"
#include "fields.h"


void create_file_structure();

double compute_CFL(double v_max, double dx, double dt);

void save_spectrum(const char *filename, double *f,  size_t N);
void save_vecfield_2_bin(RealField *f, size_t i);


#endif