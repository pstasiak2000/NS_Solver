#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>
#include "params.h"
#include "fields.h"


void create_file_structure();

void save_vecfield_2_bin(RealField *f, size_t i);


#endif