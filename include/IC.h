#ifndef IC
#define IC

// #include <stddef.h>
#include "fields.h"

void set_initial_condition(RealField *f, int init_cond);

// Initial conditions for Taylor-Green flows
void init_TG(RealField *f);
void init_ABC(RealField *f);

#endif // !
