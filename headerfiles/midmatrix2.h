#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include "util.h"
#include <gsl/gsl_blas.h>

void midmatrix2(gsl_matrix**,gsl_matrix*,int*);
//float**** midmatrix2(gsl_matrix*,int*);
void flip_matrix_col(gsl_matrix**,int,int);

