#ifndef __incl__fastnnls
#define __incl__fastnnls

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include "midmatrix2.h"

/* matlab: fastnnls(XtX, Xty, tol) */
extern void fastnnls(gsl_matrix *, gsl_vector *, gsl_vector *, gsl_vector *, double tol);

//void fastnnls(gsl_matrix *XtX, gsl_vector *Xty, gsl_vector *x, gsl_vector *y,double tol);

int row_sum(gsl_matrix*,int,int,int);

#endif /* __incl__fastnnls */
