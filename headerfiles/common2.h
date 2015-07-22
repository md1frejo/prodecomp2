/* common header file for all other files */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "midmatrix2.h"
//#include "fastnnls.h"

#include <gsl/gsl_linalg.h>

#define PRINT printf("****************************************\n");
#define B exit(0);
#define P(X) printf("\t%d\n",X);
#define F(X) printf("\t%f\n",X);
#define S(X) printf("%s\n",X);
#define C printf("\t----------------------------------\n");
#define P2(X,Y) printf("\t%d %d\n",X,Y);
#define P3(X,Y,Z) printf("\t%d %d %d\n",X,Y,Z);
#define MAT(X) gsl_matrix_fprintf(stdout,X,"%f");
#define VEC(X) gsl_vector_fprintf(stdout,X,"%f");
#define MATE(M,X,Y) printf("%f\n",gsl_matrix_get(M,X,Y));
#define VECE(M,X) printf("%f\n",gsl_vector_get(M,X)); 
#define SIZEM(X) printf("%dx%d\n",(int) X->size1,(int) X->size2);
#define SIZEV(X) printf("%dx1\n",(int) X->size);
#define MMAX(X) printf("%f\n",gsl_matrix_max(X));
#define VMAX(X) printf("%f\n",gsl_vector_max(X));
#define VMAX2(K,X) printf("%d %f\n",K,gsl_vector_max(X));
#define VMAXS2(K,X) printf("%s %f\n",K,gsl_vector_max(X));
#define VMAXDS3(I,K,X) printf("%d %s %f\n",I,K,gsl_vector_max(X));
#define MMIN(X) printf("%f\n",gsl_matrix_min(X));
#define VMIN(X) printf("%f\n",gsl_vector_min(X));
#define OK printf("ok\n");
