/* fastnnls.c implementation Jonas F */

// using fixed matrices and vectors with vector views instead of floating sizes

#ifndef MAXFLOAT
#define MAXFLOAT 3.40282347e+38F
#endif

#include "headerfiles/common2.h"
#include "headerfiles/fastnnls.h"
#include "headerfiles/util.h"

//#include <math.h>
//#include <gsl/gsl_linalg.h>
//#include "util.h"

/* Heuristic value  */
//#define EPS  1.0e-12
#define EPS 2.22045E-16 // from numpy.finfo(flaot) 2.2204460492503131e-16

// from eps.c:
//#define EPS 1.19209e-07

/* amd46 redhat: did not find it in /usr/include/math.h */
//#define MAXFLOAT 3.40282347e+38F
//#define MAXFLOAT 3.40282347e+38F

static void solveEquation(int n, gsl_vector *z, int nzp,gsl_vector *PP,
	gsl_permutation *perm, gsl_matrix *XtX, gsl_vector *Xty,
	gsl_matrix *workMat, gsl_vector *workVec) {
  
  int i, j, s, m;
  
  gsl_vector *tres=gsl_vector_calloc(nzp);
  gsl_matrix_view wm=gsl_matrix_submatrix(workMat,0,0,nzp,nzp);
  gsl_vector_view wv=gsl_vector_subvector(workVec,0,nzp);

  perm->size = nzp;
  
  int *inic=(int *) malloc(nzp*sizeof(int));
  for(i=0; i<nzp; i++)
    inic[i]=-1;

  for(i=0,m=0; i<n; i++)
    if(gsl_vector_get(PP,i)!=-1)
      inic[m++]=gsl_vector_get(PP,i);

  for (i=0; i<nzp; i++) {
    for (j=0; j<nzp; j++) 
      gsl_matrix_set(&wm.matrix,i,j,gsl_matrix_get(XtX,*(inic+i),*(inic+j)));
    gsl_vector_set(&wv.vector,i,gsl_vector_get(Xty,*(inic+i)));
  }
  
  gsl_permutation_init(perm);
  gsl_linalg_LU_decomp(&wm.matrix,perm,&s);
  gsl_linalg_LU_solve(&wm.matrix,perm,&wv.vector,tres);

  perm->size=n;

  for (i=0; i<nzp; i++) {
    if(inic[i]!=-1)
      gsl_vector_set(z,inic[i],gsl_vector_get(tres,i));
    else
      printf("\twarning: inic neg: %d\n",inic[i]);
  }
  gsl_vector_free(tres);
  free(inic);
}

// x goes to inf. 
/* static void determineW(gsl_matrix *XtX,gsl_vector *Xty,gsl_vector *x,gsl_vector *w) { */
/*   gsl_blas_dgemv(CblasNoTrans,-1.0,XtX,x,0.0,w); */
/*   gsl_vector_add(w,Xty); */
/* } */

void fastnnls(gsl_matrix *XtX, gsl_vector *Xty, gsl_vector *x, gsl_vector *w,double tol) {

  int i,n, nzz, nzp, iter, itmax, t;
  double v, vv, alpha;
  gsl_vector *z, *workVec;
  gsl_matrix *workMat;
  gsl_permutation *perm;

  //  printf("eps: %e\n",geteps());
  
  n=XtX->size1;

  if(n!=XtX->size2) {
    printf("XtX must be a square matrix (got %d x %d)\n",n,(int) XtX->size2);
    return;
  }

  if (n != Xty->size) {
    printf("XtX must be same size as Xty (got %d vs %d)\n", n,(int) Xty->size);
    return;
  }

  if (n != x->size) {
    printf("XtX must be same size as x (got %d vs %d)\n", n, (int) x->size);
    return;
  }

  if (n != w->size) {
    printf("XtX must be same size as w (got %d vs %d)\n", n, (int) w->size);
    return;
  }

//    float eps=geteps();
  float limit=1.0;
  if (tol <= limit)
    tol=10*EPS*normMatrix(XtX) * XtX->size1; // from function

  //  tol=0.000001;

  gsl_vector *P=gsl_vector_alloc(n);
  gsl_vector *PP=gsl_vector_calloc(n);
  gsl_vector *Z=gsl_vector_calloc(n);
  gsl_vector *ZZ=gsl_vector_calloc(n);

  perm = gsl_permutation_calloc(n);
  z=gsl_vector_alloc(n);
  workVec = gsl_vector_alloc(n);
  workMat = gsl_matrix_alloc(n,n);
  gsl_vector *temp1=gsl_vector_calloc(n);
  gsl_vector *temp2=gsl_vector_calloc(n);

  for(i=0; i<n; i++) {
    gsl_vector_set(Z,i,i);
    gsl_vector_set(ZZ,i,i);
  }
  
  gsl_vector_set_zero(x);
  gsl_vector_set_all(P,-1);
  gsl_vector_set_all(PP,-1);

  nzz=n;

  gsl_blas_dgemv(CblasNoTrans,1.0,XtX,x,0.0,temp1);
  gsl_vector_memcpy(temp2,Xty);
  gsl_vector_sub(temp2,temp1);
  gsl_vector_memcpy(w,temp2);

  iter=0;
  itmax=30*n; // hardcoded

/* matlab: outer loop: while any(Z) && any(w(ZZ) > tol) */
  
  while(anyV(Z,n) && anyGtTolVector3(w,n,ZZ,tol)) {

    t=maxIndexVector3(w,n,ZZ);

    t=gsl_vector_get(ZZ,t);

    gsl_vector_set(P,t,t);
    gsl_vector_set(Z,t,-1);

    gsl_vector_set_all(PP,-1);
    gsl_vector_set_all(ZZ,-1);

    findNonNeg3(PP,P,n,&nzp);
    findNonNeg3(ZZ,Z,n,&nzz);

    if(nzp)
      solveEquation(n,z,nzp,PP,perm,XtX,Xty,workMat,workVec);
    for(i=0; i<n; i++)
      if(gsl_vector_get(ZZ,i)!=-1)
	gsl_vector_set(z,gsl_vector_get(ZZ,i),0);

    /* matlab: inner loop: while any((z(PP) <= tol) & iter < itmax; iter = iter + 1 */
    for(; iter<itmax && anyVtol(z,PP,n,tol); iter++) {
      alpha=MAXFLOAT; /* a bit dangerous */
      
      for(i=0; i<n; i++) {
	vv=gsl_vector_get(z,i);
	if((gsl_vector_get(P,i)!=-1) && (vv<=tol)) {
	  v=gsl_vector_get(x,i);
	  v=v/(v-vv);
	  if(v<alpha) {
	    alpha=v;
	  }
	}
      }

      for(i=0; i<n; i++) {
	v=gsl_vector_get(x,i);
      	vv=gsl_vector_get(z,i);
      	v+=alpha*(vv-v);
      	gsl_vector_set(x,i,v);
      }
      for(i=0; i<n; i++) {
	v=gsl_vector_get(x,i);
	if(fabs(v)<tol && gsl_vector_get(P,i)!=-1) { // note abs
	  gsl_vector_set(Z,i,i);
	  gsl_vector_set(P,i,-1);
	}
      }
      gsl_vector_set_all(PP,-1);
      gsl_vector_set_all(ZZ,-1);

      findNonNeg3(PP,P,n,&nzp);
      findNonNeg3(ZZ,Z,n,&nzz);
      if(nzp)
	solveEquation(n,z,nzp,PP,perm,XtX,Xty,workMat,workVec);
      for(i=0; i<n; i++)
	if(gsl_vector_get(ZZ,i)!=-1)
	  gsl_vector_set(z,gsl_vector_get(ZZ,i),0);
    }
  
    gsl_vector_memcpy(x,z);
    
    gsl_blas_dgemv(CblasNoTrans,1.0,XtX,x,0.0,temp1);
    gsl_vector_memcpy(temp2,Xty);
    gsl_vector_sub(temp2,temp1);
    gsl_vector_memcpy(w,temp2);
  }
  
  gsl_vector_free(workVec);
  gsl_matrix_free(workMat);
  gsl_vector_free(z);
  gsl_vector_free(temp1);
  gsl_vector_free(temp2);

  gsl_permutation_free(perm);

  free(P);
  free(PP);
  free(Z);
  free(ZZ);
}
