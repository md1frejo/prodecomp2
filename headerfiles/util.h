#include <gsl/gsl_linalg.h>

#ifndef __incl__util
#define __incl__util

#define ABS(x) (((x) < 0) ? -(x) : (x))
#define MIN(x,y) (((x) < (y) ? (x) : (y))
#define VECTOR_GET(v, i) ((v)->data[(i)*(v)->stride])
#define VECTOR_SET(v, i, x) (v)->data[(i)*(v)->stride] = (x)
#define MATRIX_GET(m, i, j) ((m)->data[(i)*(m)->tda + (j)])
#define MATRIX_SET(m, i, j, x) (m)->data[(i)*(m)->tda + (j)] = (x)

float geteps(void);

// imporvement: linked list instead of ZZ list
typedef struct node {
  int x;
  struct node *next;
} linklist;

extern void sign_matrix(gsl_matrix *sign_MTM,gsl_matrix *MTM);

extern void sign_vector(gsl_vector*,gsl_vector*);

/* matlab: norm(matrix, 1) */
extern double normMatrix(gsl_matrix *matrix);

/* matlab: any(vector(indices) > tol) as boolean */
extern int anyGtTolVector(gsl_vector *vector, int n, int *indices, double tol);

extern int anyGtTolVector2(gsl_vector *vector, int n, int *indices, double tol);
extern int anyGtTolVector3(gsl_vector *vector, int n, gsl_vector *indices, double tol);

/* matlab: any(vector(indices) <= tol) as boolean */
extern int anyLeTolVector(gsl_vector *vector, int n, int *indices, double tol);

/* matlab: indices(max(vector(indices))[1]), i.e. index, not value of max */
/* note: this returns indices[maxi], not maxi */
int maxIndexVector(gsl_vector *vector, int n, int *indices);
int maxIndexVector2(gsl_vector *vector, int n, int *indices);
int maxIndexVector3(gsl_vector *vector, int n,gsl_vector *indices);

/* this finds nonzero elements in array and puts result in nonZeroArray */
/* note: nonZeroArray must be allocated memory (of length at least nonZero) */
void findNonZero(int n, int *array, int *nonZero, int *nonZeroArray);

/* this finds nonnegative elements in array and puts result in nonNegArray */
/* note: nonNegArray must be allocated memory (of length at least nonNeg) */
void findNonNeg(int n, int *array, int *nonNeg, int *nonNegArray);
void findNonNeg2(int n, int *array, int *nonNeg, int *nonNegArray);
void findNonNeg3(gsl_vector *,gsl_vector *,int,int *);

int anyL(int*,int);
int anyV(gsl_vector *,int);
int anyVtol(gsl_vector *list,gsl_vector *PP,int n,double tol);

/* matrix1 = matrix2 * matrix3 */
/* matrix1 must be already allocated */
extern void matrixMatrixMultiply(gsl_matrix *matrix1, gsl_matrix *matrix2, gsl_matrix *matrix3);

/* BLAS level 3 approach  */
extern void matrixMatrixMultiply2(gsl_matrix *matrix1, gsl_matrix *matrix2, gsl_matrix *matrix3);

/* vector1 = matrix * vector2 */
/* vector1 must already be allocated */
extern void matrixVectorMultiply(gsl_vector *vector1, gsl_matrix *matrix, gsl_vector *vector2);

int row_sum(gsl_matrix*,int,int,int);

double difference(float*** projs,gsl_matrix *defs,gsl_matrix *fdir,float*** f,gsl_matrix *ff,int,int,int,int,int);

void calcprojections(float*** f,gsl_matrix *defs,int op,int comps,gsl_matrix *ff,float*** Fout,int mp);

void old_kron(gsl_matrix* M3, gsl_matrix *fdir, gsl_matrix *M1,int mp,int comps,int np);

void kronecker(gsl_matrix* M3, gsl_matrix *fdir, gsl_matrix *M1);

void col_kronecker(gsl_matrix*, gsl_matrix*,gsl_vector*,int,int);

void kronecker_sub(gsl_matrix_view* M3, gsl_matrix_view* fdir, gsl_matrix_view* M1);

void kron1(gsl_matrix* M3, gsl_matrix *fdir, gsl_matrix *M1);

void kron2(gsl_matrix* M3, gsl_matrix *fdir, gsl_matrix *M1);

void kron3(gsl_matrix* M3, gsl_matrix *fdir, gsl_matrix *M1);

void kron4(gsl_matrix* M3, gsl_matrix *fdir, gsl_matrix *M1);

void write_result_indirect(gsl_matrix *matrix);

void write_result_direct(gsl_matrix *matrix);

void open_spectrum(gsl_matrix *spectra3,char* spec_name);

void slice(gsl_matrix* spectra3,gsl_matrix *slice,int start,int end);

void normalizem(gsl_matrix* spectra);

void normalize_spec(gsl_matrix**,int,float);

void normalizev(gsl_vector* spectra);

void normalize(float*** spectra,int x,int y,int z);

int max_check(float*** spectra,int x,int y,int z,float norm);

int max_check_spectra(gsl_matrix**,int,float);

void max_value_m(gsl_matrix* spec);

void max_value_v(gsl_vector* spec);

void max_value(float*** spectra,int x,int y,int z);

float sumf(gsl_matrix** f,int mp,int m2,int m1);

float sumfdir(gsl_matrix* fdir,int np,int m1);

void matrix_mul(float** a,int row,int col,float** b);

void generate_test_spectra(float*** projs,int mp,int direct,int op);

int M1_check(gsl_matrix*,int);

void permutation(gsl_matrix*); 

#endif /* __incl__util */
