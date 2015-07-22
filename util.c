/* useful functions */

#include "headerfiles/common2.h"
#include "headerfiles/util.h"

/* matlab: sign(matrix) returns matrix with 1,0,-1 depending on elements sign */
void sign_matrix(gsl_matrix *sign_MTM,gsl_matrix *MTM) {
  double val=0;
  int a;
  int b;

  for (a=0; a<MTM->size1; a++) {
    for (b=0; b<MTM->size2; b++) {
      val=gsl_matrix_get(MTM,a,b);
      if(val>1)  
	gsl_matrix_set(sign_MTM,a,b,1);
      else if(val<0)
	gsl_matrix_set(sign_MTM,a,b,0);
      else
	gsl_matrix_set(sign_MTM,a,b,-1);
    }
  }
}

void sign_vector(gsl_vector *sign_vector_MTM,gsl_vector *diag_MTM) {
  double val=0;
  int a;

    for (a=0; a<diag_MTM->size; a++) {
      val=gsl_vector_get(diag_MTM,a);
      if(val>1)  
	gsl_vector_set(sign_vector_MTM,a,1);
      else if(val<0)
	gsl_vector_set(sign_vector_MTM,a,0);
      else
      gsl_vector_set(sign_vector_MTM,a,-1);
  }
}

/* matlab: norm(matrix, 1) */
double normMatrix(gsl_matrix *matrix)
{
    int i, j;
    float maxVal, colSum;
    //   printf("XtX max: %f\n",gsl_matrix_max(matrix));
    maxVal = 0;
    for (i = 0; i < matrix->size2; i++)
    {
	colSum = 0;
        for (j = 0; j < matrix->size1; j++)
	{
	    colSum += ABS(gsl_matrix_get(matrix, j, i));
	}

	if (colSum > maxVal)
	    maxVal = colSum;
    }
    return maxVal;
}

/* matlab: any(vector(indices) > tol) as boolean */
int anyGtTolVector(gsl_vector *vector, int n, int *indices, double tol)
{
    int i;

    for (i = 0; i < n; i++)
    {
	if (gsl_vector_get(vector, indices[i]) > tol)
	    return 1;
    }

    return 0;
}

int anyGtTolVector2(gsl_vector *vector, int n, int *indices, double tol) {

    int i;

    for (i=0; i<n; i++) {
      if(gsl_vector_get(vector,indices[i])>tol) // note: -1 -> 0
	return 1;
    }
    return 0;
}

int anyGtTolVector3(gsl_vector *vector,int n,gsl_vector *ZZ,double tol) {

  int i,k;

    for (i=0; i<n; i++) {
      k=gsl_vector_get(ZZ,i);
      if(k!=-1 && gsl_vector_get(vector,k)>tol)
	return 1;
    }
    return 0;
}

/* matlab: any(vector(indices) <= tol) as boolean */
int anyLeTolVector(gsl_vector *vector, int n, int *indices, double tol)
{
    int i;

    for (i = 0; i < n; i++)
      {
	if (gsl_vector_get(vector, indices[i]) <= tol)
	  return 1;
      }
    
    return 0;
}

/* matlab: indices(max(vector(indices))[1]), i.e. index, not value of max */
int maxIndexVector(gsl_vector *vector, int n, int *indices)
{
    int i, maxi;
    double maxv, v;

    if (n < 1)
	return -1;

    maxi = indices[0];
    maxv = gsl_vector_get(vector, indices[0]);

    for (i = 1; i < n; i++)
    {
	v = gsl_vector_get(vector, indices[i]);
	if (v > maxv)
	{
	    maxi = indices[i];
	    maxv = v;
	}
    }

    return maxi;
}

// jf
int maxIndexVector2(gsl_vector *vector, int n, int *indices) {
    int i, maxi;
    double maxv;

    if (n < 1)
	return -1;

    maxi=indices[0];
    maxv=gsl_vector_get(vector,indices[0]);

    for (i=0; i<n; i++) {
      if(gsl_vector_get(vector,indices[i])>maxv) {
	maxi=indices[i];
	maxv=gsl_vector_get(vector,indices[i]);
      }
    }
    return maxi;
}

int maxIndexVector3(gsl_vector *vector,int n,gsl_vector *indices) {

  int i,maxi,k=0;
  double maxv;
  
  if(n<1)
    return -1;
  
  maxi=0;
  maxv=-1;

  for (i=0; i<n; i++) {
    k=gsl_vector_get(indices,i);
    if(k!=-1 && gsl_vector_get(vector,k)>maxv) {
      maxi=k;
      maxv=gsl_vector_get(vector,k);
    }
  }
  return maxi;
}

void findNonZero(int n, int *array, int *nonZero, int *nonZeroArray)
{
    int i, m;

    m = 0;
    //    #pragma omp parallel default(shared) private(i)
    //    {
    for (i = 0; i < n; i++)    
      if (array[i])
	nonZeroArray[m++] = array[i];
    //    }
    *nonZero = m;
}

void findNonNeg(int n,int *array,int *nonNeg,int *nonNegArray) {
    int i, m;

    m = 0;
    for (i = 0; i < n; i++)
    {
      if (array[i] >= 0)
	    nonNegArray[m++] = array[i];
    }

    *nonNeg = m;
}

void findNonNeg2(int n,int *array,int *nonNeg,int *nonNegArray) {
  
  int i,m;
  
  m=0;
  for(i=0; i<n; i++) {
    if (array[i]>0) // note: because we skipped -1 we have to use >0
      nonNegArray[m++]=array[i];
  }
  *nonNeg=m;
}

void findNonNeg3(gsl_vector *tt,gsl_vector *t,int n,int *nNeg) {
  
  int i,m=0;
  
  //  gsl_vector_set_zero(tt);

  for(i=0; i<n; i++)
    if(gsl_vector_get(t,i)!=-1) {
      gsl_vector_set(tt,i,i);
      m++;
    }
  *nNeg=m;
}


int anyL(int *list,int n)  {
  
  int i; 

  for(i=0; i<n; i++)
    if(list[i])
      return 1;

  return 0;
}

int anyV(gsl_vector *list,int n)  {
  
  int i; 

  for(i=0; i<n; i++)
    if(gsl_vector_get(list,i)!=-1)
      return 1;

  return 0;
}

int anyVtol(gsl_vector *z,gsl_vector *PP,int n,double tol)  {
  
  int i,k; 

  for(i=0; i<n; i++) {
    k=gsl_vector_get(PP,i);
    //    printf("i: %d inic[i]: %d z[i] %f tol: %f\n",i,inic[i],gsl_vector_get(z,inic[i]),tol);
    if(k!=-1 && gsl_vector_get(z,k)<=tol)
      return 1;
  }
  return 0;
}

/* matrix1 = matrix2 * matrix3 */
/* matrix1 must be already allocated */
/* TBD: put checks in to see that matrices have consistent sizes */

void matrixMatrixMultiply(gsl_matrix *matrix1, gsl_matrix  *matrix2, gsl_matrix *matrix3) {

  double alpha = 1.0, beta = 0.0;

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, alpha, matrix2,  
		 matrix3, beta, matrix1); 

}

/* vector1 = matrix * vector2 */
/* vector1 must already be allocated */
/* TBD: put checks in to see that have consistent sizes */

void matrixVectorMultiply(gsl_vector *vector1, gsl_matrix *matrix, gsl_vector *vector2)
{
    gsl_blas_dgemv(CblasNoTrans, 1.0, matrix, vector2, 0.0, vector1);
}

int row_sum(gsl_matrix *defs,int j,int n,int nd) {
  int i=0;
  int rsum=0;
  //  for(i=n; i<nd; i++)    
  for(i=n; i<=(nd-1); i++)    
    rsum+=ABS(gsl_matrix_get(defs,j,i));
  return rsum;   
}

double difference(float*** projs,gsl_matrix *defs,gsl_matrix *fdir,float*** f,gsl_matrix *ff,int op,int nd,int np,int mp,int comps) {

  int q=fdir->size2;
  int i,j,k,l,z=19;
  gsl_vector *sum=gsl_vector_alloc(q);

  float*** Pout;
  float**** Ptemp;

  float targetvalue=0.123123;

  Pout=calloc(mp,sizeof(float**)); 
  
  for(i=0; i<mp; i++)
    Pout[i]=calloc(mp,sizeof(float*));

  for(i=0; i<mp; i++)   
    for(j=0; j<mp; j++)   
      Pout[i][j]=calloc(z,sizeof(float));

  //  calcprojections(f,defs,op,comps,ff,Pout,mp);

  Ptemp=calloc(np,sizeof(float***)); 
  
  for(i=0; i<np; i++)
    Ptemp[i]=calloc(comps,sizeof(float**));

  for(i=0; i<np; i++)   
    for(j=0; j<comps; j++)   
      Ptemp[i][j]=calloc(op,sizeof(float*));

 for(i=0; i<np; i++)   
    for(j=0; j<comps; j++)   
      for(k=0; k<op; k++)
	Ptemp[i][j][k]=calloc(q,sizeof(float));

 for(i=0; i<np; i++)
   for(j=0; j<comps; j++)
     for(k=0; k<op; k++)
       for(l=0; l<q; l++)
	 Ptemp[i][j][k][l]=Pout[k][i][j]*gsl_matrix_get(ff,k,l);

 // ...P=sum(P,4);
 for(i=0; i<q; q++)
   gsl_vector_set(sum,i,Ptemp[0][0][0][q]);

 // ...targetvalue=sum(sum(sum(abs(P-projs))));

  return targetvalue;
}

void kronecker(gsl_matrix* CC, gsl_matrix *AA, gsl_matrix *BB) {
  int m,n,p,q,x=0,y=0;
  for(m=0; m<AA->size1; m++) {
    for(p=0; p<BB->size1; p++) {
      for(n=0; n<AA->size2; n++) { 
	for(q=0; q<BB->size2; q++) {
	  gsl_matrix_set(CC,y,x,gsl_matrix_get(AA,m,n)*gsl_matrix_get(BB,p,q));
	  x++;
	  if(x==CC->size2) {
	    x=0;
	    y++;
	  }
	}
      }
    }
  }
}

void col_kronecker(gsl_matrix *AA, gsl_matrix *BB,gsl_vector *column,int n,int q) {
  int m,p,x=0;
  for(m=0; m<AA->size1; m++) 
    for(p=0; p<BB->size1; p++) 
      gsl_vector_set(column,x++,gsl_matrix_get(AA,m,n)*gsl_matrix_get(BB,p,q));
}

void kronecker_sub(gsl_matrix_view *M3, gsl_matrix_view *fdir, gsl_matrix_view *M1) {
}


void write_result_indirect(gsl_matrix *matrix) {
  int a;
  FILE *f;

  for(a=0; a<matrix->size2; a++) {
    gsl_matrix_view indirect_shape=gsl_matrix_submatrix(matrix,0,a,matrix->size1,1);
    char string[50];
    sprintf(string,"output/indirect.%d",a);
    f=fopen(string,"w");
    gsl_matrix_fprintf(f,&indirect_shape.matrix,"%f");
    fclose(f);  
  }
}

void write_result_direct(gsl_matrix *matrix) {
  int a;
  FILE *f;

  for(a=0; a<matrix->size2; a++) {
    gsl_matrix_view direct_shape=gsl_matrix_submatrix(matrix,0,a,matrix->size1,1);
    char string[50];
    sprintf(string,"output/direct.%d",a);
    f=fopen(string,"w");
    gsl_matrix_fprintf(f,&direct_shape.matrix,"%f");
    fclose(f);  
  }
}

// beta version
void open_spectrum(gsl_matrix* spectra,char* spec_name) {
  FILE *f;
  gsl_matrix *temp=gsl_matrix_alloc(spectra->size2,spectra->size1);
  f=fopen(spec_name,"r");
  gsl_matrix_fscanf(f,temp);
  fclose(f);
  gsl_matrix_transpose_memcpy(spectra,temp);
  gsl_matrix_free(temp);
}

void slice(gsl_matrix* spectra,gsl_matrix *cut,int start,int end) {
  int a,b,c;
  for(c=0,a=start-1; a<end; a++,c++) 
    for(b=0; b<spectra->size1; b++) 
      gsl_matrix_set(cut,b,c,gsl_matrix_get(spectra,b,a));
}

void normalizem(gsl_matrix* spectra) {
  double max;
  max=gsl_matrix_max(spectra);
  max=0.1/max;
  gsl_matrix_scale(spectra,max);
}

void normalize_spec(gsl_matrix** spectra,int op,float factor) {
  //  long double max=0,temp;
  int a;
  float max=0, temp=0;
  for(a=0; a<op; a++) {
    //  temp=(gsl_matrix_max(spectra[a]));
    temp=(gsl_matrix_max(spectra[a]));
    if(max<temp) 
      max=temp;
  }
  max=(float) factor/max;
  for(a=0; a<op; a++) 
    gsl_matrix_scale(spectra[a],max);
}

void normalizev(gsl_vector* spectra) {
  double max;
  max=gsl_vector_max(spectra);
  max=1.0/max;
  gsl_vector_scale(spectra,max);
}

void normalize(float*** spectra,int x,int y,int z) {
  int a,b,c;
  double scale=0;
  for(a=0; a<x; a++) {
    for(b=0; b<y; b++) {
      for(c=0; c<z; c++) {
	if(spectra[a][b][c]>=scale)
	  scale=spectra[a][b][c];
      }
    }
  }

  scale=(1/scale);

  for(a=0; a<x; a++) 
    for(b=0; b<y; b++) 
      for(c=0; c<z; c++) 
	spectra[a][b][c]*=scale;
}  

int max_check(float*** spectra,int x,int y,int z,float norm) {
  int a,b,c;
  for(a=0; a<x; a++) 
    for(b=0; b<y; b++) 
      for(c=0; c<z; c++) 
	if(spectra[a][b][c]>norm)
	   return 0;
  return 1;
}

int max_check_spectra(gsl_matrix** spectra,int op,float norm) {
  printf("%f\n",norm);
  int a,b,c;
  for(a=0; a<op; a++)
    for(b=0; b<spectra[a]->size1; b++) 
      for(c=0; c<spectra[a]->size2; c++) 
	if(gsl_matrix_get(spectra[a],b,c)>norm+0.00001)
	  return 0;
  return 1;
}


void max_value(float*** spectra,int x,int y,int z) {
  int a,b,c;
  double val=0.0;
  for(a=0; a<x; a++) 
    for(b=0; b<y; b++) 
      for(c=0; c<z; c++) 
	if(spectra[a][b][c]>val)
	  val=spectra[a][b][c];
  printf("max value: %f\n",val);
}

void max_value_m(gsl_matrix* spec) {
  printf("%f\n",gsl_matrix_max(spec));
}

void max_value_v(gsl_vector* spec) {
  printf("%f\n",gsl_vector_max(spec));
}
   
float sumf(gsl_matrix** f,int mp,int m2,int m1) {
  int a;
  float c=0;
  for(a=0; a<mp; a++)
    c+=gsl_matrix_get(f[m1],a,m2);
  return c;
}

float sumfdir(gsl_matrix* fdir,int np,int m1) {
  int a;
  float c=0;
  for(a=0; a<np; a++)
    c+=gsl_matrix_get(fdir,a,m1);
  return c;
}

void generate_test_spectra(float*** projs,int mp,int direct,int op) {
  int a,b,c;
  srand((unsigned)time(NULL));
  for(a=0; a<mp; a++)  
    for(b=0; b<direct; b++)  
      for(c=0; c<op; c++) 
	projs[a][b][c]=(float)rand()/((float)(RAND_MAX)+(float)(1));
}

void matrix_mul(float** a,int row,int col,float** b) {
  //  int i,j,n;
  //  float sum=0;

  // alternative solution
/*   for(jjj=0; jjj<comps; jjj++) {  */
/*     for(a=0; a<mp; a++) {  */
/*       for(b=0; b<comps; b++) {  */
/* 	for(c=0; c<mp; c++) {  */
/* 	  sum+=Fout[a][c][jjj]*f[c][jj][jjj];  */
/* 	  //		printf("%f\n",sum);  */
/* 	}  */
/* 	gsl_matrix_set(ff,i,jjj,sum);  */
/*       }  */
/*     }  */
/*   }  */
  
  // alternative solution
}

int M1_check(gsl_matrix* M1,int comps) {
  int x,y,j;
  for(j=0; j<comps; j++) {
    for(y=0,x=j*M1->size2/comps; x<(j+1)*M1->size2/comps; y++,x++) {
      if(gsl_matrix_get(M1,y,x)>1.0) 
	return 0;
    }
  }
  return 1;
}

/* void permutation(gsl_matrix* MTM) { */
/*   int a,b; */
/*   float temp; */
/*   int x=MTM->size1; */
/*   gsl_vector *col=gsl_vector_alloc(x); */
/*   for(a=0; a<x-1; a++) { */
/*     gsl_matrix_get_col(col,MTM,a); */
/*     temp=gsl_vector_get(col,x-1); */
/*     for(b=1; b<x; b++) */
/*       gsl_matrix_set(MTM,b,a+1,gsl_vector_get(col,b-1)); */
/*     gsl_matrix_set(MTM,0,a+1,temp); */
/*   } */
/*   gsl_vector_free(col); */
/* } */

void permutation(gsl_matrix* MTM) {
  int a,b=0;
  float temp=0;
  int x=MTM->size1;
  for(b=0; b<x; b++) {
    temp=gsl_matrix_get(MTM,b,0);    
    for(a=0; a<x-b; a++) {
      gsl_matrix_set(MTM,a,a+b,temp);
      gsl_matrix_set(MTM,a+b,a,temp);
    }
  }
}

// void read_spectrum(char* specname, )

float geteps(void) {
   float max = 1.0, min = 0.0, test;
   int i;
   int N=100000;

   for (i = 0; i < N; i++)
   {
      float one_plus_test;

      test = (max + min) / ((float)2.0);
      one_plus_test = ((float)1.0) + test;
      if (one_plus_test == ((float)1.0))
      {
         min = test;
      }
      else
      {
         max = test;
      }
   }
   printf("The epsilon machine is %.50lf\n", max);
   return max;
}
