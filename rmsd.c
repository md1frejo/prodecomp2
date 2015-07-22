/* rmsd function for prodecomp2 */

#include <gsl/gsl_linalg.h>
#include <math.h>

float rmsd1(gsl_matrix *s,gsl_matrix *r,int range,float cutoff) {
  
  int i=0,j=0;
  int dim1=r->size1;
  int dim2=r->size2;
  float n=dim1*dim2;
  float rmsd=0;
  float cut=cutoff*cutoff;

  gsl_matrix *st=gsl_matrix_alloc(dim1,dim2);
  gsl_matrix *rt=gsl_matrix_alloc(dim1,dim2);

  gsl_matrix_memcpy(st,s);
  gsl_matrix_memcpy(rt,r);

  gsl_matrix_sub(st,rt);
  gsl_matrix_memcpy(rt,st);
  gsl_matrix_mul_elements(st,rt);
  //  printf("dim1: %d dim2: %d\n",dim1,dim2);
  //  n=0;
  cut=0;
  for(j=0; j<dim1; j++)
    for(i=0; i<dim2; i++)
      if(gsl_matrix_get(st,j,i)>cut) 
	rmsd+=gsl_matrix_get(st,j,i);    
  
  //  printf("sum: %f n: %f\n sum/n %f\n",sum,n,sum/n);
  rmsd=rmsd/n;
  rmsd=sqrt(rmsd);

  gsl_matrix_free(st);
  gsl_matrix_free(rt);

  return rmsd;
}

float rmsd2(gsl_matrix *s,gsl_matrix *r,int range,int offs) {
  
  int i=0,k=0;
  int dim1=r->size1,dim2=r->size2;
  float n=dim1*range;
  float rmsd=0,rmsd2=0;
  //  printf("dim1: %d n: %f offs: %d range: %d\n",dim1,n,offs,range);
  gsl_matrix *st=gsl_matrix_alloc(dim1,dim2);
  gsl_matrix *rt=gsl_matrix_alloc(dim1,dim2);

  gsl_matrix_memcpy(st,s);
  gsl_matrix_memcpy(rt,r);
  gsl_matrix_sub(st,rt);
  gsl_matrix_memcpy(rt,st);
  gsl_matrix_mul_elements(st,rt);

  //  gsl_matrix_view temp=gsl_matrix_submatrix(s,0,offs*range,dim1,range);
  //  printf("max: %f\n",gsl_matrix_max(&temp.matrix));

  for(i=0; i<dim1; i++)
    for(k=0; k<range; k++) 
      rmsd+=gsl_matrix_get(st,i,k+range*offs);
   
  //  printf("rmsd: %f rmsd/n %f sqrt(rmsd/n): %f\n",rmsd,rmsd/n,sqrt(rmsd/n));

  for(i=0; i<dim1; i++)
    for(k=0; k<range; k++) 
      rmsd2+=gsl_matrix_get(rt,i,k+range*offs);

  //  printf("rmsd2: %f rmsd2/n %f sqrt(rmsd2/n): %f\n",rmsd2,rmsd2/n,sqrt(rmsd2/n));
  rmsd=rmsd/n;
  rmsd=sqrt(rmsd);
  
  gsl_matrix_free(st);
  gsl_matrix_free(rt);

  return rmsd;
}
