// returns random values to a vector

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>

void setRandom1(gsl_matrix *v,int factor) {

  int i,j;
  //  printf("size1: %d size2 :%d\n",v->size1,v->size2);

  const gsl_rng_type *T;
  gsl_rng *r1;
  gsl_rng_env_setup();
  
  T=gsl_rng_default;
  r1=gsl_rng_alloc(T);

  for (i=0; i<v->size1; i++)
    for(j=0; j<v->size2; j++) 
      gsl_matrix_set(v,i,j,gsl_rng_uniform(r1)*factor);

  gsl_rng_free(r1);

}

int setRandom2(int len) {

  int v;
  float vv;

  const gsl_rng_type *T;
  gsl_rng *r1;
  gsl_rng_env_setup();
  
  T=gsl_rng_default;
  r1=gsl_rng_alloc(T);
  v=gsl_rng_uniform(r1);

  vv=gsl_rng_uniform(r1)*len;
  printf("random: %f\n",vv);
  gsl_rng_free(r1);

  return v;
}
