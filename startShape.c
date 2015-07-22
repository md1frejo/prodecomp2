/* generate start shapes when more than one component */

#include <gsl/gsl_linalg.h>
#include "headerfiles/common.h"

// find max for every peak
void startMax(gsl_vector *pp,gsl_vector* maxp, int comps,int mid,int range) {

  gsl_vector *temp3=gsl_vector_alloc(pp->size);
  int i,h;

  gsl_vector_memcpy(temp3,pp);
 
  for(h=0; h<comps; h++) {
    gsl_vector_set(maxp,h,gsl_vector_max_index(pp));
    if((gsl_vector_get(maxp,h)-range)<0 ||\
       ((gsl_vector_get(maxp,h)+range)>pp->size) ||\
       gsl_vector_max(pp)==0)
      gsl_vector_set(maxp,h,setRandom2(pp->size));
    //    printf("set random: %d %u\n",setRandom2(pp->size),pp->size);
    //    printf("h: %d maX: %f\n",h,gsl_vector_get(maxp,h));
    for(i=gsl_vector_get(maxp,h)-range; i<gsl_vector_get(maxp,h)+range; i++) {
      gsl_vector_set(pp,i,0);
    }
  }
  gsl_vector_memcpy(pp,temp3);
  gsl_vector_free(temp3);
}  

int startShape(gsl_vector *pp,gsl_vector *max,gsl_matrix *stemp,int range,int comps) {

  int h,i;

  for(h=0; h<comps; h++) {
    if((gsl_vector_get(max,h)-range)<0 || ((gsl_vector_get(max,h)+range)>pp->size)) 
      return 0;
    for(i=gsl_vector_get(max,h)-range; i<gsl_vector_get(max,h)+range; i++) {
      gsl_matrix_set(stemp,h,i,gsl_vector_get(pp,i));
      gsl_vector_set(pp,i,0);
    }
  }
  return 1;
}
