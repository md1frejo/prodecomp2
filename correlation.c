/* correlation function for prodecomp2 */

#include "headerfiles/common.h"
#include "headerfiles/common2.h"

void correlation(gsl_matrix *chains,double ***finalShapes,gsl_matrix *corrMatrix,int ind,seqs *seq,float cut) {

  int d=0,f=0,g,h,i,j,pos1=0,pos2=0;
  int size=corrMatrix->size1;
  float corr,max1,max2;
 
  double *hx=malloc(ind*sizeof(double));
  double *cx=malloc(ind*sizeof(double));
 
  for(i=0; i<size; i++) {
    for(h=0; h<ind; h++) {
      cx[h]=finalShapes[i][4][h]; // hardcoded
      hx[h]=finalShapes[i][6][h];
    }

    for(h=0; h<ind; h++) {
      cx[h]+=finalShapes[i][5][h];
      hx[h]+=finalShapes[i][7][h];
    }

    for(j=0; j<size; j++) {
      corr=0;
      corr=cosinesim(cx,finalShapes[j][2],ind);
      corr+=cosinesim(hx,finalShapes[j][3],ind);
      gsl_matrix_set(corrMatrix,i,j,corr/2.0);
    }
  }

  // set all diagonal elements to zero
  for(i=0; i<size; i++)
    gsl_matrix_set(corrMatrix,i,i,0);

  // remove negative correlation
  for(i=0; i<size; i++)
    for(j=0; j<size; j++)
      if(gsl_matrix_get(corrMatrix,i,j)<cut)
	gsl_matrix_set(corrMatrix,i,j,0);

  // remove circular correlation max 2
  for(i=0; i<size; i++)
    for(j=0; j<size; j++) {
      if(gsl_matrix_get(corrMatrix,i,j)<gsl_matrix_get(corrMatrix,j,i))
	gsl_matrix_set(corrMatrix,i,j,0);
      else
	gsl_matrix_set(corrMatrix,j,i,0);
    }
  
  // find max correlation and set the rest to zero
  for(h=0,i=0; i<size; i++) {
    max1=0;
    max2=0;
    for(j=0; j<size; j++) {
      if(gsl_matrix_get(corrMatrix,i,j)>max1) {
	max1=gsl_matrix_get(corrMatrix,i,j);
	pos1=j;
      }
    }
    for(j=0; j<size; j++) {
      if(gsl_matrix_get(corrMatrix,j,pos1)>max2) {
	max2=gsl_matrix_get(corrMatrix,j,pos1);
	pos2=j;
      }
    }
    //    printf("%d %d %d\n",i,pos1,pos2);
    if(pos2==i) {
      //      gsl_matrix_set(corrMatrix,pos1,pos2,0);
      for(j=0; j<size; j++) 
	if(j!=pos1)
	  gsl_matrix_set(corrMatrix,pos2,j,0);
      for(j=0; j<size; j++) 
	if(j!=pos2)
	  gsl_matrix_set(corrMatrix,j,pos1,0);
    }
  }
  
  gsl_vector *spoints=gsl_vector_alloc(size);
  gsl_vector_set_all(spoints,-1);

  for(g=0,i=0; i<size; i++) {
    for(j=0; j<size; j++) 
      if(gsl_matrix_get(corrMatrix,i,j)!=0)
	break;
    if(j==size) {
      h=startingPoint(corrMatrix,i);
      if(h!=-1)
	gsl_vector_set(spoints,g++,startingPoint(corrMatrix,i));
    }
  }	 
  //  f=startingPoint(corrMatrix,0);
  //  gsl_vector_fprintf(stdout,spoints,"%f");
  
  saveMatrix(corrMatrix,"corr.csv",1);
  g=0;
  h=0;
  f=gsl_vector_get(spoints,d++);
  while(1) {
    for(i=0; i<size; i++) 
      if(gsl_matrix_get(corrMatrix,f,i)!=0) {
	//	f=gsl_vector_get(spoints,d++);
	gsl_matrix_set(chains,h,g++,f);
	//	gsl_matrix_set(corrMatrix,f,i,0);
	f=i;
	break;
      }
    if(i==size) {
      //f=startingPoint(corrMatrix,0);
      f=gsl_vector_get(spoints,d++);
      g=0;
      h++;
      i=0;
    }
    if(f==-1)
      break;
  }
  
  saveMatrix(chains,"chains.csv",1);
  saveMatrix(chains,"chains.txt",0);
  
  gsl_vector_free(spoints);
  
  free(hx);
  free(cx);
}
