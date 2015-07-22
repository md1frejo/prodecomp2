/* implmentation of the midmatrix2.m used in prodecomp 
 *
 * Version 0.02 Jonas Fredriksson 11 06 2013
 * 
 */

#include "headerfiles/common2.h"
#include "headerfiles/midmatrix2.h"

void midmatrix2(gsl_matrix** Fout,gsl_matrix* Fin,int* addorsub) { 

  int i,j,dummy=0;
  int a,b,c,d;
  int m=Fin->size1, p=Fin->size2;

  // Fout=zeros(m,m,p);

/*   for (i=0; i<m; i++)  */
/*     for (j=0; j<m; j++)  */
/*       for (l=0; l<p; l++)  */
/* 	Fout[i][j][l]=0; */
  
  /* Fin(:,1,:)=[Fin((m-1)/2:-1:1,1,:);Fin(m:-1:(m+1)/2,1,:)]; */
  /* m-1 instead of (m/2)-1  */

  for(i=0,j=(m/2)-1; i<m/4; i++,j--) 
    gsl_matrix_swap_rows(Fin,i,j);
 
  for(i=m/2,j=m-1; i<(m-(m/4))-1; i++,j--) 
    gsl_matrix_swap_rows(Fin,i,j);

  // Fin(:,1,:)=[Fin(m,1,:);Fin(1:m-1,1,:)];
  // Fout(i,:,:)=[Fin(:,1,:)];
  /* put elemnt directly into Fout? */

  // alternative solution
  for(d=0; d<p; d++) {
    for(a=0; a<m; a++) {
      for(b=a,c=m-1; b<m; b++,c--)
	gsl_matrix_set(Fout[d],b,a,gsl_matrix_get(Fin,c,d));
      //	Fout[b][a][d]=gsl_matrix_get(Fin,c,d);
      for(b=0; b<a; b++,c--)
	gsl_matrix_set(Fout[d],b,a,gsl_matrix_get(Fin,c,d));
      //	Fout[b][a][d]=gsl_matrix_get(Fin,c,d);
    }
  } 

  if(addorsub[1]==1)
    dummy++;
  else if(addorsub[1]==-1) 
    flip_matrix_col(Fout,m,p);
}
  
/* flipmatrix implementation for 3D */
void flip_matrix_col(gsl_matrix** Fout,int m,int p) {

  int x,y,z;
  float temp;
  /* y<m/2? */
  for(x=0; x<m; x++) {
    for(y=0; y<m/2; y++) {
      for(z=0; z<p; z++) {
	temp=gsl_matrix_get(Fout[z],x,y);
	gsl_matrix_set(Fout[z],x,y,gsl_matrix_get(Fout[z],x,m-y-1));
	gsl_matrix_set(Fout[z],x,m-y-1,temp);
/* 	temp=Fout[x][y][z]; */
/* 	Fout[x][y][z]=Fout[x][m-y-1][z]; */
/* 	Fout[x][m-y-1][z]=temp; */
      }
    }
  }
}
