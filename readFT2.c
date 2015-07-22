/* reads FT2 into matrixes */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_statistics.h>
#include "headerfiles/common.h"

/* reads ft2 format and returns matrices with stacked data */

int intcompare(const void *a,const void *b) {

  int *x = (int *) a;
  int *y = (int *) b;
  return *x - *y;
}

/* int floatcompare(const void *a,const void *b) { */
  
/*   double *x = (double *) a; */
/*   double *y = (double *) b; */

/*   if(*x<*y)  */
/*     return -1; */
/*   else if(*x>*y)  */
/*     return 1;  */
/*   return 0; */
/* } */

double readFT2(float ***spec2,experiments *exp,int c,int x, int y,int start) {

  int j,k,bytes;
  //  char *line2=NULL;
  int header[2048];
  //  double *data=(double*) malloc(x*y*sizeof(double));
  printf("\nplanes: %d \nindirect points: %d \ndirect points: %d\n",c,x,y);
  FILE *fp;

  for(k=0; k<c; k++) {
    fp=fopen(exp[k].path,"rb");
    printf("reading spectra: %s\n",exp[k].path);
    bytes=0;
    if(fp==NULL)
      printf("could not open file\n");
    bytes+=fread(header,sizeof(int),2048,fp);
    printf("bytes (header): %d\n",bytes);
    bytes=0;
    fseek(fp,2048,SEEK_SET); // offset has to be explained
    for(j=0; j<x; j++) {
      bytes+=fread(spec2[k][j],sizeof(float),y,fp);
    }
    /* printf("spectra %s %d read\n",line2,i); */
    fclose(fp);
  }

  return 0;
}
