#include "headerfiles/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/* read and convert peaklist in xeasy format to internal format */

void ppm2pts1(header *hed,char *path3,int **peaklist,int intervals) {

  int i;
  FILE *f;
  char *t,*token1;
  float dir,ind;

  f=fopen(path3,"r");
  char buffer[200];
  for(i=0; i<3; i++) // hardcoded
    t=fgets(buffer,200,f);
  token1=strtok(buffer," ");

  for(i=0; i<intervals; i++) {
    t=fgets(buffer,200,f);
    if(t==NULL)
      break;
    token1=strtok(buffer," ");
    token1=strtok(NULL," ");
    //    printf("%s %d %d %s\n",path3,intervals,i,token1);
    dir=atof(token1);
    peaklist[i][0]=round((1-(dir-(hed[0].O1_ppm-(hed[0].SW_ppm/2.0)))/hed[0].SW_ppm)*hed[0].size);
    //    printf("%d dir: %d ",i,peaklist[i][0]);
    token1=strtok(NULL," ");
    ind=atof(token1);
    peaklist[i][1]=round((1-(ind-(hed[1].O1_ppm-(hed[1].SW_ppm/2.0)))/hed[1].SW_ppm)*hed[1].size);
    //    printf("ind: %d\n",peaklist[i][1]);
    t=t;
    peaklist[i][0]-=1;
    peaklist[i][1]-=1;
  }
  fclose(f);
}

int getComps(int** peaklist,int x1,int x2,int r,int intervals) {

  int i,comps=0;
  //  printf("%d %d %d %d\n",x1,x2,r,intervals);
  for(i=0; i<intervals; i++) {
    if(x1-r>=0 && peaklist[i][0]>=x1-r && peaklist[i][0]<=x2+r) { // note partial check
      printf("\tfound interval: %d dir: %d ind: %d\n",i,peaklist[i][0],peaklist[i][1]);
      peaklist[i][2]=1;
      //      printf("read: %d %d\n",peaklist[i][0],peaklist[i][1]);
      comps++;
    }
  }
  return comps;
}
  

int getPeaks(int** peaklist,int **peaklist2,int x1,int x2,int intervals,int comps,\
	     int inter,int random,int r) {

  int i,j=999; // hardcoded
  //  printf("x1: %d x2: %d\n",x1,x2);
  if(random==0 && comps>1) {
    for(j=1,i=0; i<intervals; i++) {
      if(peaklist2[0][0]!=peaklist[i][0] && peaklist2[0][1]!=peaklist[i][1]\
	 && peaklist[i][0]>=x1-r && peaklist[i][0]<=x2+r && j<comps) {
	peaklist2[j][0]=peaklist[i][0];
	peaklist2[j][1]=peaklist[i][1];
	j++;
      }
    }
    if(j>=comps)
      for(i=j; i<comps; i++) 
	peaklist2[i][0]=peaklist2[i][1]=-1;
  }
  else {
    for(i=1; i<comps; i++) 
      peaklist2[i][0]=peaklist2[i][1]=-1;
  }
  return j;
}

int getPeaks2(int** peaklist,int x1,int x2,int intervals,int comps,\
	     int inter,int random,int r) {

  int i,j=999; // hardcoded
  int nh=peaklist[inter][0];
  int n=peaklist[inter][1];

  peaklist[inter][2]=1;

  if(random==0 && comps>1) {
    for(j=0,i=0; i<intervals; i++) {
      if(nh!=peaklist[i][0] && n!=peaklist[i][1]\
	 && peaklist[i][0]>=x1-r && peaklist[i][0]<=x2+r) {
	peaklist[i][2]=1;
	//	printf("read: %d %d\n",peaklist[i][0],peaklist[i][1]);
	j++;
      }
    }
  }
  return j;
}

void readNMRDraw(int** peaklist,char *path4,int itervals) {
  /* read peaklist from nmrDraw */
  ;
}
