/* get different parameters */

#include "headerfiles/common.h"
//#include "main.h"

//#include "getParam.h"

void getDefs(definitions *defs) {

  /* temporary */
  defs[0].sh="Hab";
  defs[0].conv=0;
  defs[0].pos=calloc(5,sizeof(int));
  defs[0].pos[0]=3;
  defs[0].sign=calloc(5,sizeof(int));
  defs[0].sign[0]=1;

  defs[1].sh="Cab";
  defs[1].conv=0;
  defs[1].pos=calloc(5,sizeof(int));
  defs[1].pos[0]=2;
  defs[1].sign=calloc(5,sizeof(int));
  defs[1].sign[0]=1;

  defs[2].sh="Cab+Hab";
  defs[2].conv=2;
  defs[2].pos=calloc(5,sizeof(int));
  defs[2].pos[0]=2;
  defs[2].pos[1]=3;
  defs[2].sign=calloc(5,sizeof(int));
  defs[2].sign[0]=1;
  defs[2].sign[1]=1;

  defs[3].sh="Cab-Hab";
  defs[3].conv=2;
  defs[3].pos=calloc(5,sizeof(int));
  defs[3].pos[0]=2;
  defs[3].pos[1]=3;
  defs[3].sign=calloc(5,sizeof(int));
  defs[3].sign[0]=1;
  defs[3].sign[1]=-1;

  defs[4].sh="CO";
  defs[4].conv=0;
  defs[4].pos=calloc(5,sizeof(int));
  defs[4].pos[0]=1;
  defs[4].sign=calloc(5,sizeof(int));
  defs[4].sign[0]=1;

  defs[5].sh="CO+Hab";
  defs[5].conv=2;
  defs[5].pos=calloc(5,sizeof(int));
  defs[5].pos[0]=1;
  defs[5].pos[1]=3;
  defs[5].sign=calloc(5,sizeof(int));
  defs[5].sign[0]=1;
  defs[5].sign[1]=1;

  defs[6].sh="CO-Hab";
  defs[6].conv=2;
  defs[6].pos=calloc(5,sizeof(int));
  defs[6].pos[0]=1;
  defs[6].pos[1]=3;
  defs[6].sign=calloc(5,sizeof(int));
  defs[6].sign[0]=1;
  defs[6].sign[1]=-1;

  defs[7].sh="CO+Cab";
  defs[7].conv=2;
  defs[7].pos=calloc(5,sizeof(int));
  defs[7].pos[0]=1;
  defs[7].pos[1]=2;
  defs[7].sign=calloc(5,sizeof(int));
  defs[7].sign[0]=1;
  defs[7].sign[1]=1;

  defs[8].sh="CO-Cab";
  defs[8].conv=2;
  defs[8].pos=calloc(5,sizeof(int));
  defs[8].pos[0]=1;
  defs[8].pos[1]=2;
  defs[8].sign=calloc(5,sizeof(int));
  defs[8].sign[0]=1;
  defs[8].sign[1]=-1;

  defs[9].conv=2;
  defs[9].sign=calloc(5,sizeof(int));
  defs[9].sign[0]=1;

  defs[10].conv=2;
  defs[10].sign=calloc(5,sizeof(int));
  defs[10].sign[0]=-1;
  
  defs[13].sh="N";
  defs[13].conv=0;
  defs[13].pos=calloc(5,sizeof(int));
  defs[13].pos[0]=0;
  defs[13].sign=calloc(5,sizeof(int));
  defs[13].sign[0]=1;

  defs[22].sh="N+CO";
  defs[22].pos=calloc(5,sizeof(int));
  defs[22].pos[0]=0;
  defs[22].pos[1]=1;
  defs[22].conv=2;
  defs[22].sign=calloc(5,sizeof(int));
  defs[22].sign[0]=1;
  defs[22].sign[1]=1;

  defs[23].sh="N-CO";
  defs[23].pos=calloc(5,sizeof(int));
  defs[23].pos[0]=0;
  defs[23].pos[1]=1;
  defs[23].conv=2;
  defs[23].sign=calloc(5,sizeof(int));
  defs[23].sign[0]=1;
  defs[23].sign[1]=-1;

  defs[14].sh="N+Hab";
  defs[14].pos=calloc(5,sizeof(int));
  defs[14].pos[0]=0;
  defs[14].pos[1]=3;
  defs[14].conv=2;
  defs[14].sign=calloc(5,sizeof(int));
  defs[14].sign[0]=1;
  defs[14].sign[1]=1;

  defs[15].sh="N-Hab";
  defs[15].pos=calloc(5,sizeof(int));
  defs[15].pos[0]=0;
  defs[15].pos[1]=3;
  defs[15].conv=2;
  defs[15].sign=calloc(5,sizeof(int));
  defs[15].sign[0]=1;
  defs[15].sign[1]=-1;

  defs[16].sh="N+Cab";
  defs[16].pos=calloc(5,sizeof(int));
  defs[16].pos[0]=0;
  defs[16].pos[1]=2;
  defs[16].conv=2;
  defs[16].sign=calloc(5,sizeof(int));
  defs[16].sign[0]=1;
  defs[16].sign[1]=1;

  defs[17].sh="N-Cab";
  defs[17].pos=calloc(5,sizeof(int));
  defs[17].pos[0]=0;
  defs[17].pos[1]=2;
  defs[17].conv=2;
  defs[17].sign=calloc(5,sizeof(int));
  defs[17].sign[0]=1;
  defs[17].sign[1]=-1;

  defs[24].sh="N+CO+Hab";
  defs[24].pos=calloc(5,sizeof(int));
  defs[24].pos[0]=0;
  defs[24].pos[1]=1;
  defs[24].pos[2]=3;
  defs[24].conv=3;
  defs[24].sign=calloc(5,sizeof(int));
  defs[24].sign[0]=1;
  defs[24].sign[1]=1;
  defs[24].sign[2]=1;

  defs[25].sh="N-CO+Hab";
  defs[25].pos=calloc(5,sizeof(int));
  defs[25].pos[0]=0;
  defs[25].pos[1]=1;
  defs[25].pos[2]=3;
  defs[25].conv=3;
  defs[25].sign=calloc(5,sizeof(int));
  defs[25].sign[0]=1;
  defs[25].sign[1]=-1;
  defs[25].sign[2]=1;

  defs[26].sh="N+CO-Hab";
  defs[26].pos=calloc(5,sizeof(int));
  defs[26].pos[0]=0;
  defs[26].pos[1]=1;
  defs[26].pos[2]=3;
  defs[26].conv=3;
  defs[26].sign=calloc(5,sizeof(int));
  defs[26].sign[0]=1;
  defs[26].sign[1]=1;
  defs[26].sign[2]=-1;

  defs[27].sh="N-CO-Hab";
  defs[27].pos=calloc(5,sizeof(int));
  defs[27].pos[0]=0;
  defs[27].pos[1]=1;
  defs[27].pos[2]=3;
  defs[27].conv=3;
  defs[27].sign=calloc(5,sizeof(int));
  defs[27].sign[0]=1;
  defs[27].sign[1]=-1;
  defs[27].sign[2]=-1;

  defs[28].sh="N+CO+Cab";
  defs[28].pos=calloc(5,sizeof(int));
  defs[28].pos[0]=0;
  defs[28].pos[1]=1;
  defs[28].pos[2]=2;
  defs[28].conv=3;
  defs[28].sign=calloc(5,sizeof(int));
  defs[28].sign[0]=1;
  defs[28].sign[1]=1;
  defs[28].sign[2]=1;

  defs[29].sh="N-CO+Hab";
  defs[29].pos=calloc(5,sizeof(int));
  defs[29].pos[0]=0;
  defs[29].pos[1]=1;
  defs[29].pos[2]=3;
  defs[29].conv=3;
  defs[29].sign=calloc(5,sizeof(int));
  defs[29].sign[0]=1;
  defs[29].sign[1]=-1;
  defs[29].sign[2]=1;

  defs[30].sh="N+CO-Hab";
  defs[30].pos=calloc(5,sizeof(int));
  defs[30].pos[0]=0;
  defs[30].pos[1]=1;
  defs[30].pos[2]=3;
  defs[30].conv=3;
  defs[30].sign=calloc(5,sizeof(int));
  defs[30].sign[0]=1;
  defs[30].sign[1]=1;
  defs[30].sign[2]=-1;

  defs[31].sh="N-CO-Hab";
  defs[31].pos=calloc(5,sizeof(int));
  defs[31].pos[0]=0;
  defs[31].pos[1]=1;
  defs[31].pos[2]=3;
  defs[31].conv=3;
  defs[31].sign=calloc(5,sizeof(int));
  defs[31].sign[0]=1;
  defs[31].sign[1]=-1;
  defs[31].sign[2]=-1;
}

void getInterv(interval *interv) {

  interv[0].x1=437;
  interv[0].x2=439;
  interv[0].y=256;
  interv[0].comps=1;
  interv[0].conv=0;
  
  interv[1].x1=59;
  interv[1].x2=61;
  interv[1].y=256;
  interv[1].comps=1;
  interv[1].conv=0;

  interv[2].x1=180;
  interv[2].x2=182;
  interv[2].y=256;
  interv[2].comps=2;
  interv[2].conv=3;

  interv[3].x1=250;
  interv[3].x2=252;
  interv[3].y=256;
  interv[3].comps=5;
  interv[3].conv=3;

  interv[4].x1=338;
  interv[4].x2=341;
  interv[4].y=256;
  interv[4].comps=3;
  interv[4].conv=1;

  interv[5].x1=234;
  interv[5].x2=236;
  interv[5].y=256;
  interv[5].comps=1;
  interv[5].conv=0;

  interv[6].x1=270;
  interv[6].x2=272;
  interv[6].y=256;
  interv[6].comps=2;
  interv[6].conv=1;

  interv[7].x1=174;
  interv[7].x2=176;
  interv[7].y=256;
  interv[7].comps=1;
  interv[7].conv=0;

  interv[8].x1=294;
  interv[8].x2=296;
  interv[8].y=256;
  interv[8].comps=2;
  interv[8].conv=0;
}

void getProj(gsl_matrix *proj) {

  gsl_matrix_set(proj,0,0,0); // Hab
  gsl_matrix_set(proj,1,0,0); // Cab  
  gsl_matrix_set(proj,2,0,1); // Cab+Hab
  gsl_matrix_set(proj,3,0,-1); // Cab-Hab
  gsl_matrix_set(proj,4,0,0); // CO
  gsl_matrix_set(proj,5,0,0); // CO+Hab
  gsl_matrix_set(proj,6,0,0); // CO-Hab
  gsl_matrix_set(proj,7,0,0); // CO+Cab
  gsl_matrix_set(proj,8,0,0); // CO-Cab
  gsl_matrix_set(proj,9,0,0); // CO+Cab+Hab
  gsl_matrix_set(proj,10,0,0); // CO-Cab+Hab
  gsl_matrix_set(proj,11,0,0); // CO+Cab-Hab
  gsl_matrix_set(proj,12,0,0); // CO-Cab-Hab
  gsl_matrix_set(proj,13,0,0); // N
  gsl_matrix_set(proj,14,0,0); // N+Hab
  gsl_matrix_set(proj,15,0,0); // N-Hab
  gsl_matrix_set(proj,16,0,0); // N+Cab
  gsl_matrix_set(proj,17,0,0); // N-Cab
  gsl_matrix_set(proj,18,0,0); // N+Cab+Hab
  gsl_matrix_set(proj,19,0,0); // N-Cab+Hab
  gsl_matrix_set(proj,20,0,0); // N+Cab-Hab
  gsl_matrix_set(proj,21,0,0); // N-Cab-Hab
  gsl_matrix_set(proj,22,0,0); // N+CO
  gsl_matrix_set(proj,23,0,0); // N-CO
  gsl_matrix_set(proj,24,0,0); // N+CO+Hab
  gsl_matrix_set(proj,25,0,0); // N-CO+Hab
  gsl_matrix_set(proj,26,0,0); // N+CO-Hab
  gsl_matrix_set(proj,27,0,0); // N-CO-Hab
  gsl_matrix_set(proj,28,0,0); // N+CO+Cab
  gsl_matrix_set(proj,29,0,0); // N-CO+Cab
  gsl_matrix_set(proj,30,0,0); // N+CO-Cab
  gsl_matrix_set(proj,31,0,0); // N-CO-Cab

}

void getSet(int **set) {

  set[0][0]=13; // N
  set[0][1]=22; // N+CO
  set[0][2]=23; // N-CO
  set[0][3]=16; // N+Cab
  set[0][4]=17; // N-Cab
  set[0][5]=14; // N+Hab
  set[0][6]=15; // N-Hab
  set[0][7]=24; // N+CO+Hab
  set[0][8]=25; // N-CO+Hab
  set[0][9]=26; // N+CO-Hab
  set[0][10]=27; // N-CO-Hab
  
  set[1][0]=4; // CO
  set[1][1]=22; // N+CO
  set[1][2]=23; // N-CO
  set[1][3]=7; // CO+Cab
  set[1][4]=8; // CO-Cab
  set[1][5]=5; // CO+Hab
  set[1][6]=6; // CO-Hab
  set[1][7]=24; // N+CO+Hab
  set[1][8]=25; // N-CO+Hab
  set[1][9]=26; // N+CO-Hab
  set[1][10]=27; // N-CO-Hab

  set[2][0]=1; // Cab
  set[2][1]=16; // Cab+N
  set[2][2]=17; // Cab-N
  set[2][3]=7; // Cab+CO
  set[2][4]=8; // Cab-CO
  set[2][5]=2; // Cab+Hab
  set[2][6]=3; // Cab-Hab
  set[2][7]=28; // N+CO+Cab 
  set[2][8]=29; // N-CO+Cab 
  set[2][9]=30; // N+CO-Cab 
  set[2][10]=31; // N-CO-Cab

  set[3][0]=0; // Hab
  set[3][1]=14; // N+Hab
  set[3][2]=15; // N-Hab
  set[3][3]=2; // Hab+Cab
  set[3][4]=3; // Hab-Cab
  set[3][5]=5; // CO+Hab
  set[3][6]=6; // CO-Hab
  set[3][7]=24; // N+CO+Hab
  set[3][8]=25; // N-CO+Hab
  set[3][9]=26; // N+CO-Hab
  set[3][10]=27; // N-CO-Hab
}

void getMatrix(gsl_matrix *cm,experiments exp,int shape,int shapes) {

  int i;

  /* printf("getMatrix: %s \n",exp.sh); */
  for(i=0; i<exp.conv; i++)  
    gsl_matrix_set(cm,exp.pos[i],shape,(int) exp.sign[i]);
  /* if(shape==1) */
  /*   for(i=0; i<8; i++) */
  /*     //      if(gsl_matrix_get(cm,i,shape)) */
  /*     printf("%f\n",gsl_matrix_get(cm,i,shape)); */
}
