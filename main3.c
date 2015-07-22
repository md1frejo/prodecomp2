/* main file for nnls */

/* mac: gcc -g -m32 -Wall -O3 -lgsl -lgslcblas -lm main2.c readFT2.c decompose.c fastnnls.c getParam.c util.c midmatrix3.c startShape.c -o main2 (not complete) */

/* ubuntu:  gcc -g -pg -O3 -Wall -pedantic -std=c99 readFT2.c decompose2.c fastnnls.c midmatrix3.c setRandom.c startShape.c getParam.c rmsd.c util.c readInput.c readPeaklist.c util2.c correlation.c -fopenmp main2.c -lgsl -lcblas -latlas -lm -o main */

/* gcc -g -pg -Wall -pedantic -std=c99 -Ofast readFT2.c decompose2.c fastnnls.c midmatrix3.c setRandom.c startShape.c getParam.c rmsd.c readInput.c readPeaklist.c util2.c correlation.c -fopenmp main2.c -L/usr/lib/openblas-base -lgsl -lopenblas -lm -o main */

/* export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/openblas-base */

/* ex backbone: */
/* ./main b prodecomp.txt hsqc.peaks 2 */

/* ex noesy: */
/* ./main n prodecompNnoesy.txt hsqcNnoesy.tab 11 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.h"
#include "rmsd.h"
#include "main.h"
#include "headerfiles/common.h"
#include "headerfiles/fastnnls.h"

#include <gsl/gsl_blas.h>

int main(int argc, char **argv) {

  srand((unsigned)time(NULL));

  char *p,*path1=NULL,*path3=NULL,*residue,*etype,*rdir,*ldir,*token,*pname;
  time_t t0,t1,tot0=0,tot1,current_time;
  clock_t c0,c1,ctot0=0,ctot1=0;
  //  char* c_time_string;
  (void) residue;
  //  int *pcount=NULL;
  int iterations=10,sel=-1,maxcomps=10,state=0,format=0,nch=0,selected=0;
  int g,h,i,j,k,l,r;
  int x1=0,x2=0,shape,comp,comps=1,exps=0,intervals,sum,sub,all=0,len1,len2,pc=0,orgcomps=0,br=0;
  //  int nthreads, tid;
  float rmsd=0,rmsdprev=0,tres=0,cut=0.010,cut2=0.10,cut3=0.2,max; 
  // cut: for rmsd. cut2: for n check, cut3: for peakpicking 
  //  int expected[]={1,1,3,3,2,2,2,3};
  int compare[]={2,4,5};
  int minus[]={0,1,1,1,0,0,0,0};
  // Hab: 0, Cab: 1, CO: 4, N: 13, N+CO: 22 N-CO:23)
  // struct for all intervals
  FILE *fs,*peakl,*flog,*rlist,*latestf;
  int x,y;

  //  (void) pname; // get rid of gcc compiler warnings

  time_t t = time(NULL);
  struct tm tm = *localtime(&t);
  char *tbuff=calloc(100,sizeof(char));
  char *lbuff=calloc(100,sizeof(char));
  char *pbuff=calloc(100,sizeof(char));
  
  printf("now: %d-%d-%d %d:%d:%d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

  if(argc<4 || argc==1) {
    printf("not enough arguments\nex: ./main b prodecomp.txt hsqc.peaks\n");
    return 0;
  }  
  else {
    etype=argv[1];
    path1=argv[2];
    path3=argv[3];
  }

  gsl_vector *order=gsl_vector_calloc(2);

  format=checkPeakFormat(path3,order);

  if(format==-10) {
    printf("wrong peak format\n");
    return 1;
  }
  else
    printf("peak format: %d\n",format);

  int nshapes=getnshapes(path1);

  int *setup=(int *) malloc(nshapes*sizeof(int));
  pname=NULL;
  ldir=NULL;
  mkdir("results",0776);
  if(strcmp(etype,"b")==0) {
    sprintf(tbuff,"%s-%d-%d-%d_%d_%d","results/backboneResults",tm.tm_year + 1900,\
	    tm.tm_mon + 1, tm.tm_mday, tm.tm_hour,tm.tm_min);
    sprintf(lbuff,"%s-%d-%d-%d_%d_%d","results/backboneLog.txt",tm.tm_year + 1900,\
	    tm.tm_mon + 1, tm.tm_mday, tm.tm_hour,tm.tm_min);
    sprintf(pbuff,"%s-%d-%d-%d_%d_%d","results/backbonePeakList.txt",tm.tm_year + 1900,\
	    tm.tm_mon + 1, tm.tm_mday, tm.tm_hour,tm.tm_min);
    mkdir(tbuff,0776);
    ldir=lbuff;
    rdir=tbuff;
    pname=pbuff;
  }
  if(strcmp(etype,"t")==0) {
    sprintf(tbuff,"%s-%d-%d-%d_%d_%d","results/tocsyResults",tm.tm_year + 1900,\
	    tm.tm_mon + 1, tm.tm_mday, tm.tm_hour,tm.tm_min);
    sprintf(lbuff,"%s-%d-%d-%d_%d_%d","results/tocsyLog.txt",tm.tm_year + 1900,\
	    tm.tm_mon + 1, tm.tm_mday, tm.tm_hour,tm.tm_min);
    sprintf(pbuff,"%s-%d-%d-%d_%d_%d","results/tocsyPeakList.txt",tm.tm_year + 1900,\
	    tm.tm_mon + 1, tm.tm_mday, tm.tm_hour,tm.tm_min);
    mkdir(tbuff,0776);
    ldir=lbuff;
    rdir=tbuff;
    pname=pbuff;
  }
  if(strcmp(etype,"n")==0) {
    sprintf(tbuff,"%s-%d-%d-%d_%d_%d","results/noesyResults",tm.tm_year + 1900,\
	    tm.tm_mon + 1, tm.tm_mday, tm.tm_hour,tm.tm_min);
    sprintf(lbuff,"%s-%d-%d-%d_%d_%d","results/noesyLog.txt",tm.tm_year + 1900,\
	    tm.tm_mon + 1, tm.tm_mday, tm.tm_hour,tm.tm_min);
    sprintf(pbuff,"%s-%d-%d-%d_%d_%d","results/noesyPeakList.txt",tm.tm_year + 1900,\
	    tm.tm_mon + 1, tm.tm_mday, tm.tm_hour,tm.tm_min);
    mkdir(tbuff,0776);
    ldir=lbuff;
    rdir=tbuff;
    pname=pbuff;
  }

  flog=fopen(ldir,"w");
  peakl=fopen(pname,"w");
  rlist=fopen("rlist.txt","w");
  
  /* else */
  /*   printf("experiment can either be backbone (b), tocsy (t) or noesy (n)\n"); */

  char line[100];
  current_time = time(NULL);
  fprintf(flog,"\tlog file prodecmp2: %s\n",ctime(&current_time));
  fprintf(flog,"\trmsd cut: %0.4f\n",cut);
  fprintf(flog,"\tcut for n check: %2.2f\n\n",cut2*100);

  printf("Reading %s %s\n",path1,path3);
  exps=getnexperiments(path1);
  intervals=getIntervals2(path3);

  char *buffer1=calloc(100,sizeof(char*));
  char *buffer2=calloc(100,sizeof(char*));

  char **nuclei=(char **) malloc(nshapes*sizeof(char*));
  for(l=0; l<nshapes; l++)
    nuclei[l]=(char *) malloc(20*sizeof(char));

  printf("Prodecomp2:\nOptions: b backbone t tocsy n noesy p prodecomp.txt d input.txt i intervall\n");

  /* printf("Prodecomp2:\nOptions: input (i) Decompose (d) plot (p) correlate (c)\n\ */
  /*         read FT2 (r) slide (s) quit (q)\n"); */

  char filename1[100];
  char filename2[100];

  interval *interv=malloc(intervals*sizeof(interval));
  definitions *defs=malloc(exps*sizeof(definitions));
  experiments *exp=malloc(exps*sizeof(experiments));

  for(l=0; l<exps; l++) {
    exp[l].path=(char *) malloc(256*sizeof(char));
    exp[l].defs=(int *) calloc(nshapes,sizeof(int));
    exp[l].sh=(char *) malloc(25*sizeof(char)); // tilt have longer names
    exp[l].pos=(int *) calloc(5,sizeof(int));
    exp[l].sign=(int *) calloc(5,sizeof(int));
  }

  header *hed=malloc(nshapes*sizeof(header));
  for(l=0; l<nshapes; l++) {
    hed[l].format=(char *) malloc(20*sizeof(char));
    hed[l].nuc=(char *) malloc(20*sizeof(char));
  }

  readHeader(path1,hed,0,1);

  for(l=0; l<nshapes; l++) {
    if(l==0)
      readHeader(path1,hed,l,1);
    else
      readHeader(path1,hed,l,0);
    strcpy(nuclei[l],hed[l].nuc);
  }

  x=hed[1].size;
  y=hed[0].size;

  if(argc==5)
    sel=atoi(argv[4]);
  if(argc==6) { 
    sel=atoi(argv[4]);
    maxcomps=atoi(argv[5]);
  }
  if(argc==7) {
    sel=atoi(argv[4]);
    maxcomps=atoi(argv[5]);
    br=atoi(argv[6]);
  }
  if(argc>7) {
    x1=atoi(argv[6]);
    x2=atoi(argv[7]);
  }

  double ***finalShapes=(double ***) malloc(intervals*sizeof(double **));

  for(l=0; l<intervals; l++) { 
    finalShapes[l]=(double **) malloc((nshapes-1)*sizeof(double*));
    for(k=0; k<nshapes-1; k++) 
      finalShapes[l][k]=(double *) malloc(x*sizeof(double));
  }

  int **peaklist=(int **) malloc(intervals*sizeof(int *)); 
  for(i=0; i<intervals; i++)
    peaklist[i]=(int *) malloc(3*sizeof(int));

  //  printf("intervalls: %d\n",makeIntervalls1(peaklist,path3));

  getIntervals3(path3,peaklist,interv,intervals,hed,format,flog,order);
  //  ppm2pts1(hed,path3,peaklist,intervals);

  int **tpeaklist=(int **) calloc(1000,sizeof(int *)); // hardcoded
  for(i=0; i<1000; i++)
    tpeaklist[i]=(int *) calloc(2,sizeof(int));

  //  path3="tocsyAll.peaks";

  //  ppm2pts1(hed,path3,tpeaklist,1000);  // hardcoded 

  //  getIntervals3(path3,tpeaklist,interv,intervals,1);

  len1=getFileLength("bmrb.txt");
  len2=getFileLength("hbd6.seq");

  stats *stat=malloc(len1*sizeof(stats)); 

  for(i=0; i<len1; i++) {
    stat[i].res=(char *) malloc(20*sizeof(char));
    stat[i].nucs=(char *) malloc(20*sizeof(char));
  }

  seqs *seq=malloc(len2*sizeof(seqs));
  
  component *components=malloc(intervals*sizeof(component));

  for(i=0; i<intervals; i++) {
    components[i].nuclei=nuclei;
    components[i].shifts=(float **) malloc(nshapes*sizeof(float *));
    components[i].points=(int **) malloc(nshapes*sizeof(int *));
    components[i].intens=(float **) malloc(nshapes*sizeof(float *));
    components[i].compare=(int *) malloc(3*sizeof(int)); //hardcoded
    components[i].compare=compare;
    for(j=0; j<nshapes; j++) {
      components[i].shifts[j]=(float *) calloc(x,sizeof(float));
      components[i].points[j]=(int *) calloc(x,sizeof(int));
      components[i].intens[j]=(float *) calloc(x,sizeof(float));
    }
  }

  gsl_matrix *corrmat=gsl_matrix_alloc(intervals,intervals);
  gsl_matrix *chains=gsl_matrix_alloc(intervals,intervals);
  gsl_matrix *corrMatrix=gsl_matrix_calloc(intervals,intervals);

  gsl_matrix_set_all(chains,-1);

  for(i=0; i<len2; i++) 
    seq[i].residue=(char *) malloc(20*sizeof(char));
    
  peaklists *fpeaklist=malloc(1000*sizeof(peaklists)); // hardcoded

  for(i=0; i<1000; i++)  // hardcoded
    fpeaklist[i].nuc=(char *) malloc(20*sizeof(char));
 
  //  pcount=&count;
  fpeaklist[0].counter=0;

  /* for(i=0; i<20; i++) { */
  /*   stat[i].bmrb=(nstat *) malloc(20*sizeof(nstat)); */
  /*   for(j=0; j<20; j++) */
  /*     stat[i].bmrb[j].nucs=(char *) malloc(20*sizeof(char)); */
  /* } */

  initStats(stat,"bmrb.txt"); // both hardcoded
  getSequence(seq,"hbd6.seq");
  
  float ***spec=(float ***) malloc(exps*sizeof(float **));
  for(i=0; i<exps; i++) {
    spec[i]=(float **) malloc(x*sizeof(float *));
    for(j=0; j<x; j++)
      spec[i][j]=(float *) malloc(y*sizeof(float));
  }
  
  gsl_vector *nshape=gsl_vector_alloc(x);

  readProdecomptext(path1,exp,nshapes,nuclei,exps);
  readFT2(spec,exp,exps,x,y,0);

  //  tiltCheck(exp,hed,nuclei,exps);
  gsl_vector *cutoffs=gsl_vector_alloc(exps);

  // normalize spectra
  //  if(iteration==0 && shape==0) {
  for(max=0,k=0; k<exps; k++) 
    for(j=0; j<y; j++)
      for(l=0; l<x; l++)
	if(max<spec[k][l][j])
	  max=spec[k][l][j];
  
  for(k=0; k<exps; k++) 
    for(j=0; j<y; j++)
      for(l=0; l<x; l++)
	spec[k][l][j]/=max;
  
  // calculate cutoffs and remove noise
  for(k=0; k<exps; k++)
    gsl_vector_set(cutoffs,k,setSpectra(spec[k],x,0,y-1));
  //  }

  //centerShift(spec,x,y,exps,x/2);

  all=nshapes-1;
  int **set=(int **) calloc(all,sizeof(int *));
  for(i=0; i<all; i++) 
   set[i]=(int *) calloc(exps,sizeof(int));
  for(k=0; k<all; k++)
    for(h=0,l=0; l<exps; l++)
      if(exp[l].defs[k+1])
	set[k][h++]=l;

  gsl_vector *temp1s=gsl_vector_alloc(x);
  gsl_vector *temp3s=gsl_vector_alloc(x);
  // setting up kf
  gsl_vector *kf=gsl_vector_alloc(exps);
  gsl_vector_set_all(kf,1.0);

  for(k=0; k<exps; k++) {
    sum=0;
    h=0;
    for(l=1; l<nshapes; l++)
      if(exp[k].defs[l]) {
	exp[k].pos[h]=l-1;
	exp[k].sign[h]=exp[k].defs[l];
	sum++;
	h++;
      }
    if(sum>1)
      exp[k].conv=sum;
    else
      exp[k].conv=0;
  }

  for(k=0; k<exps; k++) {
    printf("exp: %d %s %d ",k,exp[k].sh,exp[k].conv);
    printf(":");
    for(h=0; h<5; h++)
      printf("%d ",exp[k].pos[h]);
    printf(":");
    for(h=0; h<5; h++)
      printf("%d ",exp[k].sign[h]);
    printf(":");
    for(l=0; l<nshapes; l++)
      printf("%d ",exp[k].defs[l]);
    printf("\n");
  }

  while(1) {
    printf("> ");
    p=fgets(line,intervals,stdin);
    token=strtok(line," ");
    if(strcmp(token,"i")==0) {
      readHeader(path1,hed,0,1);
      for(l=0; l<nshapes; l++) {
  	if(l==0)
  	  readHeader(path1,hed,0,1);
  	else
  	  readHeader(path1,hed,l,0);
  	strcpy(nuclei[l],hed[l].nuc);
      }
      readProdecomptext(path1,exp,nshapes,nuclei,exps); // hardcoded
      for(l=0; l<exps; l++) {
      	printf("exp: %d %s",l,exp[l].sh);
      	for(h=0; h<nshapes; h++) {
      	  printf("%d ",exp[l].defs[h]);
      	}
      	printf("\n");
      }
    }
    /* printf("pts->ppm: %f\n",pts2ppm(85.760010,hed,1,0)); */
    /* printf("ppm->pts: %f\n",ppm2pts(3.98,hed,7,0)); */
    if(strcmp(token,"c\n")==0) {
      tot0=time(NULL);
      ctot0=clock();
      printf("correlation...\n");
      getCorrmat(corrmat,components,intervals,nshapes,x);

    }
    if(strcmp(token,"r\n")==0)
      printf("calling read FT2...\n");
    if(strcmp(token,"d\n")==0) {
      latestf=fopen("results/latestPeakList.txt","w");
      printf("calling decompose...\n");
      tot0=time(NULL);
      ctot0=clock();
      if(sel!=-1)
      	i=sel;
      else
      	i=0;

      //#pragma omp parallel 
      for(; i<intervals; i++) {

	rmsdprev=0;
	if(x1!=0 && x2!=0) {
	  interv[i].x1=x1;
	  interv[i].x2=x2;
	}

	r=abs(interv[i].x1-interv[i].x2);

	double ***shapes=(double ***) calloc(maxcomps,sizeof(double **));
	double ***shapes2=(double ***) calloc(maxcomps,sizeof(double **));
	
	for(g=0; g<maxcomps; g++) {
	  shapes[g]=(double **) calloc((nshapes-1),sizeof(double*));
	  shapes2[g]=(double **) calloc((nshapes-1),sizeof(double*));
	  for(k=0; k<nshapes-1; k++) {
	    shapes[g][k]=(double *) calloc(x,sizeof(double));
	    shapes2[g][k]=(double *) calloc(x,sizeof(double));
	  }
	}
	
	gsl_matrix *fdir=gsl_matrix_calloc(maxcomps,r+1);
	gsl_matrix *fdir2=gsl_matrix_calloc(maxcomps,r+1);
	      	
	if(comps && sel!=-1)
	  interv[i].comps=1; //comps;
	else
	  interv[i].comps=1;//getComps(tpeaklist,interv[i].x1,interv[i].x2,1);

	for(j=0; j<intervals; j++)
	  peaklist[j][2]=-1; // initialize for further peakpicking

	orgcomps=getComps(peaklist,interv[i].x1,interv[i].x2,1,intervals);
	printf("\ttotal identified components: %d\n",orgcomps);

	for(g=0; g<maxcomps; g++) {

	  comps=interv[i].comps;
	  int **peaklist2=(int **) malloc(comps*sizeof(int *));
	  for(j=0; j<comps; j++) {
	    peaklist2[j]=(int *) malloc(2*sizeof(int));
	    peaklist2[j][0]=peaklist2[j][1]=-1;
	  }

	  peaklist2[0][0]=peaklist[i][0];
	  peaklist2[0][1]=peaklist[i][1];

	  for(k=1,j=0; j<intervals; j++)
	    if(peaklist[j][2]==1 && k<comps && j!=i) {
	      peaklist2[k][0]=peaklist[j][0];
	      peaklist2[k][1]=peaklist[j][1];
	      k++;
	    }

	  t0=time(NULL);
	  c0=clock(); 
	
	  gsl_matrix *slice=gsl_matrix_calloc(x,r);
	  gsl_matrix *slice2=gsl_matrix_alloc(x,r);
	  gsl_matrix *S=gsl_matrix_calloc(x,r*exps);
	  gsl_matrix *R=gsl_matrix_calloc(x,r*exps);
  
	  gsl_vector *corr=gsl_vector_calloc(100); // hardcoded
	  gsl_matrix *cshapes=gsl_matrix_alloc(comps*exps,x);
	  double **ppmd=(double **) malloc(comps*sizeof(double *));
	  double ***ppm=(double ***) malloc(comps*sizeof(double **));
	  char ***snames=(char ***) malloc(comps*sizeof(char **));
	  char **dnames=(char **) malloc(comps*sizeof(char *));
	  
	  for(j=0; j<comps; j++) {
	    ppmd[j]=(double *) malloc((r+1)*sizeof(double));
	    ppm[j]=(double **) malloc((nshapes-1)*sizeof(double*));
	    snames[j]=(char **) malloc((nshapes-1)*sizeof(char*));
	    for(k=0; k<nshapes-1; k++) {
	      ppm[j][k]=(double *) malloc(x*sizeof(double));
	      snames[j][k]=(char *) malloc(257*sizeof(char));
	    }
	    dnames[j]=(char*) malloc(257*sizeof(char));
	  }
	  printf("\ttarget component: %d\n\tpts: %d %d\n\tppm: %2.2f %3.2f\n", \
		 0,peaklist[i][0],peaklist[i][1],\
		 pts2ppm(peaklist[i][0],hed,0,0),pts2ppm(peaklist[i][1],hed,1,1));
	  fprintf(flog,"\ttarget component: %d %d pts: %d %d ppm: %2.2f %3.2f\n",\
		  0,peaklist[i][2],peaklist[i][0],peaklist[i][1],pts2ppm(peaklist[i][0],\
		  hed,0,0),pts2ppm(peaklist[i][1],hed,1,1));
	  for(sub=0,k=1; k<comps; k++) {
	    if(peaklist2[k][1]!=-1) {
	      printf("\tadded component  : %d pts: %d %d ppm: %2.2f %3.2f\n", \
		     k,peaklist2[k][0],peaklist2[k][1],			\
		     pts2ppm(peaklist2[k][0],hed,0,0),pts2ppm(peaklist2[k][1],hed,1,1));
	      fprintf(flog,"\tadded component  : %d pts: %d %d ppm: %2.2f %3.2f\n",\
		      k,peaklist2[k][0],peaklist2[k][1],pts2ppm(peaklist2[k][0],\
		      hed,0,0),pts2ppm(peaklist2[k][1],hed,1,1));
	    }
	    else {
	      sub+=1;
	      printf("\trandom comp %d added\n",k);
	    }	    
	  }
	  for(h=0; h<iterations; h++) {
	    for(shape=0; shape<nshapes-1; shape++) {
	      decompose(spec,interv[i],shape,nshapes-1,\
			     exp,h,set,x,exps,corr,shapes,shapes2,fdir,fdir2,\
			cutoffs,peaklist2,cshapes,y,nshape);
	      for(comp=0; comp<comps; comp++) {
  		snprintf(filename1,100,"%s%s%d%s%d%s%d%s%d%s%d",rdir,"/interval_",i,\
			 "_shape_",interv[i].x1,"_",interv[i].x2,"_",comp,"_",shape);
		strcpy(snames[comp][shape],filename1);
  		snprintf(filename2,100,"%s%s%d%s%d%s%d%s%d",rdir,"/interval_",i,\
			 "_fdir_",interv[i].x1,"_",interv[i].x2,"_",comp);
		strcpy(dnames[comp],filename2);
	      }
	      if(h==0 && comps==1) {
		fprintf(rlist,"%s",snames[0][shape]);
		fprintf(rlist,"\n");
	      }
	    }
	    if(h==0 && comps==1) {
	      fprintf(rlist,"%s",dnames[0]);
	      fprintf(rlist,"\n");
	    }
	    reconstruct(interv[i],exps,exp,temp1s,temp3s,slice,slice2,R,S,comps,x,shapes,\
	    		0,r,spec,fdir,cutoffs,cshapes);
	    rmsd=rmsd1(S,R,r,0);
	    tres=rmsd*cut;
	    
	    printf("\tinterval: %d (%d) range: %d-%d comps: %d iteration: %d grand rmsd: %f correlation: %1.3f\n",\
		   i,interv[i].seq,interv[i].x1,interv[i].x2,comps,h,rmsd,gsl_vector_get(corr,h));
	    fprintf(flog,"\tinterval: %d %d range: %d-%d comps: %d iteration: %d grand rmsd: %f correlation: %1.3f\n",\
		    i,interv[i].seq,interv[i].x1,interv[i].x2,comps,h,rmsd,gsl_vector_get(corr,h));
	    
	    if(h && (tres>fabs(rmsdprev-rmsd) || h==br-1)) {
	      break;
	    }
	    else
	      rmsdprev=rmsd;
	  }
	  for(j=0; j<comps; j++) { 
	    for(k=0; k<nshapes-1; k++) {
	      for(l=0; l<x; l++) {
		ppm[j][k][l]=pts2ppm(l,hed,k+1,0);
	      }
	      if(strcmp(seq[i].residue,"PRO")==0)
		pc++;
	      if(minus[k])
		residue=seq[i+pc-1].residue; // one intervall for every residue
	      else
		residue=seq[i+pc].residue; // one intervall for every residue

	      /* peakpicker(shapes2[j][k],ppm[j][k],k+1,x,len1,i+pc,fpeaklist,stat,\ */
	      /* 		 seq,hed,expected[k],pcount,residue,minus[k]); */
	      fs=fopen(snames[j][k],"w");
	      for(l=0; l<x; l++) {
	      	sprintf(buffer1,"%f ",ppm[j][k][l]);
	      	sprintf(buffer2,"%f",shapes[j][k][l]);
		fprintf(fs,buffer1,"%f");
		fprintf(fs,buffer2,"%f");
		fprintf(fs,"\n");
	      }
	      fclose(fs);
	    }
	    fs=fopen(dnames[j],"w");

	    for(l=0; l<r+1; l++) 
	      ppmd[j][l]=pts2ppm(l+interv[i].x1,hed,0,1);
	    
	    if(j==0) {
	      fpeaklist[0].dpoint=getDir(fdir2,ppmd,interv[i].x1,r+1,1);
	      fpeaklist[0].dshift=getDir(fdir2,ppmd,interv[i].x1,r+1,0);
	      gsl_vector_view tm=gsl_matrix_row(fdir2,0);
	      fpeaklist[0].intens=gsl_vector_max(&tm.vector);
	    }

	    for(l=0; l<r+1; l++) {
	      sprintf(buffer1,"%f ",ppmd[j][l]);
	      sprintf(buffer2,"%f ",gsl_matrix_get(fdir2,j,l));
	      fprintf(fs,buffer1,"%f");
	      fprintf(fs,buffer2,"%f");
	      fprintf(fs,"\n");
	    }
	    /* gsl_vector_view tr=gsl_matrix_row(fdir2,j); */
	    /* gsl_vector_fprintf(fs,&tr.vector,"%f"); */
	    fclose(fs);
	  }
	  nch=0;
	  for(j=0; j<1; j++) {
	    //	  for(j=0; j<comps; j++) {
	    /* gsl_vector_view rt=gsl_vector_view_array(shapes[j][0],x); */
	    /* gsl_vector_memcpy(nshape,&rt.vector); */
	    nch=ncheck2(peaklist[i][1],nshape,x,cut2,j,flog);
	    if(nch==0) {
	      selected=j;
	      printf("\ttarget comp %d\n",selected);
	      break;
	    }
	  }
	  if(nch) { // hardcoded
	    ++interv[i].comps;
	    //	    fprintf(flog,"\n");
	    t1=time(NULL);
	    c1=clock();
	    printf("\twall clock time: %ld CPU time: %2.2f\n",(long) (t1 - t0),(float) (c1 - c0)/CLOCKS_PER_SEC);
	  }
	  else {
	    fprintf(peakl,"%d NH: %2.3f %d %2.2f\n",i,fpeaklist[0].dshift,fpeaklist[0].dpoint,\
		    fpeaklist[0].intens);
	    for(j=0; j<nshapes-1; j++) {
	      strcpy(fpeaklist[0].nuc,hed[j+1].nuc);
	      for(k=0; k<x; k++)
		finalShapes[selected][j][k]=shapes[selected][j][k];
	      l=simplePeakpicker(shapes[selected][j],ppm[selected][j],shape,x,fpeaklist,cut3); // hardcoded
	      if(l) 
		for(h=0; h<l; h++) {
		  fprintf(peakl,"%d %s: %3.2f %d %2.3f\n",i,fpeaklist[0].nuc,\
			  fpeaklist[h].shift,fpeaklist[h].point,fpeaklist[h].intens);
		  fprintf(latestf,"%d %s: %3.2f %d %2.3f\n",i,fpeaklist[0].nuc,\
			  fpeaklist[h].shift,fpeaklist[h].point,fpeaklist[h].intens);
		} 
	      else {
		fprintf(peakl,"%d %s: %d %d %d\n",i,fpeaklist[0].nuc,999,999,0);	      
		fprintf(latestf,"%d %s: %d %d %d\n",i,fpeaklist[0].nuc,999,999,0);	      
	      }
	    }
	    // glycine detection
	    if(strcmp(etype,"b")==0 &&\
	       glycineDetection(ppm,finalShapes[i],gsl_vector_get(cutoffs,getPlane(exp,exps,4)),\
				gsl_vector_get(cutoffs,getPlane(exp,exps,6)),stat,0,4,5,6,7,x,len1)) {
	      for(k=0,j=0; j<len2; j++)
		if(strcmp(seq[j].residue,"GLY")==0 && seq[j].component==0) {
		  seq[j].component=1;
		  k++;
		}
	      if(k==0)
		printf("\tall glycines detected\n");
      	    }
	    else
	      printf("\tnot a backbone experiment so no glycines detected\n");

	    t1=time(NULL);
	    c1=clock();
	    printf("\twall clock time: %ld CPU time: %f\n\t----------------------------------------------------------------------------------------\n",(long) (t1 - t0),(float) (c1 - c0)/CLOCKS_PER_SEC);
	    state=1;
	  }

	  gsl_vector_free(corr);
	  /* gsl_matrix_free(fdir); */
	  /* gsl_matrix_free(fdir2); */

	  gsl_matrix_free(slice);
	  gsl_matrix_free(slice2);
	  gsl_matrix_free(S);
	  gsl_matrix_free(R);
	  gsl_matrix_free(cshapes);
	  
	  for(l=0; l<comps; l++) { 
	    for(k=0; k<nshapes-1; k++) {
	      free(ppm[l][k]);
	      /* free(shapes[l][k]); */
	      /* free(shapes2[l][k]); */
	      free(snames[l][k]);
	    }
	    free(ppm[l]);
	    /* free(shapes[l]); */
	    /* free(shapes2[l]); */
	    free(snames[l]);
	  }
	  
	  free(ppm);
	  /* free(shapes); */
	  /* free(shapes2); */
	  free(snames);

	  for(l=0; l<comps; l++) {
	    free(peaklist2[l]);
	    free(ppmd[l]);
	    free(dnames[l]);
	  }
	  
	  free(ppmd);
	  free(peaklist2);
	  free(dnames);
	  if(comps==interv[i].comps)
	    break;
	}

	for(l=0; l<maxcomps; l++) { 
	  for(k=0; k<nshapes-1; k++) {
	    free(shapes[l][k]);
	    free(shapes2[l][k]);
	  }
	  free(shapes[l]);
	  free(shapes2[l]);
	}      
	free(shapes);
	free(shapes2);
	
	gsl_matrix_free(fdir);
	gsl_matrix_free(fdir2);
	
	if(i==sel)
	  break;
      }
            
      /* FILE *pipe=popen("gnuplot -persist","w"); */
      /* fprintf(pipe,"plot ./results/shape_257_260_0_1 with lines"); */
      /* close(pipe);  */
      // writing peaklist

      for(i=0; i<1000; i++)
	if(fpeaklist[i].counter)
	  printf("%d %f %f %s %d\n",shape,fpeaklist[i].shift,	\
		 fpeaklist[i].zeros,fpeaklist[i].nuc,fpeaklist[i].resdn);      

      fclose(latestf);
    }
    if(strcmp(token,"p\n")==0) {
      printf("calling plot ...\n");
      if(state) {
	;
	/* FILE *pipe=popen("gnuplot -persist","w"); */
	/* fprintf(pipe,"plot ./results/shape_257_260_0_1 with lines"); */
	/* close(pipe);  */
      }
      else 
	printf("nothing to plot\n");
    }
    if(strcmp(token,"s\n")==0)
      printf("calling sliding...\n");
    if(strcmp(token,"q\n")==0)
      tot0=time(NULL);
      ctot0=clock();
      break;
      //      break; // remove this when making options like q for quit
  }

  p=p;
  tot1=time(NULL);
  ctot1=clock();
  printf ("\ttotal wall clock time: %ld s %2.2f min\n",\
	  (long) (tot1 - tot0),(long) (tot1 - tot0)/60.0);
  printf ("\ttotal CPU time:        %f\n",\
	  (float) (ctot1 - ctot0)/CLOCKS_PER_SEC);
  //  fprintf(flog,"\n");
  fprintf(flog,"\n\ttotal wall clock time: %ld s %2.2f min\n",\
	  (long) (tot1 - tot0),(long) (tot1 - tot0)/60.0);
  fprintf(flog,"\ttotal CPU time:        %f\n\n",	\
	  (float) (ctot1 - ctot0)/CLOCKS_PER_SEC);

  // correlation
  if(strcmp(etype,"b")==0 && sel==-1) {
    correlation(chains,finalShapes,corrMatrix,x,seq,0.25);
    gsl_matrix *statsh=gsl_matrix_alloc(1000,x); // temp
    statShapes(statsh,stat,x,len1,hed,seq,len2);
    for(i=0; i<x; i++)
      gsl_vector_set(temp1s,i,gsl_matrix_get(statsh,2,i));
    for(i=0; i<x; i++)
      gsl_vector_set(temp1s,i,gsl_matrix_get(statsh,3,i));
    for(i=0; i<x; i++)
      gsl_vector_set(temp1s,i,gsl_matrix_get(statsh,4,i));
    // sliding
    // peakpicking
    gsl_matrix_free(statsh);
  }

  fclose(rlist);
  fclose(flog);
  fclose(peakl);
  //  fclose(latestf);
  gsl_vector_free(kf);
  gsl_vector_free(temp1s);
  gsl_vector_free(temp3s);
  gsl_vector_free(cutoffs);
  gsl_vector_free(nshape);  

  gsl_matrix_free(corrmat);
  gsl_matrix_free(corrMatrix);
  gsl_matrix_free(chains);
  gsl_vector_free(order);

  for(i=0; i<exps; i++) { 
    for(j=0; j<x; j++)
      free(spec[i][j]);
    free(spec[i]);
  }

  free(spec);

  free(interv);
  
  for(l=0; l<exps; l++) {
    free(exp[l].path);
    free(exp[l].defs);
    free(exp[l].sh);
    free(exp[l].pos);
    free(exp[l].sign);
  }

  free(defs);    
  free(exp);

  for(l=0; l<nshapes; l++) {
    free(hed[l].format);
    free(hed[l].nuc);      
  }

  for(l=0; l<intervals; l++) {
    for(k=0; k<nshapes-1; k++) 
      free(finalShapes[l][k]);
    free(finalShapes[l]);
  }
  
  free(finalShapes);
  free(hed);

  for(l=0; l<nshapes; l++)
    free(nuclei[l]);

  free(nuclei);

  for(i=0; i<all; i++)
    free(set[i]);

  free(set);

  for(i=0; i<intervals; i++)
    free(peaklist[i]);

  free(peaklist);

  for(i=0; i<1000; i++)
    free(tpeaklist[i]);

  free(tpeaklist);


  for(i=0; i<len1; i++) {
    free(stat[i].res);
    free(stat[i].nucs);
  }

  free(stat);

  for(i=0; i<len2; i++)
    free(seq[i].residue);

  free(seq);

  for(i=0; i<1000; i++)
    free(fpeaklist[i].nuc);

  free(fpeaklist);
  free(setup);
  free(buffer1);
  free(buffer2);
  free(tbuff);
  free(lbuff);
  free(pbuff);

  return 1;
}
