/* decomposition of spectra */

#include "headerfiles/common2.h"
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include "headerfiles/fastnnls.h"
#include "headerfiles/setRandom.h"
#include "headerfiles/startShape.h"
#include "headerfiles/common.h"
#include <time.h>
#include <stdlib.h>

int decompose(float ***spec,interval interv,\
	       int shape,int nshapes,\
	       experiments *exp,int iteration,\
	       int **set,int ind,int exps,\
	       gsl_vector *corr,\
	       double ***shapes,double ***shapes2,\
	       gsl_matrix *fdir,gsl_matrix *fdir2,\
	       gsl_vector *cutoffs,int **peaklist2,\
	       gsl_matrix *cshapes,int direct,\
	      gsl_vector *nshape) {
  
  int i,j,m;
  int k=0,t=1,l=0,br=0;
  int sign1=1,pl=0,pr,comp,plane=0;
  int x1=interv.x1,x2=interv.x2,comps=interv.comps;
  int index1,conv,mid=0,sign,sel=-1,ns=0;
  int r=abs(x2-x1)+1;

  float rfdir=0;
  float max=0,sum=0;
  float tres,value,cut=0;
  float tichonov=0.1;
  float a=0,b=0;

  double tol2=1.0e-15;

  gsl_vector *temp1=gsl_vector_alloc(ind);
  gsl_vector *temp2=gsl_vector_alloc(ind);
  gsl_vector *temp3=gsl_vector_alloc(ind);
  gsl_vector *temp4=gsl_vector_alloc(ind);

  cut=cut; // to get rid of [-Wunused-but-set-variable]
  max=max; // to get rid of [-Wunused-but-set-variable]

  // iteration==0
  if(iteration==0) {
    for(i=0; i<exps; i++) { 
      //      printf("%d %d %d\n",shape,exp[i].defs[shape+1],exp[i].conv==0);
      if(exp[i].defs[shape+1]>=1 && exp[i].conv==0) {
	sel=i;
  	break;
      }
    }
  }
  // iteration>0
  else  
    for(t=0,j=0; j<exps; j++) 
      if(set[shape][j]!=0 || j==0) {
  	t++;
	//	printf("T: %d j: %d shape: %d set[shape][j]: %d\n",t,j,shape,set[shape][j]);
      }

  double *vec2=malloc(ind*sizeof(double));

  gsl_vector *P=gsl_vector_calloc(r*t*ind);
  gsl_vector *Ptot=gsl_vector_calloc(r*t*ind);

  gsl_matrix *S2=gsl_matrix_calloc(exps*ind,comps);
  gsl_matrix *Pm=gsl_matrix_alloc(t*ind,r);
  gsl_matrix *P2m=gsl_matrix_alloc(exps*ind,r);

  gsl_matrix *M=gsl_matrix_calloc(t*ind,comps*ind);
  gsl_matrix *Mtot=gsl_matrix_alloc(r*t*ind,comps*ind);
  gsl_matrix *MtM=gsl_matrix_alloc(comps*ind,comps*ind);
  gsl_matrix *MtM2=gsl_matrix_alloc(comps,comps);

  gsl_vector *MtP2=gsl_vector_alloc(comps);
  gsl_vector *X=gsl_vector_calloc(comps*ind);
  gsl_vector *Y=gsl_vector_alloc(comps*ind);
  gsl_vector *X2=gsl_vector_alloc(comps);
  gsl_vector *Y2=gsl_vector_alloc(comps);
  gsl_vector *P2=gsl_vector_calloc(t*ind);
  gsl_vector *MtP=gsl_vector_alloc(comps*ind);
  
  int planes[2]={1,-1};

  if(sel==-1 && iteration==0) {
    printf("\twarning: need a single plane for initial setup of shape %d\n",shape);
    printf("\tfalling back to first plane (sel=0)\n");
    sel=0;
  }

  for(i=0; i<exps; i++) {
    sum=0;
    for(j=0; j<r; j++) 
      for(k=0; k<ind; k++) 
	sum+=spec[i][k][j+x1];
    if(sum==0)
      printf("\tslice i %d set[shape][i]: %d has no signal\n",i,set[shape][i]);
  }
  sum=0;

  if(iteration==0) 
    for(i=0; i<t; i++)
      for(j=0; j<r; j++) 
	for(k=0; k<ind; k++) 
	  gsl_matrix_set(Pm,k+i*ind,j,spec[sel][k][j+x1]);
  else {
    for(i=0; i<t; i++)
      for(j=0; j<r; j++) 
	for(k=0; k<ind; k++) 
	  gsl_matrix_set(Pm,k+i*ind,j,spec[set[shape][i]][k][j+x1]);
    m=0;
    for(j=0; j<r; j++) 
      for(i=0; i<t; i++)
	for(k=0; k<ind; k++) 
	  gsl_vector_set(Ptot,m++,spec[set[shape][i]][k][j+x1]);
  }

  m=0;
  if(iteration>=0) 
    for(i=0; i<t; i++) { 
      gsl_matrix_view tt=gsl_matrix_submatrix(Pm,i*ind,0,ind,r);
      if(gsl_matrix_isnull(&tt.matrix))
	printf("\twarning: shape %d in plane %d has no signals\n",shape,set[shape][i]);
      else
	;
	//  	printf("\tmax intesity in plane %d: %f\n",set[shape][i],gsl_matrix_max(&tt.matrix));
    }

  for(i=0; i<exps; i++)
    for(j=0; j<r; j++) 
      for(k=0; k<ind; k++)
	gsl_matrix_set(P2m,k+i*ind,j,spec[i][k][j+x1]);

  if(iteration==0 && shape==0) {
    for(j=0; j<nshapes; j++) 
      for(k=0; k<ind; k++)
	shapes2[comps-1][j][k]=shapes[comps-1][j][k]=0;
    for(j=0; j<r; j++) {
      gsl_matrix_set(fdir,comps-1,j,0);
      gsl_matrix_set(fdir2,comps-1,j,0);
    }
  }

  if(iteration && shape==0) {
    for(i=0; i<comps; i++)
      for(j=0; j<nshapes; j++) 
	  for(k=0; k<ind; k++)
	    shapes2[i][j][k]=shapes[i][j][k];
    gsl_matrix_memcpy(fdir2,fdir);
  }

  /* if(comps==1 && iteration==0 && shape==0) { */
  /*   gsl_vector_view te=gsl_vector_view_array(shapes[0][1],ind); */
  /*   writeTemp(&te.vector,"T1a"); */
  /* } */
  /* if(comps==2 && iteration==0 && shape==0) { */
  /*   gsl_vector_view te=gsl_vector_view_array(shapes[0][1],ind);  */
  /*   writeTemp(&te.vector,"T2a"); */
  /* } */
  /* if(comps==3 && iteration==0 && shape==0) {  */
  /*   gsl_vector_view te=gsl_vector_view_array(shapes[0][1],ind); */
  /*   writeTemp(&te.vector,"T3a"); */
  /* }   */

  // iteration 0 initial gueses
  if(iteration==0) {
    for(comp=comps-1; comp<comps; comp++) { 
      pr=peaklist2[comp][1];
      tres=gsl_vector_get(cutoffs,sel);
      if(tres==0)
	tres=1;
      mid=abs(peaklist2[comp][0]-x1);
      if((x1+mid>x2) || (x2-mid<x1) || (peaklist2[comp][0]==-1))
	mid=abs(x2-x1)/2;
      if(strncmp(exp[sel].sh,"+N",2)==0 && pr!=-1 && exp[sel].defs[0]==0 && exp[sel].defs[1]>=1 &&\
	 exp[sel].conv==0) {
	ns=shape;
	for(j=0; j<ind; j++) {
	  if((pr-3)<j && j<(pr+3))
	    shapes[comp][shape][j]=gsl_matrix_get(Pm,j,mid);
	  else
	    shapes[comp][shape][j]=0;
	}
	//	printf("\tcomp: %d slice %d setting up %s and fdir\n",comp,pr,exp[sel].sh);
	gsl_vector_view tt=gsl_vector_view_array(shapes[comp][shape],ind);
	if(gsl_vector_max(&tt.vector)<tres) {
	  printf("\twarning: comp: %d %d N shape has weak signal\n\
               replacing it with random values\n",comp,shape);
	  for(j=0; j<ind; j++) 
	    shapes[comp][shape][j]=((double) rand() / (double) RAND_MAX)*tres;
	  }
	for(j=0; j<r; j++) 
	  gsl_matrix_set(fdir,comp,j,gsl_matrix_get(Pm,pr,j));
	gsl_vector_view ts=gsl_matrix_row(fdir,comp);
	if(gsl_vector_max(&ts.vector)<tres) {
	  printf("\twarning: comp: %d fdir has a weak signal replacing it with random values\n",comp);
	  for(j=0; j<r; j++)
	    gsl_matrix_set(fdir,comp,j,((double) rand() / (double) RAND_MAX)*tres);	
	}
      }
      else {
	if(pr!=-1) {
	  //	  printf("\tcomp: %d slice %d setting up %s\n",comp,pr,exp[sel].sh);
	  for(j=0; j<ind; j++) 
	    shapes[comp][shape][j]=gsl_matrix_get(Pm,j,mid);
	  gsl_vector_view zc=gsl_vector_view_array(shapes[comp][shape],ind);
	  if(gsl_vector_isnull(&zc.vector)) {
	    printf("\twarning: comp: %d shape %d is zero replacing with random values\n",\
		   comp,shape);
	    for(j=0; j<ind; j++) 
	      gsl_vector_set(&zc.vector,j,((double) rand() / (double) RAND_MAX)*tres);
	    /* for(j=0; j<r; j++) */
	    /*   gsl_matrix_set(fdir,comp,j,gsl_matrix_get(Pm,pr,j)); */
	  //gsl_matrix_set(fdir,comp,j,((double) rand() / (double) RAND_MAX)*tres);
	  }
	}
	else {
	  printf("\tcomp: %d initialising random values for %s\n",comp,exp[sel].sh);
	  for(j=0; j<ind; j++) 
	    shapes[comp][shape][j]=((double) rand() / (double) RAND_MAX)*tres;
	  if(shape==0) {
	    printf("\tcomp: %d initialising random values for fdir\n",comp);
	    for(j=0; j<r; j++)
	      gsl_matrix_set(fdir,comp,j,((double) rand() / (double) RAND_MAX)*tres);
	  }
	}
      }
    }
    
    /* reverse convolution */
    if(shape==nshapes-1 && comps>1) {
      for(k=comps-1; k<comps; k++) {
	if(peaklist2[k][1]!=-1) {
	  for(m=1; m<nshapes; m++) {
	    printf("\treverse convolution for comp: %d shape %d\n",k,m); 
	    mid=peaklist2[k][0];
	    for(sign=0; sign<2; sign++) {
	      sign1=planes[sign];
	      pl=getPlanes(exp,exps,m,sign1);
	      if(pl) { // note: if conv>0 and pl==0
		gsl_vector_view sn=gsl_vector_view_array(shapes[k][0],ind); // N hardcoded
		for(j=0; j<ind; j++)
		  gsl_vector_set(temp1,j,spec[pl][j][mid]);
		if(sign==0) {
		  gsl_vector_set_zero(temp2);
		  rconvolution(&sn.vector,temp1,temp2,sign1*-1,ind%2?ind/2:ind/2-1);
		}
		if(sign==1) {
		  gsl_vector_set_zero(temp3);
		  gsl_vector_set_zero(temp4);
		  rconvolution(&sn.vector,temp1,temp3,sign1*-1,ind%2?ind/2:ind/2-1);
		  for(j=0; j<ind; j++) {
		    a=gsl_vector_get(temp2,j);
		    b=gsl_vector_get(temp3,j);
		    if(a && b)
		      gsl_vector_set(temp4,j,(a+b)/2);
		  }
		}
	      }
	    }
	    plane=getPlane(exp,exps,m);
	    tres=gsl_vector_get(cutoffs,plane);
	    if(tres==0)
	      tres=1;
	    if(gsl_vector_isnull(temp4)==0 && pl)
	      for(j=0; j<ind; j++)
		shapes[k][m][j]=gsl_vector_get(temp4,j);
	    else {
	      printf("\tzero signal in comp %d shape %d replacing with random values\n",k,m);
	      for(j=0; j<ind; j++) 
		shapes[k][m][j]=((double) rand() / (double) RAND_MAX)*tres;
	      for(j=0; j<r; j++)
	      	gsl_matrix_set(fdir,k,j,((double) rand() / (double) RAND_MAX)*tres);
	    }
	  }
	}
      }
    }
  }

  // iteration > 0
  else {

    /* gsl_vector_view te=gsl_vector_view_array(shapes[0][1],ind); */
    /* if(comps==1)  */
    /*   writeTemp(&te.vector,"T1b"); */
    /* if(comps==2)  */
    /*   writeTemp(&te.vector,"T2b"); */
    /* if(comps==3)  */
    /*   writeTemp(&te.vector,"T3b"); */
    
    // direct dim
    if(shape==0) {
      for(i=0; i<exps; i++) {
	for(l=0; l<comps; l++) {
	  gsl_vector_view st=gsl_matrix_subcolumn(S2,l,i*ind,ind);
	  gsl_vector_view tt=gsl_matrix_row(cshapes,l+i*comps);
	  gsl_vector_memcpy(&st.vector,&tt.vector);
	}
      }
      gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,S2,S2,0.0,MtM2);	
      for(j=0; j<r; j++) {	
	gsl_vector_view P2r=gsl_matrix_column(P2m,j);
	if(gsl_vector_isnull(&P2r.vector))
	  printf("warning: slice %d is zero\n",j+x1);
      	gsl_blas_dgemv(CblasTrans,1.0,S2,&P2r.vector,0.0,MtP2);
	fastnnls(MtM2,MtP2,X2,Y2,tol2);
	for(l=0; l<comps; l++)
	  gsl_matrix_set(fdir,l,j,gsl_vector_get(X2,l));
      }
      for(l=0; l<comps; l++) {
      	gsl_vector_view se=gsl_matrix_row(fdir,l);
      	if(gsl_vector_isnull(&se.vector)) {
      	  printf("\tdirect: warning: comp %d fdir is zero, initialize to 1\n",l);
	  for(j=0; j<r; j++)
	    gsl_vector_set(&se.vector,j,1);
	}
      }
    }
    // indirect shapes   
    for(l=0; l<comps; l++) {
      rfdir=0;
      gsl_vector_view td=gsl_matrix_row(fdir,l);
      if(gsl_vector_isnull(&td.vector)==0)
	for(j=0; j<r; j++) 
	  rfdir+=gsl_vector_get(&td.vector,j);
      else {
	printf("\tcomp %d: fdir zero, initialize to 1\n",l);
	rfdir=1;
      }
      for(i=0; i<t; i++) {
	index1=set[shape][i];
	conv=exp[index1].conv;
	gsl_matrix_view di=gsl_matrix_submatrix(M,i*ind,l*ind,ind,ind);
	if(conv==0)
	  gsl_matrix_set_identity(&di.matrix);
	else
	  convolve4a(temp1,temp2,exp,&di.matrix,shapes2[l],shape,index1,ind,nshapes);
	if(gsl_matrix_isnull(&di.matrix)) {
	  printf("\terror: comp: %d shape %d: midmatrix %d is a zero matrix\
                   using a unit matrix\n",l,shape,i);
	  gsl_matrix_set_identity(&di.matrix);
	}
      }
    }

    for(j=0; j<r; j++) {
	gsl_matrix_view subt=gsl_matrix_submatrix(Mtot,j*t*ind,0,t*ind,comps*ind);
	gsl_matrix_view subm=gsl_matrix_submatrix(M,0,0,t*ind,comps*ind);
	gsl_matrix_memcpy(&subt.matrix,&subm.matrix);
      }

    for(l=0; l<comps; l++)
      for(j=0; j<r; j++) {
	gsl_matrix_view subt=gsl_matrix_submatrix(Mtot,j*t*ind,l*ind,t*ind,ind);
	gsl_matrix_scale(&subt.matrix,gsl_matrix_get(fdir,l,j));
      }
			 
    for(j=0; j<r; j++) {
      gsl_vector_view pd=gsl_matrix_column(Pm,j);
      gsl_vector_add(P2,&pd.vector);
    }
    //    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,M,M,0.0,MtM);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Mtot,Mtot,0.0,MtM);
    //    gsl_blas_dgemv(CblasTrans,1.0,M,P2,0.0,MtP);
    gsl_blas_dgemv(CblasTrans,1.0,Mtot,Ptot,0.0,MtP);

    for(i=0; i<comps*ind; i++) {
      value=gsl_matrix_get(MtM,i,i);
      gsl_matrix_set(MtM,i,i,value+tichonov*tichonov);
    }
    fastnnls(MtM,MtP,X,Y,tol2);
    for(l=0; l<comps; l++) {
      gsl_vector_view re=gsl_vector_subvector(X,l*ind,ind);
      if(gsl_vector_isnull(&re.vector)) {
	printf("\terror: comp: %d shape: %d is zero using the previous result\n",l,shape);
	for(j=0; j<ind; j++)
	  gsl_vector_set(&re.vector,j,shapes2[l][shape][j]);
      }
      else
	;
      if(shape==ns && l==0) {
	gsl_vector_memcpy(nshape,&re.vector);
	gsl_vector_view rt=gsl_vector_view_array(shapes2[l][shape],ind);
	max=gsl_vector_max(&rt.vector);
	for(j=0; j<ind; j++) {
	  vec2[j]=gsl_vector_get(&re.vector,j);
	}
	gsl_vector_set(corr,iteration,cosinesim(shapes2[l][shape],vec2,ind));
	//	for(j=0; j<ind; j++)
	//	printf("similarity: %d %f ",iteration,gsl_vector_get(corr,iteration));
	/* for(j=0; j<ind; j++) */
	/*   shapes[l][shape][j]=gsl_vector_get(&re.vector,j); */
      }
      else
	for(j=0; j<ind; j++)
	  shapes[l][shape][j]=gsl_vector_get(&re.vector,j);
    } 
  }

  if(shape==nshapes-1) {
    gsl_matrix_set_zero(cshapes);
    for(i=0; i<exps; i++)
      if(exp[i].conv)
    	for(j=0; j<comps; j++) {
    	  gsl_vector_view st=gsl_matrix_row(cshapes,j+i*comps);
	  //	  printf("comp: %d zero? %d\n",j,gsl_vector_isnull(&st.vector)?1:0);
    	  convolve4b(temp1,&st.vector,exp,shapes[j],i,ind,nshapes);
	  //	  printf("comp: %d zero? %d\n",j,gsl_vector_isnull(&st.vector)?1:0);
    	}
      else
    	for(j=0; j<comps; j++) {
    	  gsl_vector_view st=gsl_matrix_row(cshapes,j+i*comps);
    	  for(k=0; k<ind; k++)
    	    gsl_vector_set(&st.vector,k,shapes[j][exp[i].pos[0]][k]);
    	}
  }
  
  free(vec2);
  
  gsl_matrix_free(Pm);
  gsl_matrix_free(P2m);
  gsl_matrix_free(S2);  

  gsl_matrix_free(Mtot);
  gsl_vector_free(Ptot);

  gsl_vector_free(MtP);
  gsl_vector_free(MtP2);
  gsl_vector_free(P);
  gsl_vector_free(P2);
  gsl_vector_free(temp1);
  gsl_vector_free(temp2);
  gsl_vector_free(temp3);
  gsl_vector_free(temp4);

  gsl_matrix_free(M);
  gsl_matrix_free(MtM2);
  gsl_matrix_free(MtM);

  gsl_vector_free(X);
  gsl_vector_free(Y);
  gsl_vector_free(Y2);
  gsl_vector_free(X2);

  return br;
}
