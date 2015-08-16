/* utility functions */

#include "headerfiles/fastnnls.h"
#include "headerfiles/common.h"
#include "headerfiles/common2.h"
#include <gsl/gsl_statistics.h>
#include <regex.h>

void writeShape(int comp,int shape,int ind,double ***shapes,char *name) {

  gsl_vector *temp=gsl_vector_calloc(ind);
  
  FILE *file; 
  int j;

  file=fopen(name,"w");

  for(j=0; j<ind; j++)
    gsl_vector_set(temp,j,shapes[comp][shape][j]);
  gsl_vector_fprintf(file,temp,"%f");
  fclose(file);
  
  gsl_vector_free(temp);

}

void writeTemp(gsl_vector *temp,char *name) {
  
  FILE *file;
  file=fopen(name,"w");
  gsl_vector_fprintf(file,temp,"%f");

  fclose(file);
}

float minmaxShape(double *shape,int ind,int max) {
  
  gsl_vector_view stemp=gsl_vector_view_array(shape,ind); 

  if(max)
    return gsl_vector_max(&stemp.vector);
  else
    return gsl_vector_min(&stemp.vector);
}

int minmaxIndexShape(double *shape,int ind,int max) {

  gsl_vector_view te=gsl_vector_view_array(shape,ind);

  if(max)
    return gsl_vector_max_index(&te.vector);
  else
    return gsl_vector_min_index(&te.vector);
}

void writeTrace(gsl_vector *P,int offset,int ind,char *name) {

  FILE *file;

  file=fopen(name,"w");
  gsl_vector_view t=gsl_vector_subvector(P,offset,ind);
  gsl_vector_fprintf(file,&t.vector,"%f");
  fclose(file);
}

void writeSpec(float ***spec,int plane,int ind,int offset,char *name) {

  FILE *file;
  int j;
  gsl_vector *temp=gsl_vector_alloc(ind);

  file=fopen(name,"w");

  for(j=0; j<ind; j++)
    gsl_vector_set(temp,j,spec[plane][j][offset]);
  gsl_vector_fprintf(file,temp,"%f");
  fclose(file);
  gsl_vector_free(temp);
}


float maxTrace(gsl_vector *P,int offset,int ind) {

  gsl_vector_view t=gsl_vector_subvector(P,offset,ind);

  return gsl_vector_max(&t.vector);
}

float avgTrace(gsl_vector *P,int offset,int ind) {

  int i;
  float m=0;

  gsl_vector_view t=gsl_vector_subvector(P,offset,ind);
  
  for(i=0; i<ind; i++)
    m+=gsl_vector_get(&t.vector,i);
  
  return m/ind;
}

void convolve1(gsl_vector *temp1,gsl_vector *temp2,experiments *exp,\
	       gsl_matrix *midmatrix,double **shapes2,int shape,int index1,int ind) {

  int h,j,i,t;
  int sign1=1;
  int p,q,mid,mida;
  float v;

  /* mid=ind%2?(ind+1)/2-1:ind/2-1; */
  /* mida=mid; */

  if(ind%2) {
    mid=(ind+1)/2;
    mida=mid-1;
  }
  else {
    mid=ind/2;
    mida=mid;
  }

  /* if(ind%2) */
  /*   mida=mid=ind/2; */
  /* else */
  /*   mida=mid=ind/2-1; */
  /* P(mida) */
  /*   P(mid) */
  gsl_matrix_set_identity(midmatrix);
  q=exp[index1].conv;
  int *order1=malloc(sizeof(int)*q);
  int *order2=malloc(sizeof(int)*q);

  for(j=0,h=0; h<q; h++) {
    p=exp[index1].pos[h];
    if(p!=shape) {
      order1[j]=p;
      order2[j]=exp[index1].sign[h];
      j++;
    }
    else
      sign1=exp[index1].sign[h];
  }

  order1[q-1]=shape;
  order2[q-1]=sign1;

  for(h=0; h<q-1; h++) {
    gsl_vector_set_zero(temp1);
    gsl_vector_set_zero(temp2);
    p=order1[h];
    for(j=0; j<ind; j++)
      gsl_vector_set(temp1,j,shapes2[p][j]);
    if(h==0 && order2[h]<0)
      gsl_vector_reverse(temp1);
    gsl_blas_dgemv(CblasNoTrans,1.0,midmatrix,temp1,0.0,temp2);
    gsl_vector_memcpy(temp1,temp2);
    sign1=order2[h+1];
    //midmatrix4(midmatrix,temp1,sign1,sign1);
    //midmatrix5(midmatrix,temp1,sign1,mid,ind);
    midmatrix6(midmatrix,temp1,sign1,mid,mida,ind);
    //    midmatrix7(midmatrix,temp1,sign1,ind,mid);
    //    printf("no tilt:  %d %d %d %d %d %s\n",shape,p,index1,sign1,order2[h],exp[index1].sh);
    if(0 && order2[h]==-1 && shape!=0 && exp[index1].pos[h]) { // tilt check
      printf("tilt -1: %d %d %d %d %d %s\n",shape,p,index1,sign1,order2[h],exp[index1].sh);
      gsl_matrix *tempMid=gsl_matrix_alloc(ind,ind);
      for(i=0; i<ind; i++)
	for(j=0; j<ind; j++) {
	  v=gsl_matrix_get(midmatrix,i,j);
	  t=j-2*j;
	  if(t<0)
	    t=ind+t;
	  gsl_matrix_set(tempMid,i,t,v);
	}
      gsl_matrix_memcpy(midmatrix,tempMid);
      //      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,midmatrix,midmatrix,0.0,tempMid);
      gsl_matrix_free(tempMid);
    }
    if(0 && order2[h]==1 && shape!=0 && exp[index1].pos[h]) { // tilt check
      printf("tilt 1: %d %d %d %d %d %s\n",shape,p,index1,sign1,order2[h],exp[index1].sh);
      gsl_matrix *tempMid=gsl_matrix_alloc(ind,ind);
      for(i=0; i<ind; i++)
	for(j=0; j<ind; j++) {
	  v=gsl_matrix_get(midmatrix,i,j);
	  t=2*j;
	  if(t<0)
	    t=ind+t;
	  gsl_matrix_set(tempMid,i,t%ind,v);
	}
      gsl_matrix_memcpy(midmatrix,tempMid);
      //      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,midmatrix,midmatrix,0.0,tempMid);
      gsl_matrix_free(tempMid);
    }
  }
  free(order1);
  free(order2);
}

void convolve2(gsl_vector *temp1,gsl_vector *temp2,experiments *exp,\
	       gsl_matrix *midmatrix,double ***shapes2,int comp,int index1) {

  int h,j;
  int ind=temp1->size;
  int p,q,mid;

  mid=ind%2?(ind+1)/2-1:ind/2-1;    

  gsl_matrix_set_identity(midmatrix);

  q=exp[index1].conv;

  for(h=0; h<q; h++) {
    p=exp[index1].pos[h];
    for(j=0; j<ind; j++)
      gsl_vector_set(temp1,j,shapes2[comp][p][j]);
    gsl_blas_dgemv(CblasNoTrans,1.0,midmatrix,temp1,0.0,temp2);
    gsl_vector_memcpy(temp1,temp2);
    if(h+1<q) 
      midmatrix5(midmatrix,temp1,exp[index1].sign[h+1],mid,ind);
    else
      break;
  }
}

void convolve3(gsl_vector *temp1,gsl_vector *temp2,gsl_vector *temp3,\
	       experiments *exp,double **shapes,int index1) {

  int h,j;
  int ind=temp1->size;
  int p,q,mid,mida,sign;

  //  mid=ind%2?(ind+1)/2:ind/2;

  if(ind%2) {
    mid=(ind+1)/2;
    mida=mid-1;
  }
  else {
    mid=ind/2;
    mida=mid;
  }

  q=exp[index1].conv;
  p=exp[index1].pos[0];

  for(j=0; j<ind; j++)
    gsl_vector_set(temp1,j,shapes[p][j]);
    
  for(h=1; h<q; h++) {
    p=exp[index1].pos[h];
    gsl_vector_view t2=gsl_vector_view_array(shapes[p],ind);
    sign=exp[index1].sign[h];
    //    convolution(temp1,&t2.vector,temp3,sign,mid,ind);
    convolution2(temp1,&t2.vector,temp3,sign,mid,mida,ind);
    gsl_vector_memcpy(temp1,temp3);
  }
}

void convolve4a(gsl_vector *temp1,gsl_vector *temp2,experiments *exp,gsl_matrix *midmatrix,\
	       double **shapes2,int shape,int index1,int ind,int nshapes) {

  int h,i,j,p;
  int mid;

  if(ind%2) // odd
    mid=ind/2;
  else  // even
    mid=ind/2-1;

  gsl_matrix_set_identity(midmatrix);

  for(h=0; h<nshapes; h++) {
    p=exp[index1].defs[h+1];
    if(p!=0 && h!=shape) {
      //      printf("shape: %d index1: %d h: %d p: %d\n",shape,index1,h,p);
      for(i=0; i<ind; i++)
	gsl_vector_set(temp1,i,shapes2[h][i]);
      //      gsl_vector_swap_elements(temp1,mid,0);
      if(p==-1)
      	gsl_vector_reverse(temp1);
      gsl_blas_dgemv(CblasNoTrans,1.0,midmatrix,temp1,0.0,temp2);

      for(i=0; i<ind; i++)
	for(j=0; j<ind; j++) 
	  gsl_matrix_set(midmatrix,i,j,gsl_vector_get(temp2,modX((mid-j+i),ind)));
    }
  }

  if(exp[index1].defs[shape+1]==-1)
    for(i=0; i<mid+1; i++) // note: changed 21-1-15
      gsl_matrix_swap_columns(midmatrix,i,ind-1-i);
}

void convolve4b(gsl_vector *temp1,gsl_vector *temp2,experiments *exp,double **shapes,\
		int index1,int ind,int nshapes) {

  int h,i,j,p;
  int mid;

  if(ind%2) // odd
    mid=ind/2;
  else  // even
    mid=ind/2-1;

  gsl_matrix *midma=gsl_matrix_alloc(ind,ind);
  gsl_matrix_set_identity(midma);

  for(h=0; h<nshapes; h++) {
    p=exp[index1].defs[h+1];
    if(p!=0) {
      for(i=0; i<ind; i++)
	gsl_vector_set(temp1,i,shapes[h][i]);
      if(p==-1)
      	gsl_vector_reverse(temp1);
      gsl_blas_dgemv(CblasNoTrans,1.0,midma,temp1,0.0,temp2);
      for(i=0; i<ind; i++)
	for(j=0; j<ind; j++) 
	  gsl_matrix_set(midma,i,j,gsl_vector_get(temp2,modX((mid-j+i),ind)));
    }
  }
  gsl_matrix_free(midma);
}

void convolution(gsl_vector *shape1,gsl_vector *shape2,gsl_vector *cshape,int sign,\
		 int mid,int le) {

  int k,j;
  float sum;

  if(sign>0) {
    gsl_vector_reverse(shape1);
   #pragma omp parallel default(shared) private(k,j)
   {
   #pragma omp for reduction(+:sum)
    for(k=0; k<le; k++) {
      sum=0;
      for(j=0; j<le; j++) 
	sum+=gsl_vector_get(shape1,modX(j+mid-k,le))*gsl_vector_get(shape2,j);
      gsl_vector_set(cshape,k,sum);
    }
   }
    gsl_vector_reverse(shape1);
  }
  if(sign<0) {
  #pragma omp parallel default(shared) private(k,j)
  {
  #pragma omp for reduction(+:sum)
    for(k=0; k<le; k++) {
      sum=0;
      for(j=0; j<le; j++) 
	sum+=gsl_vector_get(shape1,modX(j-mid+k,le))*gsl_vector_get(shape2,j);
      gsl_vector_set(cshape,k,sum);
    }
  }
  }
}

void convolution2(gsl_vector *shape1,gsl_vector *shape2,gsl_vector *cshape,int sign,\
		 int mid,int mida,int le) {

  int i,j,k;
  int midb=mida-1;
  float sum;

  if(sign>0) {
    for(i=0; i<le; i++) {
      sum=0;
      for(j=0; j<le; j++) {
	k=midb-j+i;
	if(k>=le) 
	  k=k-le;
	if(k<0) 
	  k=k+le;
	sum+=gsl_vector_get(shape1,k)*gsl_vector_get(shape2,j);
	}
      gsl_vector_set(cshape,i,sum);
      }
    } 
 
  if(sign<0) {
    for(i=0; i<le; i++) {
      sum=0;
      for(j=0; j<le; j++) { 
	k=j+i-mid;
	if(k>=le) 
	  k=k-le;
	if(k<0) 
	  k=k+le;
	sum+=gsl_vector_get(shape1,k)*gsl_vector_get(shape2,j);
      }
      gsl_vector_set(cshape,i,sum);
    }
  }
} 

void midmatrix4(gsl_matrix *M,gsl_vector *temp1,int sign,int swap) {

  int i,j;
  int le=temp1->size;
  int mid;//=le/2-1; // NOTE! changed from le/2 to le/2-1 because of 0-191 not 1-192 (start at zero)

  mid=le%2?(le+1)/2-1:le/2-1;
  
  if(sign>0) {
    gsl_vector_reverse(temp1);
    for(i=0; i<le; i++) 
      for(j=0; j<le; j++) 
	gsl_matrix_set(M,i,j,gsl_vector_get(temp1,(j+mid-i+le)%le));  
    gsl_vector_reverse(temp1);
  }
  else { 
    for(i=0; i<le; i++) 
      for(j=0; j<le; j++) 
	gsl_matrix_set(M,i,j,gsl_vector_get(temp1,(j+mid-i+le)%le));
    if(swap<0)
      for(i=0; i<mid; i++)
	gsl_matrix_swap_rows(M,i,le-1-i);
  } 
}

void midmatrix5(gsl_matrix *M,gsl_vector *temp1,int sign,int mid,int le) {

  int i,j;

  gsl_vector_reverse(temp1);
#pragma omp parallel default(shared) private(i,j)
  {
  //pragma omp for reduction()
  for(i=0; i<le; i++)
    for(j=0; j<le; j++) 
      gsl_matrix_set(M,i,j,gsl_vector_get(temp1,modX(j+mid-i,le)));
  }
  gsl_vector_reverse(temp1);

  if(sign<0) 
    for(i=0; i<mid; i++)
      gsl_matrix_swap_columns(M,i,le-1-i);
  
}

void midmatrix6(gsl_matrix *M,gsl_vector *temp1,int sign,int mid,int mida,int le) {

  int i,j,k;

  gsl_vector_reverse(temp1);
  for(i=0; i<le; i++)
    for(j=0; j<le; j++) {
      if(abs(j-i)<mida)
	//     if(j<le-mid+i && j>=i-mid)
	k=j-i+mid;
      else {
	k=j-i-mida;
	if(k<-le)
	//	if(j>=i-mid)
	//  k=j-i-mida;
	  k=k+2*le;
        if(k<0) 
	  k=k+le;
	/* else */
	/*   k=j-i-mida+2*le; */
      }
      gsl_matrix_set(M,i,j,gsl_vector_get(temp1,k));
    }
  gsl_vector_reverse(temp1);
  
  if(sign<0)
    for(i=0; i<mid; i++)
      gsl_matrix_swap_columns(M,i,le-1-i);
}

void midmatrix7(gsl_matrix *M,gsl_vector *temp1,int sign,int le,int mid) {

  int i,j;

  // midmatix from supporting information

  /* if(sign<0) */
  /*   gsl_vector_reverse(temp1); */

  //gsl_vector_reverse(temp1);

  for(i=0; i<le; i++)
    for(j=0; j<le; j++) 
      gsl_matrix_set(M,i,j,gsl_vector_get(temp1,modX((mid-j+i),le)));
  
  //gsl_vector_reverse(temp1);

  /* if(sign<0) */
  /*   for(i=0; i<mid; i++) */
  /*     gsl_matrix_swap_columns(M,i,le-1-i); */

  /* if(sign<0) */
  /*   for(i=0; i<mid; i++) */
  /*     gsl_matrix_swap_rows(M,i,le-1-i); */
}

void rconvolution(gsl_vector *nshape,gsl_vector *temp1,gsl_vector *temp2,int sign,int mid) {

  int k,j;
  int le=temp1->size;
  float sum=0;

  for(k=0; k<le; k++) {
    sum=0;
    for(j=0; j<le; j++) {
      sum+=gsl_vector_get(nshape,j)*gsl_vector_get(temp1,modX(j+mid-k,le));
    }
    gsl_vector_set(temp2,k,sum);
  }
  if(sign<0)
    gsl_vector_reverse(temp2);
}

int getPlane(experiments *exp,int exps2,int shape) {

  int i;
  
  for(i=0; i<exps2; i++)
    if(exp[i].conv==0 && exp[i].pos[0]==shape) {
      return i;
    }
  return 0;     
}

int getPlanes(experiments *exp,int exps2,int shape,int sign) {

  int i;
  
  for(i=0; i<exps2; i++)
    if(exp[i].conv==2 && exp[i].defs[1] && exp[i].defs[shape+1]==sign)
      return i;

  return 0;
}

float kfactor(gsl_matrix *M,gsl_vector *nshape,gsl_vector *P,float cutoff) {
  
  int i=0,j=0;
  int n=M->size1;
  float sum1=0,sum2=0; //,s1=0,s2=0;
  //  printf("min: %f\n",gsl_vector_min(P));
  gsl_vector *temp=gsl_vector_alloc(n);
  gsl_blas_dgemv(CblasNoTrans,1.0,M,nshape,0.0,temp);
  //  normalizeTemp(temp);
  for(j=0; j<n; j++) {
     /* if(1) { */
     if(gsl_vector_get(temp,j)>cutoff && gsl_vector_get(P,j)>cutoff) {
      i++;
      sum1+=gsl_vector_get(temp,j);
      sum2+=gsl_vector_get(P,j);
    }
  }
  
  gsl_vector_free(temp);
  
  if(sum1==0 || sum2==0)
    sum1=sum2=1.0;
  //  printf("%d\n",i);
  /* sum1/=j; */
  /* sum2/=j; */
  
  return sum2/sum1;
}

float subFactor1(gsl_vector *temp3, float ***shapes2,int comp,int shape,int ind,float cutoff) {

  // subtraction factor taking into account width of shape
  int i,j;
  int maxi=0;
  float m,max=0;
  gsl_vector *temp=gsl_vector_alloc(3);

  for(j=0; j<ind; j++)
    if(max<shapes2[comp][shape][j]) {
      max=shapes2[comp][shape][j];
      maxi=j;
    }

  for(j=0,i=maxi-1; i<maxi+2; i++)
    if(i>0 && i<temp3->size && shapes2[comp][shape][i]>cutoff && gsl_vector_get(temp3,i)>cutoff)
      gsl_vector_set(temp,j++,gsl_vector_get(temp3,i)/shapes2[comp][shape][i]);
    else
      gsl_vector_set(temp,j++,1);

  m=gsl_vector_max(temp);
  gsl_vector_free(temp);

  return m;
}

float subFactor2(gsl_vector *temp3, gsl_vector *st,float cutoff) {

  int i,j;
  float m;
  int maxi;
  gsl_vector *temp=gsl_vector_alloc(3);
  maxi=gsl_vector_max_index(temp3);

  for(j=0,i=maxi-1; i<maxi+2; i++)
    if(i>0 && i<temp3->size && gsl_vector_get(st,i)>cutoff &&  gsl_vector_get(temp3,i)>0)
      gsl_vector_set(temp,j++,gsl_vector_get(st,i)/gsl_vector_get(temp3,i));
    else 
      gsl_vector_set(temp,j++,1);
 
  m=gsl_vector_max(temp);
  gsl_vector_free(temp);

  return m;
}

void normalizeTemp(gsl_vector *temp) {

  float max=gsl_vector_max(temp);
  gsl_vector_scale(temp,1/max);

}

void normalizeShape(float *** shapes,int ind,int comp,int shape) {

  int j;
  float max=0;

  for(j=0; j<ind; j++)
    if(shapes[comp][shape][j]>max)
      max=shapes[comp][shape][j];

  for(j=0; j<ind; j++)
    shapes[comp][shape][j]/=max;

}

void normalizeShape2(gsl_vector *f,int ind,int norm) {

  int i;
  float sum=0;

  if(norm==0)
    sum=gsl_vector_max(f);
  else
    for(i=0; i<ind; i++)
      sum+=gsl_vector_get(f,i);

  //  sum*=2;

  if(sum)
    gsl_vector_scale(f,1.0/sum);
  else
    printf("zero sum, no normalization\n");
}

float nCutoff(gsl_vector *P) {

  int i,j,top;
  int t=P->size;
  double mean,std;
  double slices[t];

  for(i=0; i<t; i++)
    slices[i]=gsl_vector_get(P,i);
   
  mean=0;

  /* for(i=0; i<t; i++) */
  /*   if(gsl_vector_get(P,i)>0) { */
  /*     mean+=gsl_vector_get(P,i); */
  /*     j++; */
  /*   } */
 
  qsort(slices,t,sizeof(double),floatcompare);

  top=(int) t*0.95;
  double data2[top];

  for(j=0; j<top; j++)
    data2[j]=slices[j];

  mean=gsl_stats_mean(data2,1,top);
  std=gsl_stats_sd_m(data2,1,top,mean);

  std=std; // just to avoid compiler messages
  if(mean==0)
    mean=0.01;
    
  return mean;
}

void lineBroadening(gsl_vector *X,int f,float cutoff) {

  int i,j;
  int n=X->size;
  float sum;
  gsl_vector *temp=gsl_vector_alloc(n);
  gsl_vector_memcpy(temp,X);

  for(i=0; i<n; i++) {
    sum=0;
    for(j=-f+i; j<=f+i; j++) {
      //      printf("%d: %d %f\n",i,j,sum);
      if(j>=0 && j<n && (j==-f+i || j==f+i))
	sum+=gsl_vector_get(temp,j)*0.25;
      if(j>=0 && j<n && j==i)
	sum+=gsl_vector_get(temp,j)*0.50;
    }
    //    printf("\n-------------------------------\n");
    if(sum && gsl_vector_get(X,i)>cutoff)
      gsl_vector_set(X,i,sum);
    else
      gsl_vector_set(X,i,0);
  }
  gsl_vector_free(temp);
}

float normFactor(gsl_matrix* midmatrix,gsl_vector* temp) {

  return pow(gsl_matrix_max(midmatrix)*gsl_vector_max(temp),0.5);
}

void copyScale(double *shape,gsl_vector *temp3,float fdir) {

  int i;
  int n=temp3->size;

  for(i=0; i<n; i++)
    gsl_vector_set(temp3,i,shape[i]);

  gsl_vector_scale(temp3,fdir);
  //  setLineWidth(temp3,1.5);
}

float setSpectra(float **plane,int ind,int x1,int x2) {
  
  // remove noise and calculate a cutoff
  int i,j,h,top;
  int x=abs((x2+1)-x1);
  double data[ind*x];
  float mean,std,cut;

  for(h=0,i=0; i<ind; i++)
    for(j=x1; j<=x2; j++)
      data[h++]=plane[i][j];

  qsort(data,x*ind,sizeof(double),floatcompare);
  top=(int) x*ind*0.95;
  double data2[top];

  for(j=0; j<top; j++)
      data2[j]=data[j];

  mean=gsl_stats_mean(data2,1,top);
  std=gsl_stats_sd_m(data2,1,top,mean);

  cut=mean+2*std;
  //cut=0;

  for(i=0; i<ind; i++)
    for(j=x1; j<=x2; j++) 
      if(plane[i][j]<cut)
	  plane[i][j]=0;

  //  printf("mean: %f  std: %f cut: %f\n",mean,std,mean+2*std);
  return mean+2*std;
}

float lorentzian(int x,int x0,float h,float l) {

  return l*(h/(pow((x-x0),2)+h*h));
}

void setLineWidth(gsl_vector *shape,float h,float in) {

  int n=shape->size;
  int i,k,l,f;
  int center=n/2;
  float sum,norm,max;

  gsl_vector *lorentz=gsl_vector_alloc(n);
  gsl_vector *temp=gsl_vector_alloc(n);
  gsl_vector_memcpy(temp,shape);

  for(i=0; i<n; i++)
    gsl_vector_set(lorentz,i,lorentzian(i,center,h,in));
  
  max=gsl_vector_max(lorentz)*0.05;

  for(f=0,i=center; i<n; i++)
    if(gsl_vector_get(lorentz,i)>max)
      f++;
    else
      break;

  sum=0;
  for(k=center-f; k<=center+f; k++) 
      sum+=gsl_vector_get(lorentz,k);

  norm=1.0/sum;
  for(i=0; i<n; i++)
    gsl_vector_set(lorentz,i,gsl_vector_get(lorentz,i)*norm);

  for(i=0; i<n; i++) {
    sum=0;
    for(l=-f,k=-f+i; k<=f+i; k++,l++)
      if(k>0 && k<n)
	sum+=gsl_vector_get(temp,k)*gsl_vector_get(lorentz,center+l);
    gsl_vector_set(shape,i,sum);
   
  }

  //  writeTemp(shape,"T1");  
  gsl_vector_free(lorentz);
  gsl_vector_free(temp);
}

float pts2ppm(float x,header *hed,int nuc,int direct) {

  float offset=hed[nuc].O1_ppm-hed[nuc].SW_ppm/2.0;
  
  return (1.0-x/hed[nuc].size)*hed[nuc].SW_ppm+offset;
}

float ppm2pts(float ppm,header *hed,int nuc,int direct) {

  float offset=hed[nuc].O1_ppm-hed[nuc].SW_ppm/2.0;
  float obs=ppm-offset;
  
  if(direct)
    return hed[0].size-obs/hed[0].SW_ppm*hed[0].size;
  else
    return hed[nuc].size-obs/hed[nuc].SW_ppm*hed[nuc].size;
}

void reconstruct(interval interv,int exps,experiments *exp,gsl_vector *temp1,\
		 gsl_vector *temp2,gsl_matrix *slice,gsl_matrix *slice2,gsl_matrix *R,\
		 gsl_matrix *S,int comps,int ind,double ***shapes2,int shape,int r,\
		 float ***spec,gsl_matrix *fdir2,gsl_vector *cutoffs,gsl_matrix *cshapes) {

  int g,h,i,j,k,l;
  int x1=interv.x1;
  gsl_matrix *midma=gsl_matrix_alloc(ind,ind);
  gsl_vector *temp3=gsl_vector_alloc(ind);

  // S: reccorded spectra
  for(g=0,h=0; h<exps; h++) {
    for(k=0; k<r; k++)
      for(l=0; l<ind; l++) {
	gsl_matrix_set(S,l,k+(h*r),spec[h][l][k+x1]);
      }
      
    // R: reconstructed spectra
    if(exp[h].conv) {
      gsl_matrix_set_zero(slice2);
      gsl_vector_set_zero(temp1);
      gsl_vector_set_zero(temp2);
      //        #pragma omp parallel for
      for(i=0; i<comps; i++) {
	gsl_matrix_set_zero(slice);
	gsl_vector_view tt=gsl_matrix_row(cshapes,h*comps+i);
	gsl_matrix_view tshape=gsl_matrix_view_vector(&tt.vector,ind,1);
	//	convolve3(temp1,temp2,temp3,exp,shapes2,i,h);
	//	gsl_matrix_view tshape=gsl_matrix_view_vector(temp3,ind,1);
	gsl_matrix_view fdirt=gsl_matrix_submatrix(fdir2,i,0,1,r);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&tshape.matrix,&fdirt.matrix,0.0,slice);
	gsl_matrix_add(slice2,slice);
      }
    }
    else {
      gsl_matrix_set_zero(slice2);
      for(i=0; i<comps; i++) {
	gsl_vector_set_zero(temp1);
	gsl_matrix_set_zero(slice);
	for(l=0; l<ind; l++) 
	  gsl_vector_set(temp1,l,shapes2[i][exp[h].pos[0]][l]);
	gsl_matrix_view tshape=gsl_matrix_view_vector(temp1,ind,1);
	gsl_matrix_view tfdir=gsl_matrix_submatrix(fdir2,i,0,1,r);
	//	  gsl_matrix_fprintf(stdout,&tshape.matrix,"%f");
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&tshape.matrix,&tfdir.matrix,0.0,slice);
	gsl_matrix_add(slice2,slice);
      }
    }
    
    for(i=0; i<r; i++)
      for(j=0; j<ind; j++)
	gsl_matrix_set(R,j,i+(g*r),gsl_matrix_get(slice2,j,i));
    g++;
  }
  gsl_matrix_free(midma);
  gsl_vector_free(temp3);
}

void centerShift(float ***spec, int ind,int dir,int planes,int shift) {

  // shift spectra in order to simulate center at 0 position

  float *temp=malloc(sizeof(float)*ind);

  int h,i,j;
  
  for(h=0; h<planes; h++)
    for(i=0; i<dir; i++) {
      for(j=0; j<ind; j++) 
	temp[j]=spec[h][j][i];
      for(j=0; j<ind; j++)
	spec[h][j][i]=temp[(j+shift+ind)%ind];
    }

  free(temp);
}

unsigned modX(int a,int c) {
    int b=a%c;
    return b<0?b+c:b;
}

void lsSolver(gsl_matrix *MtM,gsl_vector *MtP,gsl_vector *X,int ind) {

  int s,j;

  gsl_permutation *np=gsl_permutation_alloc(ind);

  gsl_linalg_LU_decomp(MtM,np,&s);
  gsl_linalg_LU_solve (MtM,np,MtP,X);
  gsl_permutation_free(np);

  for(j=0; j<ind; j++)
    if(gsl_vector_get(X,j)<0)
      gsl_vector_set(X,j,0);
}

void regularization(gsl_vector *diag,int ind) {

  int i;

  for(i=0; i<ind; i++) {
    if(gsl_vector_get(diag,i)>0)
      gsl_vector_set(diag,i,1);
    if(gsl_vector_get(diag,i)==0)
      gsl_vector_set(diag,i,0);
    if(gsl_vector_get(diag,i)<0)
      gsl_vector_set(diag,i,-1);
  }
}

void subtract1(gsl_matrix *spec,gsl_matrix *shape,gsl_matrix *dir,int sub) {
  
  int i,j;
  int ind=shape->size1,r=dir->size2;
  float neg=0;

  gsl_matrix *recon=gsl_matrix_alloc(ind,r);
  
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,shape,dir,0.0,recon);
  gsl_matrix_sub(spec,recon);

  for(i=0; i<ind; i++)
    for(j=0; j<r; j++)
      if(gsl_matrix_get(spec,i,j)<0)
	//	gsl_matrix_set(spec,i,j,0);
	neg+=1;

  neg=100*neg/(ind*r);
  gsl_matrix_free(recon);
}

void writeMatrix(gsl_matrix *M,char *name) {

  FILE *file;
  file=fopen(name,"w");
  gsl_matrix_fprintf(file,M,"%f");
  fclose(file);
}

float gfactors1(double *gfactor,double ***shapes2,gsl_matrix *fdir2,int comp,\
	     int shapes,int shape,int ind) {

  int i,j;
  int rr=fdir2->size2;
  float sum=0,A=1;

  double *norms=malloc(sizeof(double)*(shapes+1));

  for(i=0; i<shapes; i++) {
    sum=0;
    gsl_vector_view st=gsl_vector_view_array(shapes2[comp][i],ind);
    for(j=0; j<ind; j++)
      sum+=gsl_vector_get(&st.vector,j);
    norms[i]=sum;
  }

  sum=0;
  for(i=0; i<rr; i++)
    sum+=gsl_matrix_get(fdir2,comp,i);
  
  norms[shapes]=sum;

  /* for(i=0; i<shapes+1; i++) */
  /*   printf("gfactor: %f\n",gfactor[i]); */

  for(i=0; i<shapes+1; i++)
    A*=norms[i];
  
  A=pow(A,1.0/(shapes+1));

  for(i=0; i<shapes+1; i++)
    gfactor[i]=A/norms[i];

  for(i=0; i<shapes+1; i++)
    printf("gfactor: %f\n",gfactor[i]);
  
  free(norms);

  return 0;
}

float gfactors2(double *gfactor,double ***shapes2,gsl_matrix *fdir2,int comp,\
	     int shapes,int shape,int ind) {

  int i;
  float A=1;

  double *norms=malloc(sizeof(double)*(shapes+1));

  for(i=0; i<shapes; i++) {
    gsl_vector_view st=gsl_vector_view_array(shapes2[comp][i],ind);
    norms[i]=gsl_blas_dnrm2(&st.vector);
  }

  gsl_vector_view st=gsl_matrix_row(fdir2,comp);
  norms[shapes]=gsl_blas_dnrm2(&st.vector);

  /* for(i=0; i<shapes+1; i++) */
  /*   printf("norms:   %f\n",norms[i]); */

  for(i=0; i<shapes+1; i++)
    A*=norms[i];
  
  A=pow(A,1.0/(shapes+1));

  for(i=0; i<shapes+1; i++)
    gfactor[i]=A/norms[i];

  for(i=0; i<shapes+1; i++)
    printf("gfactor: %f\n",gfactor[i]);
  
  free(norms);

  return 0;
}

float anormalization(double ***shapes,int nshapes,gsl_matrix *fdir,int comp,int ind) {

  int i;
  float f,a=1;

  for(i=0; i<nshapes; i++) {
    gsl_vector_view st=gsl_vector_view_array(shapes[comp][i],ind);
    f=gsl_blas_dnrm2(&st.vector);
    a*=f;
    gsl_vector_scale(&st.vector,1.0/f);
  }

  gsl_vector_view st=gsl_matrix_row(fdir,comp);
  f=gsl_blas_dnrm2(&st.vector);
  a*=f;
  gsl_vector_scale(&st.vector,1.0/f);

  return a;
}

int simplePeakpicker(double *shape,double *ppm,int nshape,int ind,peaklists *fpeaklist,float cut) {

  int i,j;
  int top=0;
  float max=0;
  float tres;
  double data[2*ind];
  double var,mean,std;

  for(i=0; i<ind; i++) {
    if(shape[i])
      data[top++]=shape[i];
    
    if(shape[i]>max)
      max=shape[i];
  }

  qsort(data,ind,sizeof(double),floatcompare);
  top=(int) ind*0.95;
  /* double data2[top]; */

  /* for(j=0; j<top; j++) */
  /*     data2[j]=data[j]; */
 
  mean=gsl_stats_mean(data,1,top);
  //  mean=val/top;
  std=gsl_stats_sd(data,1,top);
  /* std=gsl_stats_sd_with_fixed_mean(shapes,1,ind,mean); */
  var=gsl_stats_variance(data,1,top);

  std=std;
  var=var;
  mean=mean;
  tres=max*cut;

  for(i=0,j=0; i<ind; i++)
    if(i>0 && i<ind-1 && shape[i]>shape[i-1] && shape[i]>shape[i+1] && shape[i]>tres) {
      fpeaklist[j].shift=ppm[i];
      fpeaklist[j].point=i;
      fpeaklist[j++].intens=shape[i];
    }
  return j;
}

void peakpicker(double *shapes, double *ppm,int shape,int ind,int len1,int pc,\
		peaklists *fpeaklist,stats *stat,seqs *seq,header *hed,\
		int expect,int *pcount,char *residue,int minus) {

  int i,j,k,top;
  float cut=0,mean,std,max=0,minp,maxp;

  double data[ind];
  double tshape[ind];

  for(i=0; i<ind; i++) {
    data[i]=shapes[i];
    tshape[i]=shapes[i];
  }

  for(i=0; i<ind; i++)
    if(max<shapes[i])
      max=shapes[i];

  qsort(data,ind,sizeof(double),floatcompare);

  top=(int) ind*0.95;
  double data2[top];

  cut=cut; // to get rid of [-Wunused-but-set-variable]

  for(j=0; j<top; j++)
      data2[j]=data[j];

  mean=gsl_stats_mean(data2,1,ind);
  std=gsl_stats_sd_m(data2,1,ind,mean);
  cut=mean+2*std;

  for(j=0; j<len1; j++) {
    printf("%s %s %s %s\n",stat[j].res,residue,hed[shape].nuc,stat[j].nucs);
    if(strcmp(stat[j].res,residue)==0 && strncmp(hed[shape].nuc,stat[j].nucs,1)==0) {
      minp=stat[j].avg-3*stat[j].stddev;
      maxp=stat[j].avg+3*stat[j].stddev;
      for(k=0; k<ind; k++)
	if(0<k && k<ind-2 && tshape[k-1]<tshape[k] && tshape[k]>tshape[k+1] && \
	   tshape[k]>max*0.25 && ppm[k]>minp && ppm[k]<maxp) {
	  fpeaklist[*pcount].counter++;
	  fpeaklist[*pcount].shift=ppm[k];
	  fpeaklist[*pcount].zeros=0.0000;
	  fpeaklist[*pcount].nuc=stat[j].nucs;
	  fpeaklist[*pcount].resdn=pc;
	  (*pcount)++;
	}
    }
  }
}

int getFileLength(char *path) {

  FILE *file;
  char buffer[200];
  char *t;
  int i,count=0;

  file=fopen(path,"r");
  for(i=0; i<350; i++) {
    t=fgets(buffer,200,file);
    if(t!=NULL)
      count++;
  }
  fclose(file);
  return count;
}

void initStats(stats *stat,char *path) {

  FILE *file;
  //  char *url="http://www.bmrb.wisc.edu/ftp/pub/bmrb/statistics/chem_shifts/selected/statsel_prot.txt";
  int getText=0;
  char buffer[200];
  char *t,*token1;
  int i,j=0,k=0;

  // bmrb.txt
  if(getText)
    ; // get stats from BMRB
  else {
    // get from ascii file
    file=fopen(path,"r");
    for(i=0; i<350; i++) {
      t=fgets(buffer,200,file);
      if(i>13 && strlen(buffer)>10) { // hardcoded solution
	token1=strtok(buffer," ");
	while(token1!=NULL) {
	  if(k==0) 
	    strcpy(stat[j].res,token1);
	  if(k==1)
	    strcpy(stat[j].nucs,token1);
	  if(k==4)
	    stat[j].min=atof(token1);
	  if(k==5)
	    stat[j].max=atof(token1);
	  if(k==6)
	    stat[j].avg=atof(token1);
	  if(k==7)
	    stat[j].stddev=atof(token1);
 	  if(k>7) 
	    break;
	  k++;
	  //	  printf("%s\n",token1);
	  token1=strtok(NULL," ");
	}
	/* printf("%d %s %s %f %f %f %f\n",j,stat[j].res,stat[j].nucs, \ */
	/*        stat[j].min,stat[j].max,stat[j].stddev,\ */
	/*        stat[j].avg); */
	j++;
	k=0;
      }
    }
    t=t;
    fclose(file);
  }
}

void getSequence(seqs *seq,char *path) {

  int i;
  char buffer[100];
  FILE *file;
  file=fopen(path,"r");
  char *t,*token1;

  for(i=0; i<300; i++) {
    t=fgets(buffer,100,file);
    token1=strtok(buffer," ");
    if(t!=NULL)
      strcpy(seq[i].residue,token1);
  }

  /* for(i=0; i<45; i++) */
  /*   printf("%s\n",seq[i].res); */

  fclose(file);
}

int floatcompare(const void *a,const void *b) {
  
  double *x = (double *) a;
  double *y = (double *) b;

  if(*x<*y) 
    return -1;
  else if(*x>*y) 
    return 1; 
  return 0;
}

void svd2(gsl_matrix *MtM,gsl_vector *MtP,gsl_vector *X,int size1) {

  int i;

  gsl_matrix *V=gsl_matrix_alloc(size1,size1);
  gsl_vector *work=gsl_vector_alloc(size1);
  gsl_vector *S=gsl_vector_alloc(size1);

  gsl_linalg_SV_decomp(MtM,V,S,work);
  gsl_linalg_SV_solve (MtM,V,S,MtP,X);

  // remove zeros
  for(i=0; i<size1; i++)
    if(gsl_vector_get(X,i)<0)
      gsl_vector_set(X,i,0);

  gsl_matrix_free(V);
  gsl_vector_free(work);
  gsl_vector_free(S);

}

int selector(interval interv,int **peaklist) {

  int j,mindex=0,min=0;
  int xrange[6][2],allcomps[6];

  xrange[0][0]=interv.x1;
  xrange[0][1]=interv.x2;
  xrange[1][0]=interv.x1+2;
  xrange[1][1]=interv.x2+1;
  xrange[2][0]=interv.x1+3;
  xrange[2][1]=interv.x2+2;
  xrange[3][0]=interv.x1-1;
  xrange[3][1]=interv.x2-2;
  xrange[4][0]=interv.x1-2;
  xrange[4][1]=interv.x2-3;
  xrange[5][0]=interv.x1+1;
  xrange[5][1]=interv.x2-1;

  for(j=0; j<6; j++) {
    allcomps[j]=getComps(peaklist,xrange[j][0],xrange[j][1],1,0);
    if(min>allcomps[j]) {
      min=allcomps[j];
      mindex=j;
    }
  }
  
  return mindex;
}


int ncheck(int peak,double*** shapes2,int comps,int ind,float cut) {
  
  int i,c,j=0,maxi=0,b=1;
  float max;
  
  printf("peak: %d\n",peak);
  
  for(i=0; i<comps; i++) {
    gsl_vector_view st=gsl_vector_view_array(shapes2[i][0],ind); //hardcoded
    max=gsl_vector_max(&st.vector);
    max*=cut;
    //    maxi=gsl_vector_max_index(&st.vector);
    printf("shape matches: ");
    for(c=0,j=0; j<ind; j++) {
      if(j>0 && j<ind-1 &&						\
	 gsl_vector_get(&st.vector,j-1)<gsl_vector_get(&st.vector,j) && \
	 gsl_vector_get(&st.vector,j)>gsl_vector_get(&st.vector,j+1) && \
	 gsl_vector_get(&st.vector,j)>max) { 
	maxi=j;
	c++;
	printf("\tmatched peaks: %d %d %2.2f",c,maxi,gsl_vector_get(&st.vector,maxi));
      }
    }
    printf("\n");
    if(c==1 && abs(maxi-peak)<2)
      b=0;
  }
  return b;
}

int ncheck2(int peak,gsl_vector *nshape,int ind,float cut2,int comp,FILE *flog) {

  int j,maxi=0,c=0,b=1;
  float max;
  //  writeTemp(nshape,"C1");
  
  max=gsl_vector_max(nshape);  
  max*=cut2;

  for(j=0; j<ind; j++) 
      if(j>0 && j<ind-1 &&\
	 gsl_vector_get(nshape,j-1)<gsl_vector_get(nshape,j) && \
	 gsl_vector_get(nshape,j)>gsl_vector_get(nshape,j+1) && \
	 gsl_vector_get(nshape,j)>max) { 
	maxi=j;
	c++;
	printf("\tcomp: %d matched peak: %d %d %2.2f\n",comp,c,maxi,gsl_vector_get(nshape,maxi));
	fprintf(flog,"\tcomp: %d matched peak: %d %d %2.2f\n",comp,c,maxi,gsl_vector_get(nshape,maxi));
      }
  
  printf("\n");

  if(c==1 && abs(maxi-peak)<3) // hardcoded
      b=0;

  return b;
}
  

double cosinesim(double *a,double *b,int ind) {

  int i;
  double t=0,n1=0,n2=0;
  
  for(i=0; i<ind; i++) {
    t+=a[i]*b[i];
    n1+=a[i]*a[i];
    n2+=b[i]*b[i];
  }

  return t/(sqrt(n1)*sqrt(n2));
}

void checksim(gsl_vector *temp3,gsl_vector *temp2,int ind) {

  int i;
  double a,b; 

  for(i=0; i<ind; i++) {
    a=gsl_vector_get(temp3,i);
    b=gsl_vector_get(temp2,i);
    if(a>0 && b>0 && b>a)
      gsl_vector_set(temp3,i,b);
    if(a>0 && b==0)
      gsl_vector_set(temp3,i,0);
  }
}

int glycineDetection(double ***ppm,double **component,float cac,float hac,stats *stat, \
		     int N,int Ca,int Cb,int Ha,int Hb,int ind,int lenstat) {

  int i,nc=0,cc=0,hc=0,n1=0,n2=0;
  float minn=0,maxn=0,minca=0,maxca=0,minha=0,maxha=0,ncut,cacut,hacut,max;

  for(i=0; i<lenstat; i++) {
    if(strcmp(stat[i].res,"GLY")==0 && strcmp(stat[i].nucs,"N")==0) {
      minn=stat[i].avg-stat[i].stddev*2;
      maxn=stat[i].avg+stat[i].stddev*2;
    }
    if(strcmp(stat[i].res,"GLY")==0 && strcmp(stat[i].nucs,"CA")==0) {
      minca=stat[i].avg-stat[i].stddev*2;
      maxca=stat[i].avg+stat[i].stddev*2;
    }
    if(strcmp(stat[i].res,"GLY")==0 && strcmp(stat[i].nucs,"HA2")==0) {
      minha=stat[i].avg-stat[i].stddev*2; // should be an average of HA2 and HA3
      maxha=stat[i].avg+stat[i].stddev*2;
    }
  }
  
  max=0;
  for(i=0; i<ind; i++)
    if(max<component[N][i])
      max=max<component[N][i];
  
  ncut=max*0.2;

  max=0;
  for(i=0; i<ind; i++)
    if(max<component[Cb][i])
      max=max<component[Cb][i];
  
  cacut=max*0.2;

  max=0;
  for(i=0; i<ind; i++)
    if(max<component[Hb][i])
      max=max<component[Hb][i];
  
  hacut=max*0.2;

  for(i=1; i<ind-1; i++) {
    if(component[N][i]>ncut && component[N][i-1]<component[N][i] &&\
       component[N][i+1]<component[N][i] && ppm[0][N][i]>minn && ppm[0][N][i]<maxn)
      nc++;
    if(component[Cb][i]>cacut && component[Cb][i-1]<component[Cb][i] &&\
       component[Cb][i+1]<component[Cb][i] && ppm[0][Cb][i]>minca && ppm[0][Cb][i]<maxca)
      cc++;
    if(component[Hb][i]>hacut && component[Hb][i-1]<component[Hb][i] &&\
       component[Hb][i+1]<component[Hb][i] && ppm[0][Hb][i]>minha && ppm[0][Hb][i]<maxha)
      hc++;
    if(component[Ca][i]>cac)
      n1++;
    if(component[Ha][i]>hac)
      n2++;
  }
  printf("\t%d %d %d %d %d | %f %f %f %f %f %f | %f %f\n",nc,cc,hc,n1,n2,minn,maxn,minca,maxca,minha,maxha,cac,hac);
  if(nc==1 && cc<2 && hc<3 && n1==0 && n2==0)
    return 1;
  else
    return 0;
}

void saveMatrix(gsl_matrix *m,char* name,int csv) {

  int i,j;

  int size1=m->size1;
  int size2=m->size2;

  FILE *f;

  f=fopen(name,"w");

  if(csv)
    for(i=0; i<size1; i++) {
      for(j=0; j<size2; j++) 
	fprintf(f,"%f,",gsl_matrix_get(m,i,j));
      fprintf(f,"\n");
    }
  else
    for(i=0; i<size1; i++) {
      if(gsl_matrix_get(m,i,0)!=-1)
	for(j=0; j<size2; j++) {
	  if(gsl_matrix_get(m,i,j)!=-1)
	    fprintf(f,"%d ",(int) gsl_matrix_get(m,i,j));
	  else {
	    fprintf(f,"\n");
	    break;
	  }
	}
    }
  
  fclose(f);
}

int startingPoint(gsl_matrix *corrMatrix,int s) {

  int f,g,i,j;
  int size=corrMatrix->size1;
  
  for(f=0,g=0,i=s; i<size; i++) {
    for(j=0; j<size; j++) 
      if(gsl_matrix_get(corrMatrix,i,j)!=0)
	break;
    if(j==size)
      for(j=0; j<size; j++) {
	for(f=0; f<size; f++)
	  if(gsl_matrix_get(corrMatrix,f,j)!=0)
	    break;
	if(f==size) 
	  for(f=i+1; f<size; f++)  
	    for(g=0; g<size; g++)
	      if(gsl_matrix_get(corrMatrix,g,f)!=0)  
		return g;
      }
  }
  return -1;
}

float gausdist(float x,float a,float b,float c,float d) {

  // a: height, b: position, c: standard deviation, d: asymptot (0)

  return a*exp(-(pow(x-b,2)/(2*pow(c,2)))); 
}

void statShapes(gsl_matrix *statsh,stats *stat,int ind,int lenstat,header *hed,seqs *seq,int len2) {
  
  int h=0,i,j,k;
  float mean,std,temp,f,offset,sw,x;
  
  gsl_vector *temp2=gsl_vector_alloc(ind);

  for(i=0; i<len2; i++) {
    gsl_matrix_set(statsh,h++,0,i); // odd solution
    gsl_matrix_set(statsh,h++,1,999);
    for(j=0; j<lenstat; j++) {
      if(strcmp(stat[j].res,seq[i].residue)==0 && strcmp(stat[j].nucs,"N")==0) {
	//	printf("%s %s\n",stat[j].res,seq[i].residue);
	offset=hed[1].O1_ppm;
	sw=hed[1].SW_ppm;
	x=offset-sw/2.0;
	f=sw/(float) ind;
	mean=stat[j].avg;
	std=stat[j].stddev;
	for(k=0; k<ind; k++) {
	  x+=f;
	  temp=gausdist(x,1,mean,std,0);
	  gsl_matrix_set(statsh,h,k,temp);
	}
	//	mean=ppm2pts(mean,hed,1,0); // hardcoded
	for(k=0; k<ind; k++) {
	  temp=ppm2pts(gsl_matrix_get(statsh,h,k),hed,1,0);
	  gsl_matrix_set(statsh,h,k,temp);
	  gsl_vector_set(temp2,k,temp);
	}
      }
      if(strcmp(stat[j].res,seq[i].residue)==0 && strcmp(stat[j].nucs,"C")==0) {
	h++;
	mean=stat[j].avg;
	mean=ppm2pts(mean,hed,2,0); // hardcoded
	std=stat[j].stddev;
	for(k=0; k<ind; k++) {
	  gsl_matrix_set(statsh,h,k,gausdist(k,1,mean,std,0));
	}
      }
      if(strcmp(stat[j].res,seq[i].residue)==0 && strcmp(stat[j].nucs,"CA")==0) {
	h++;
	mean=stat[j].avg;
	mean=ppm2pts(mean,hed,5,0); // hardcoded
	std=stat[j].stddev;
	for(k=0; k<ind; k++) {
	    gsl_matrix_set(statsh,h,k,gausdist(k,1,mean,std,0));
	}
      }
      if(strcmp(stat[j].res,seq[i].residue)==0 && strcmp(stat[j].nucs,"CB")==0) {
	h++;
	mean=stat[j].avg;
	mean=ppm2pts(mean,hed,6,0); // hardcoded
	std=stat[j].stddev;
	for(k=0; k<ind; k++) {
	  gsl_matrix_set(statsh,h,k,gausdist(k,1,mean,std,0));
	}
      }
      if(strcmp(stat[j].res,seq[i].residue)==0 && strcmp(stat[j].nucs,"HA")==0) {
	h++;
	mean=stat[j].avg;
	mean=ppm2pts(mean,hed,7,0); // hardcoded
	std=stat[j].stddev;
	for(k=0; k<ind; k++) {
	  gsl_matrix_set(statsh,h,k,gausdist(k,1,mean,std,0));
	}
      }
      if(strcmp(stat[j].res,seq[i].residue)==0 && strcmp(stat[j].nucs,"HB")==0) {
	h++;
	mean=stat[j].avg;
	mean=ppm2pts(mean,hed,8,0); // hardcoded
	std=stat[j].stddev;
	for(k=0; k<ind; k++) {
	  gsl_matrix_set(statsh,h,k,gausdist(k,1,mean,std,0));
	}
      }
    }
  }

  gsl_vector_free(temp2);

}

int makeIntervalls1(int **peaklist,char *path) {
  // from nmrDraw peaklist

  FILE *p;
  char buffer[200];
  p=fopen(path,"r");
  char *t,*token;

  // get rid of gcc warnings
  (void) t;
  (void) token;

  int i,x,k=0;
  int len=20;

  for(i=0; i<len; i++)  
    if(i>17) { // hardcoded
      t=fgets(buffer,200,p);
      token=strtok(NULL," ");
      token=strtok(NULL," ");
      token=strtok(buffer," ");
      printf("buffer: %s\n",buffer);
      x=round(atof(buffer));
      peaklist[k][0]=x;
      token=strtok(buffer," ");
      x=round(atof(buffer));
      peaklist[k++][1]=x;
    } 
  fclose(p);

  return k;
}

int getIntervals2(char *path3) {

  // reads number of intervalls from either nmrDraw or xeasy

  FILE *ff;
  int i;
  int intervals=0;
  char *p,*token1,*ptr;
  char buffer[1000];

  ff=fopen(path3,"r");
  for(i=0; i<200; i++) { // hardcoded
    p=fgets(buffer,1000,ff);
    if(p==NULL)
      break;
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strtol(token1,&ptr,10)!=0) {
  	intervals++;
  	break;
      }
      else
  	break;
    }
  }
  p=p;
  fclose(ff);

  return intervals;
}

int checkPeakFormat(char *path3,gsl_vector *order) {

  char *token,*p;
  char buffer[100];
  char *name=malloc(100*sizeof(char));
  FILE *f;
  int i,format=0;

  //  printf("path3: %s\n",path3);
  f=fopen(path3,"r");
  if(f==NULL)
    exit(0);

  strcpy(name,path3);

  token=strtok(name,".");
  token=strtok(NULL," ");

  /* key:
     NH N
     [] [] */

  if(strcmp(token,"tab")==0) {
    gsl_vector_set(order,0,2);
    gsl_vector_set(order,1,1);
    format=1;
  }
  if(strcmp(token,"peaks")==0) {
    format=2;
    for(i=0; i<3; i++) {
      p=fgets(buffer,100,f);
      if(strstr(p,"#INAME") && strstr(p,"1") && (strstr(p,"15N") || strstr(p,"15N"))) {
	gsl_vector_set(order,0,2);
	gsl_vector_set(order,1,1);
	break;
      }
      if(strstr(p,"#INAME") && strstr(p,"1") && (strstr(p,"1H") || strstr(p,"NH"))) {
	gsl_vector_set(order,0,1);
	gsl_vector_set(order,1,2);
	break;
      }
    }
  }
  
  p=p; // supress compile errors
  fclose(f);
  free(name);
  return format;
}

int getIntervals3(char *path,int **peaklist,interval *interv,int intervals,header *hed,\
		  int format,FILE *flog,gsl_vector *order) {

  // get intervall points from either nmrDraw or neasy

  FILE *ff;
  int i,j,intpart,b=1,n=0,nh=0;
  char *p,*token1;
  char buffer[1000];
  float fpart=0.0;

  ff=fopen(path,"r");
  if(format==2) { // peak format
    for(i=0; i<3; i++) 
      p=fgets(buffer,1000,ff);
    nh=gsl_vector_get(order,0);
    n=gsl_vector_get(order,1);
    //    printf("nh: %d n: %d\n",nh,n);
    for(i=3,j=0; i<200; i++) { // hardcoded
      p=fgets(buffer,1000,ff);
      if(p==NULL)
	break;
      //      printf("p: %s\n",p);
      token1=strtok(p," ");
      interv[j].seq=atoi(token1);
      token1=strtok(NULL," ");
      //      printf("1 token1: %s\n",token1);
      if(nh==1) {
	fpart=atof(token1);
	interv[j].nh=fpart;
	//	fpart=ppm2pts(atof(token1),hed,0,1); // direct
	intpart=(int) roundf(ppm2pts(fpart,hed,0,1));
	intpart--; // for index start at 0
	intpart--; // for start x value in 3 points fdir
	interv[j].x1=intpart;
	interv[j].x2=intpart+2;
      }
      if(n==1) {
	fpart=atof(token1);
	interv[j].n=fpart;
	//	printf("fpart: %f\n",fpart);
	intpart=(int) roundf(ppm2pts(fpart,hed,1,0));
	//	printf("intpart: %d\n",intpart);
	intpart--;
	interv[j].y=intpart;
      }
      token1=strtok(NULL," ");
      //      printf("2 token1: %s\n",token1);
      if(nh==2) {
	fpart=atof(token1);
	interv[j].nh=fpart;
	//	printf("token1: %s\n",token1);
	intpart=(int) roundf(ppm2pts(fpart,hed,0,1));
	//	printf("intpart: %d\n",intpart);
	intpart--;
	intpart--;
	interv[j].x1=intpart;
	interv[j].x2=intpart+2;
      }
      if(n==2) {
	fpart=atof(token1);
	interv[j].n=fpart;
	intpart=(int) roundf(ppm2pts(fpart,hed,1,0));
	intpart--;
	interv[j].y=intpart;
      }
      peaklist[j][0]=interv[j].x1+b;
      peaklist[j][1]=interv[j].y;
      peaklist[j][2]=interv[j].seq;
      printf("interval: %d %d %d %d %f %f\n",j,peaklist[j][2],peaklist[j][0],peaklist[j][1],\
	     interv[j].nh,interv[j].n);
      fprintf(flog,"interval: %d %d %d %d %f %f\n",j,peaklist[j][2],peaklist[j][0],peaklist[j][1], \
	     interv[j].nh,interv[j].n);
      j++;
    }
  }

  p=p;
  fclose(ff);

  return intervals;
}

float getDir(gsl_matrix *fdir2,double **ppmd,int x1,int r,int b) {

  int i,c=0;
  float max=0;

  for(i=0; i<r; i++)
    if(gsl_matrix_get(fdir2,0,i)>max) {
      max=gsl_matrix_get(fdir2,0,i);
      c=i;
    }

  if(b)
    return x1+c;
  else
    return ppmd[0][c];
}

void getCorrmat(gsl_matrix *corrmat,component *components,int intervals,int nshapes,int ind) {

  int i=0,n,j,k=0;
  int i1,f2,counta,countb;
  float f1,f3;
  char s1[10],buffer[100];
  char *s2;
  FILE *fp;
  float *cali=NULL,*cai=NULL,*cbi=NULL;

  /* float *cam=calloc(intervals,sizeof(float)); */
  /* float *cbm=calloc(intervals,sizeof(float)); */
  /* float *caim=calloc(intervals,sizeof(float)); */
  /* float *cbim=calloc(intervals,sizeof(float)); */

  gsl_vector *cam=gsl_vector_calloc(intervals);
  gsl_vector *cbm=gsl_vector_calloc(intervals);
  gsl_vector *caim=gsl_vector_calloc(intervals);
  gsl_vector *cbim=gsl_vector_calloc(intervals);

  fp=fopen("results/latestPeakList.txt","r");

  printf("getCorrmat\n");
  while (fgets(buffer,100,fp)!=NULL) {
    sscanf(buffer,"%i%s%f%d%f",&i1,s1,&f1,&f2,&f3);
    n=strlen(s1);
    s2=(char *) calloc(n-1,sizeof(char));
    strncpy(s2,s1,n-1);
    for(j=0; j<nshapes; j++) {
      if(strcmp(s2,components[i1].nuclei[j])==0) {
	if(components[i1].shifts[j][0]==0)
	  k=0;
	components[i1].shifts[j][k]=f1;
	components[i1].points[j][k]=f2;
	components[i1].intens[j][k]=f3;
	printf("%s %s %f %d %f %d\n",s2,components[i1].nuclei[j],components[i1].shifts[j][k],\
	       components[i1].points[j][k],components[i1].intens[j][k],k);
	k++;
      }
    }
    free(s2);

    cali=components[i1].shifts[components[i1].compare[0]]; // hardcoded
    cai=components[i1].shifts[components[i1].compare[1]];
    cbi=components[i1].shifts[components[i1].compare[2]];

    counta=0;
    countb=0;

    for(i=0; i<ind; i++)
      if(cai[i])
	printf("controll 1: %f\n",cai[i]);

    for(j=0; j<ind; j++) {
      if(cai[j])
	for(k=0; k<ind; k++) 
	  if(cali[k] && fabs(cai[j]-cali[k])<0.5) // hardcoded
	    counta++;
	  //	    printf("cali: %f ca(-1) %f %f %d\n",cali[k],cai[j],fabs(cai[j]-cali[k]),j);
      if(counta==0)
	cai[j]=0;	
      if(cbi[j])
	for(k=0; k<ind; k++) 
	  if(cali[k] && fabs(cbi[j]-cali[k])<0.5) // hardcoded
	    //	    printf("cali: %f cb(-1) %f %f\n",cali[k],cbi[j],fabs(cbi[j]-cali[k]));
	    countb++;
      if(countb==0)
	cbi[j]=0;	
    }

    for(i=0; i<ind; i++)
      if(cai[i])
	printf("controll 2: %f\n",cai[i]);
    P(counta)

    if(counta==1)
      gsl_vector_set(caim,i1,cai[0]);
    if(countb==1)
      gsl_vector_set(cbim,i1,cbi[0]);
  }

  fclose(fp);

  fp=fopen("caim","w");
  gsl_vector_fprintf(fp,caim,"%f");
  fclose(fp);

  fp=fopen("cbim","w");
  gsl_vector_fprintf(fp,cbim,"%f");
  fclose(fp);

  // select cali(i-1) and i
  /* for(i=0; i<1; i++) { */
  /*   cali=components[i].shifts[components[i].compare[0]]; */
  /*   cai=components[i].shifts[components[i].compare[1]]; */
  /*   cbi=components[i].shifts[components[i].compare[2]]; */
  /* } */

  /* for(i=0; i<ind; i++) */
  /*   if(cali[i]) */
  /*     printf("cali: %f\n",cali[i]); */

  /* for(i=0; i<ind; i++) */
  /*   if(cai[i]) */
  /*     printf("ca: %f\n",cai[i]); */

  /* for(i=0; i<ind; i++) */
  /*   if(cbi[i]) */
  /*     printf("cb: %f\n",cbi[i]); */

  /* // detect i */
  /* for(i=0; i<1; i++) */
  /*   for(j=0; j<ind; j++) { */
  /*     if(cai[j]) */
  /* 	for(k=0; k<ind; k++) { */
  /* 	  if(cali[k] && fabs(cai[j]-cali[k])<0.5) // hardcoded */
  /* 	    printf("cali: %f ca(-1) %f %f %d\n",cali[k],cai[j],fabs(cai[j]-cali[k]),j); */
  /* 	  else */
  /* 	    cai[j]=0; */
  /* 	} */
  /*     if(cbi[j]) */
  /* 	for(k=0; k<ind; k++) { */
  /* 	  if(cali[k] && fabs(cbi[j]-cali[k])<0.5) // hardcoded */
  /* 	    printf("cali: %f cb(-1) %f %f\n",cali[k],cbi[j],fabs(cbi[j]-cali[k])); */
  /* 	  else */
  /* 	    cbi[j]=0; */
  /* 	} */
  /*   } */

  /* counta=0; */
  /* countb=0; */

  /* for(i=0; i<ind; i++) { */
  /*   if(cai[i]) */
  /*     counta++; */
  /*   if(cbi[i]) */
  /*     countb++; */
  /* } */
  /* for(i=0; i<intervals; i++) */
  /*   for(j=0; j<intervals; j++) */
  /*     ; */

  gsl_vector_free(cam); 
  gsl_vector_free(cbm); 
  gsl_vector_free(caim); 
  gsl_vector_free(cbim); 
}

/* void convCheck(gsl_vector *v1,gsl_vector *v2,gsl_vector *v3,int ind) { */

/*   printf("") */
  
/* } */

/* int tiltCheck(experiments exp,header hed,char **nuclei,int exps) { */

/*   int i; */

/*   for(i=0; i<exps; i++) */
/*     if(strchr(exp[i].sh,43)) { // 43 = + */
      

/*   return 1; */
/* } */
