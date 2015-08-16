/* parse prodecomp.txt and input.txt files */

#include "headerfiles/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void readHeader(char *path,header *hed,int shape,int direct) {

  FILE *f;
  char buffer[1000],b[100];
  int i=-1;
  char *token1;
  char *t;

  if(direct) {
    // direct dim
    f=fopen(path,"r");
    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"FORMAT=")!=0) {
	strcpy(hed[shape].format,token1);
	break;
      }
      else
	token1=strtok(NULL," ");
    }

    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"NUCLEI=")!=0) {
	strncpy(b,token1,strlen(token1)-1);
	b[strlen(token1)-1]='\0';
	strcpy(hed[shape].nuc,b);
	break;
      }
      else
	token1=strtok(NULL," ");
    }

    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"SW_ppm=")!=0) {
	strncpy(b,token1,strlen(token1)-2);
	hed[shape].SW_ppm=atof(b+1);
	break;
      }
      else
	token1=strtok(NULL," ");
    }

    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"SW_hz=")!=0) {
	strncpy(b,token1,strlen(token1)-2);
	b[strlen(token1)-2]='\0';
	hed[shape].SW_hz=atof(b+1);
	break;
      }
      else
	token1=strtok(NULL," ");
    }
    
    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"O1_ppm=")!=0) {
	strncpy(b,token1,strlen(token1)-2);
	b[strlen(token1)-2]='\0';
	hed[shape].O1_ppm=atof(b+1);
	break;
      }
      else
	token1=strtok(NULL," ");
    }

    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"SIZE=")!=0) {
	strncpy(b,token1,strlen(token1)-2);
	b[strlen(token1)-2]='\0';
	hed[shape].size=atoi(b+1);
	break;
      }
      else
	token1=strtok(NULL," ");
    }
    fclose(f);
  }

  // indirect dim

  else {
    f=fopen(path,"r");
    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"FORMAT=")!=0) {
	strncpy(b,token1,strlen(token1));
	strcpy(hed[shape].format,token1);
	break;
      }
      else
	token1=strtok(NULL," ");
    }
    
    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"NUCLEI=")!=0 && i==shape) {
	strncpy(b,token1,strlen(token1)-1);
	b[strlen(token1)-1]='\0';
	strcpy(hed[shape].nuc,b);
	break;
      }
      else
	token1=strtok(NULL," ");
      i++;
    }
    i=-1;
    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"SW_ppm=")!=0 && i==shape) {
	strcpy(b,token1);
	b[strlen(token1)-1]='\0';
	hed[shape].SW_ppm=atof(b);
	break;
      }
      else
	token1=strtok(NULL," ");
      i++;
    }

    i=-1;
    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"SW_hz=")!=0 && i==shape) {
	strcpy(b,token1);
	b[strlen(token1)-1]='\0';
	hed[shape].SW_hz=atof(b);
	break;
      }
      else
	token1=strtok(NULL," ");
      i++;
    }

    i=-1;
    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"O1_ppm=")!=0 && i==shape) {
	strcpy(b,token1);
	b[strlen(token1)-1]='\0';
	hed[shape].O1_ppm=atof(b);
	break;
      }
      else
	token1=strtok(NULL," ");
      i++;
    }

    i=0;
    t=fgets(buffer,200,f);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"SIZE=")!=0 && i==2) {
	strcpy(b,token1);
	b[strlen(token1)-1]='\0';
	hed[shape].size=atoi(b);
	break;
      }
      else
	token1=strtok(NULL," ");
      i++;
    }
    fclose(f);
  }
  t=t;
}

void readProdecomptext(char *path,experiments *exp,int nucs,char **nuclei,int exps) {
  
  FILE *f;
  char buffer[1000],buf[40];
  char *t,*token1;
  int i,j,k;
  int s;

  j=0;
  // paths first
  f=fopen(path,"r");
  for(i=0; i<200; i++) { // hardcoded
    t=fgets(buffer,300,f);
    if(i>5) {
      token1=strtok(buffer," ");
      while(token1!=NULL && strcmp(token1,"#")!=0) {
	if(strcmp(token1,"PATH:")==0) {
	  token1=strtok(NULL," ");
	  strcpy(exp[j].path,token1);
	  j++;
	  break;
	}
	else
	  token1=strtok(NULL," ");
      }
    }
    t=t;
  }
  fclose(f);

  j=0;
  // then definitions
  f=fopen(path,"r");
  for(i=0; i<200; i++) {  // hardcoded
    t=fgets(buffer,1000,f);
    if(i>5) {
      token1=strtok(buffer," ");
      while(token1!=NULL && strcmp(token1,"#")!=0) {
	//	printf("2a token1: %s %zu %d %d\n",token1,strlen(token1),j,exps);
	if(strcmp(token1,"DEFINITION:")==0 && j<exps) { // why exps?
	  //	  printf("2b token1: %s %d\n",token1,nucs);
	  for(k=0; k<nucs; k++) {
	    token1=strtok(NULL,",");
	    if(k==0) {
	      printf("3a token1 %s %d\n",token1,k);
 	      strcpy(buf,"");
	    }
	    else { 
	      //	      strcat(buf,nuclei[k]);
	      //	      printf("3b token1 %s %d\n",token1,k);
	      s=atoi(token1);
	      if(s) { 
		exp[j].defs[k]=s;
		if(s<0)
		  strcat(buf,"-");
		if(s>0)
		  strcat(buf,"+");
		strcat(buf,nuclei[k]);
	      }
	    }
	    //	    printf("defs: %d\n",exp[j].defs[k]);
	    if(token1==NULL)
	      break;
	  }
	  strcpy(exp[j].sh,buf);
	  j++;
	}
	else 
	  token1=strtok(NULL," ");
      }
    }
  }
  fclose(f);
}

int getIntervals(char *path) {

  FILE *ff;
  int i;
  int intervals=0;
  char *p,*token1;
  char buffer[1000];

  ff=fopen(path,"r");
  for(i=0; i<200; i++) {
    p=fgets(buffer,1000,ff);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"PATH:")==0) {
  	intervals++;
  	break;
      }
      else
  	token1=strtok(NULL," ");
    }
  }
  p=p;
  fclose(ff);

  return intervals;
}

int getnshapes(char *path) {

  // only indirect dim

  FILE *ff;
  int i;
  int nshapes=0;
  char *p,*token1;
  char buffer[1000];

  ff=fopen(path,"r");
  for(i=0; i<200; i++) {
    p=fgets(buffer,1000,ff);
    token1=strtok(buffer," ");
    while(token1!=NULL) {
      if(strcmp(token1,"NUCLEI=")==0) {
	token1=strtok(NULL," "); 
	while(token1!=NULL) {
	  //	  printf("token1: %s\n",token1);
	  nshapes++;
	  token1=strtok(NULL," ");
	}
      }
      else
  	token1=strtok(NULL," ");
    }
  }

  p=p;
  fclose(ff);

  return nshapes;
}

int getnexperiments(char *path) {

  FILE *ff;
  int i;
  int exps=0;
  char *p,*token1;
  char buffer[1000];
  
  ff=fopen(path,"r");
  printf("%s\n",path);
  for(i=0; i<200; i++) { // hardcoded
    p=fgets(buffer,1000,ff);
    token1=strtok(buffer," ");
    while(token1!=NULL && strcmp(token1,"#")!=0) {
      if(strcmp(token1,"PATH:")==0) {
  	exps++;
  	break;
      }
      else
  	token1=strtok(NULL," ");
    }
  }
  p=p;
  fclose(ff);
  
  return exps;
}

void getInterv2(char *path,interval *interv) {

  FILE *f;  
  char buffer[200];
  char *t,*token1;
  int i,j,x1,x2,ct;

  f=fopen(path,"r");

  j=0;
  for(i=0; i<100; i++) {
    t=fgets(buffer,200,f);
    if(i>3) { // hardcoded
      token1=strtok(buffer,"/");
      //      printf("1: token1: %s\n",token1);
      while(token1!=NULL) {
	if(strcmp(token1,"ProdecompOutput")==0) {
	  token1=strtok(NULL,"/");
	  x1=atoi(strtok(token1,"_"));
	  token1=strtok(NULL,"_");
	  x2=atoi(strtok(token1,"c"));
	  token1=strtok(NULL,"c");
	  ct=atoi(strtok(token1,"c"));
	  interv[j].x1=x1;
	  interv[j].x2=x2;
	  interv[j].comps=ct; // unecessary?
	  j++;
	  break;
	}
	else
	token1=strtok(NULL,"/");
      }
    }
    t=t;
  }
  fclose(f);
}
