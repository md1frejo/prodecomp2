/* common prototypes */

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <math.h>

typedef struct {

  char *sh;
  int *pos;
  int conv;
  int *sign;

} definitions;

typedef struct {

  int x1,x2,y,dir,conv,sel,comps,seq;
  float nh,n;
  char residue[20];
  gsl_vector *f;
  gsl_vector *fdir;

} interval;

typedef struct {

  char *format;
  char *nuc;
  float SW_ppm;
  float SW_hz;
  float O1_ppm;
  int size;

} header;

typedef struct {

  char *sh;
  int *pos;
  int *sign;
  char *path;
  int *defs;
  int conv;
  int t;
  int tilt;
} experiments;

typedef struct {

  char *nucs;
  float min;
  float max;
  float avg;
  float stddev;

} nstat;

typedef struct {

  char *res;
  char *nucs;
  float min;
  float max;
  float avg;
  float stddev;

} stats;

typedef struct {

  int counter;
  float dshift;
  int dpoint;
  float shift;
  int point;
  float intens;
  float zeros;
  char *nuc;
  int resdn;
  int seq;

} peaklists;
  
typedef struct {

  char *residue;
  int component;

} seqs;

typedef struct { 

  int comp;
  int shapes;
  char **nuclei;
  float **shifts;
  int **points;
  float **intens;
  int *compare;

} component;

void setRandom1(gsl_matrix*,int);
int setRandom2(int);

float rmsd1(gsl_matrix *,gsl_matrix *,int,float);
float rmsd2(gsl_matrix *,gsl_matrix *,int,int);

void readHeader(char *, header *,int,int);

void midmatrix3(gsl_matrix *,gsl_vector *,int,int);

void ppm2pts1(header *,char *,int **,int);
int getPeaks(int**,int **,int,int,int,int,int,int,int);
int getPeaks2(int**,int,int,int,int,int,int,int);
int getComps(int **,int,int,int,int);
void readNMRDraw(int**,char *,int);

void readProdecomptext(char *,experiments *,int,char**,int);
int intcompare(const void *,const void *);
double readFT2(float ***,experiments *,int,int,int,int);
int getIntervals(char *);
int getIntervals2(char *);
int getnshapes(char *);
int getnexperiments(char *);
void getInterv2(char *,interval *);

void getDefs(definitions *); 
void getInterv(interval *); 
void getProj(gsl_matrix *);
void getSet(int **);
void getMatrix(gsl_matrix *,experiments,int,int);
int getPlane(experiments *,int,int);
int getPlanes(experiments *,int,int,int);

void writeTemp(gsl_vector *,char *);
void writeShape(int,int,int,double ***,char *);
float minmaxShape(double *,int,int); 
int minmaxIndexShape(double *,int,int);
void writeTrace(gsl_vector *,int,int,char *); 
float maxTrace(gsl_vector *,int,int);
void writeSpec(float ***,int,int,int,char *); 
float avgTrace(gsl_vector *,int,int);
void convolve1(gsl_vector *,gsl_vector *,experiments *,gsl_matrix *,double **,int,int,int);
void convolve2(gsl_vector *,gsl_vector *,experiments *,gsl_matrix *,double ***,int,int);
void convolve3(gsl_vector *,gsl_vector *,gsl_vector *,experiments *,double **,int);
void convolve4a(gsl_vector *,gsl_vector *,experiments *,gsl_matrix *,double **,int,int,int,int);
void convolve4b(gsl_vector *,gsl_vector *,experiments *,double **,int,int,int);
void midmatrix4(gsl_matrix *,gsl_vector *,int,int);
void midmatrix5(gsl_matrix *,gsl_vector *,int,int,int);
void midmatrix6(gsl_matrix *,gsl_vector *,int,int,int,int);
void rconvolution(gsl_vector *,gsl_vector *,gsl_vector *,int,int);
void convolution(gsl_vector *,gsl_vector *,gsl_vector *,int,int,int);
void convolution2(gsl_vector *,gsl_vector *,gsl_vector *,int,int,int,int);
float kfactor(gsl_matrix *,gsl_vector *,gsl_vector *,float);
float subFactor1(gsl_vector *, float ***,int,int,int,float);
float subFactor2(gsl_vector *, gsl_vector *,float);
void normalizeTemp(gsl_vector *);
void normalizeShape(float *** ,int,int,int);
float nCutoff(gsl_vector *);
void lineBroadening(gsl_vector *,int,float);
float normFactor(gsl_matrix*,gsl_vector*);
void copyScale(double *,gsl_vector *,float);
float setSpectra(float **,int,int,int);
float lorentzian(int,int,float,float);
void setLineWidth(gsl_vector *,float,float);
float pts2ppm(float,header *,int,int);
float ppm2pts(float,header *,int,int);
void reconstruct(interval,int,experiments *,gsl_vector *,\
		 gsl_vector *,gsl_matrix *,gsl_matrix *,\
		 gsl_matrix *,gsl_matrix *,int,int,double ***,\
		 int,int,float ***,gsl_matrix *,gsl_vector *,\
		 gsl_matrix *);

void centerShift(float ***,int,int,int,int);
unsigned modX(int,int);
void lsSolver(gsl_matrix *,gsl_vector *,gsl_vector *,int);
void regularization(gsl_vector *,int);
void subtract1(gsl_matrix *,gsl_matrix *,gsl_matrix *,int);
void writeMatrix(gsl_matrix *,char *);
float gfactors1(double *,double ***,gsl_matrix *,int,int,int,int);
float gfactors2(double *,double ***,gsl_matrix *,int,int,int,int);
float anormalization(double ***,int,gsl_matrix *,int,int);
void normalizeShape2(gsl_vector *,int,int);
int simplePeakpicker(double *,double *,int,int,peaklists *,float);
void peakpicker(double *,double *,int,int,int,int,peaklists *,stats *,\
		seqs *,header *,int,int *,char *,int);
int getFileLength(char *);
void initStats(stats *,char *);
void getSequence(seqs *,char *);
int floatcompare(const void *,const void *);
void svd2(gsl_matrix *,gsl_vector *,gsl_vector *,int);
int selector(interval,int **);
int ncheck(int,double ***,int,int,float);
int ncheck2(int,gsl_vector *,int,float,int,FILE *);
double cosinesim(double *,double *,int);
void checksim(gsl_vector *,gsl_vector *,int);
int glycineDetection(double ***,double **,float,float,stats *,int,int,int,int,int,int,int);

int decompose(float ***,interval,\
		int,int,experiments*,int,int**,int,\
	       int,gsl_vector *,double ***,double ***,gsl_matrix *,gsl_matrix *, \
	       gsl_vector *,int **,gsl_matrix *,int,gsl_vector *);

void correlation(gsl_matrix *,double ***,gsl_matrix *,int,seqs *,float);
void saveMatrix(gsl_matrix *,char *,int);
int startingPoint(gsl_matrix *,int);
float gausdist(float,float,float,float,float);
void statShapes(gsl_matrix *,stats *,int,int,header *,seqs *,int);
int makeIntervalls1(int **,char *);
int checkPeakFormat(char *,gsl_vector *);
int getIntervals3(char *,int **,interval *,int,header *,int,FILE *,gsl_vector *);
float getDir(gsl_matrix *,double **,int,int,int);
void midmatrix7(gsl_matrix *,gsl_vector *,int,int,int);
void getCorrmat(gsl_matrix *,component *,int,int,int);
