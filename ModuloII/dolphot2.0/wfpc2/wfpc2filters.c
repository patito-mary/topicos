#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

typedef struct {
   char name[11],color;
   double zp[4];
   char xorder[8];
   int xformc[8];
   double xform[8][5];
   //double idc[2][3][33];
} WFPC2filttype;

int WFPC2_NFILTERS=-1;
WFPC2filttype *WFPC2filters;

#define NFILTERS_DEFINED 38

// New CTE Parameters (cold only);
//background "contamination";
double bgcorr[NFILTERS_DEFINED][2]={
   {0.000000e-06,2.869069e-04}, //F122W
   {0.000000e-06,2.869069e-04}, //F160BW
   {0.000000e-06,2.869069e-04}, //F170W
   {0.000000e-06,2.869069e-04}, //F185W
   {0.000000e-06,2.869069e-04}, //F218W
   {0.000000e-06,2.869069e-04}, //F255W
   {8.540038e-07,1.010904e-04}, //F300W
   {8.540038e-07,1.171397e-04}, //F336W
   {8.540038e-07,1.171397e-04}, //F343N
   {8.540038e-07,9.556493e-05}, //F375N
   {8.540038e-07,9.556493e-05}, //F380W
   {8.540038e-07,9.556493e-05}, //F390N
   {1.085865e-05,3.178520e-05}, //F410M
   {1.085865e-05,3.178520e-05}, //F437N
   {1.085865e-05,3.178520e-05}, //F439W
   {1.085865e-05,3.178520e-05}, //F450W
   {1.085865e-05,3.178520e-05}, //F467M
   {1.085865e-05,3.178520e-05}, //F469N
   {1.085865e-05,3.178520e-05}, //F487N
   {1.077151e-05,3.013676e-05}, //F502N
   {1.077151e-05,3.013676e-05}, //F547M
   {1.077151e-05,3.013676e-05}, //F555W
   {1.077151e-05,3.013676e-05}, //F569W
   {1.077151e-05,3.013676e-05}, //F588N
   {1.077151e-05,3.013676e-05}, //F606W
   {7.776181e-06,2.758164e-05}, //F622W
   {7.776181e-06,2.758164e-05}, //F631N
   {7.776181e-06,2.758164e-05}, //F656N
   {7.776181e-06,2.758164e-05}, //F658N
   {7.776181e-06,2.758164e-05}, //F673N
   {7.776181e-06,2.758164e-05}, //F675W
   {7.776181e-06,2.758164e-05}, //F702W
   {7.917809e-06,2.726974e-05}, //F785LP
   {7.917809e-06,2.726974e-05}, //F791W
   {7.917809e-06,2.726974e-05}, //F814W
   {7.917809e-06,2.726974e-05}, //F850LP
   {7.917809e-06,2.726974e-05}, //F953N
   {7.917809e-06,2.726974e-05}  //F1042M
};

//zero point corrections;
double fcorr[NFILTERS_DEFINED][2][4]={
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F122M
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F160BW
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F170W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F185W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F218W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F255W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F300W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F336W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F343N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F375N
   {{-0.142423,-0.144297,-0.137441,-0.150253},{-0.122420,-0.141541,-0.118225,-0.140175}}, //F380W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F390N
   {{-0.116887,-0.118761,-0.111904,-0.124717},{-0.096884,-0.116004,-0.092689,-0.114639}}, //F410M
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F437N
   {{-0.065428,-0.067302,-0.060445,-0.073258},{-0.045425,-0.064546,-0.041230,-0.063180}}, //F439W
   {{-0.058529,-0.060403,-0.053547,-0.066359},{-0.038526,-0.057647,-0.034331,-0.056281}}, //F450W
   {{0.003011,0.001137,0.007994,-0.004818},{0.023014,0.003894,0.027210,0.005260}}, //F467M
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F469N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F487N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F502N
   {{-0.040534,-0.039352,-0.035151,-0.043770},{-0.035417,-0.029731,-0.031918,-0.046449}}, //F547M
   {{-0.044708,-0.043526,-0.039325,-0.047944},{-0.039591,-0.033905,-0.036092,-0.050623}}, //F555W
   {{-0.023414,-0.022232,-0.018031,-0.026651},{-0.018297,-0.012612,-0.014799,-0.029329}}, //F569W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F588N
   {{-0.027214,-0.026032,-0.021831,-0.030451},{-0.022097,-0.016412,-0.018599,-0.033129}}, //F606W
   {{-0.047777,-0.048085,-0.039618,-0.057412},{-0.028274,-0.041024,-0.033258,-0.033968}}, //F622W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F631N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F656N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F658N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F673N
   {{-0.060353,-0.060661,-0.052194,-0.069988},{-0.040850,-0.053600,-0.045834,-0.046545}}, //F675W
   {{-0.065856,-0.066164,-0.057696,-0.075490},{-0.046352,-0.059102,-0.051337,-0.052047}}, //F702W
   {{-0.020142,-0.042453,-0.036936,-0.038852},{-0.031887,-0.042999,-0.031070,-0.036445}}, //F785LP
   {{0.004806,-0.017505,-0.011987,-0.013904},{-0.006939,-0.018051,-0.006122,-0.011497}}, //F791W
   {{-0.013481,-0.035791,-0.030274,-0.032190},{-0.025225,-0.036338,-0.024409,-0.029783}}, //F814W
   {{-0.021065,-0.043375,-0.037858,-0.039774},{-0.032809,-0.043921,-0.031992,-0.037367}}, //F850LP
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F953N
   {{-0.084016,-0.106326,-0.100809,-0.102725},{-0.095760,-0.106872,-0.094943,-0.100318}} //F1042M
};

// Old CTE Parameters (cold and warm) from Dolphin 2000
//background "contamination";
double bgcorr_warm[NFILTERS_DEFINED][2][2]={
   {{0.000000e-06,2.869069e-04},{0.000000e-06,2.292851e-04}}, //F122W
   {{0.000000e-06,2.869069e-04},{0.000000e-06,2.292851e-04}}, //F160W
   {{0.000000e-06,2.869069e-04},{0.000000e-06,2.292851e-04}}, //F170W
   {{0.000000e-06,2.869069e-04},{0.000000e-06,2.292851e-04}}, //F185W
   {{0.000000e-06,2.869069e-04},{0.000000e-06,2.292851e-04}}, //F218W
   {{0.000000e-06,2.869069e-04},{0.000000e-06,2.292851e-04}}, //F255W
   {{8.540038e-07,1.010904e-04},{2.347928e-06,8.078760e-05}}, //F300W
   {{8.540038e-07,1.171397e-04},{2.347928e-06,1.000053e-04}}, //F336W
   {{8.540038e-07,1.171397e-04},{2.347928e-06,1.000053e-04}}, //F343N
   {{8.540038e-07,9.556493e-05},{2.347928e-06,1.000053e-04}}, //F375N
   {{8.540038e-07,9.556493e-05},{2.347928e-06,1.000053e-04}}, //F380W
   {{8.540038e-07,9.556493e-05},{2.347928e-06,1.000053e-04}}, //F390N
   {{1.085865e-05,3.178520e-05},{0.000000e-06,5.432855e-05}}, //F410M
   {{1.085865e-05,3.178520e-05},{0.000000e-06,5.432855e-05}}, //F437N
   {{1.085865e-05,3.178520e-05},{0.000000e-06,5.432855e-05}}, //F439W
   {{1.085865e-05,3.178520e-05},{0.000000e-06,5.432855e-05}}, //F450W
   {{1.085865e-05,3.178520e-05},{0.000000e-06,5.432855e-05}}, //F467M
   {{1.085865e-05,3.178520e-05},{0.000000e-06,5.432855e-05}}, //F469N
   {{1.085865e-05,3.178520e-05},{0.000000e-06,5.432855e-05}}, //F487N
   {{1.077151e-05,3.013676e-05},{5.497646e-06,8.393854e-05}}, //F502N
   {{1.077151e-05,3.013676e-05},{5.497646e-06,8.393854e-05}}, //F547M
   {{1.077151e-05,3.013676e-05},{5.497646e-06,7.363844e-05}}, //F555W
   {{1.077151e-05,3.013676e-05},{5.497646e-06,7.331492e-05}}, //F569W
   {{1.077151e-05,3.013676e-05},{5.497646e-06,7.331492e-05}}, //F588N
   {{1.077151e-05,3.013676e-05},{5.497646e-06,6.200388e-05}}, //F606W
   {{7.776181e-06,2.758164e-05},{1.270904e-06,1.109481e-04}}, //F622W
   {{7.776181e-06,2.758164e-05},{1.270904e-06,1.109481e-04}}, //F631N
   {{7.776181e-06,2.758164e-05},{1.270904e-06,1.109481e-04}}, //F656N
   {{7.776181e-06,2.758164e-05},{1.270904e-06,1.109481e-04}}, //F658N
   {{7.776181e-06,2.758164e-05},{1.270904e-06,1.109481e-04}}, //F673N
   {{7.776181e-06,2.758164e-05},{1.270904e-06,1.109481e-04}}, //F675W
   {{7.776181e-06,2.758164e-05},{1.270904e-06,7.000852e-05}}, //F702W
   {{7.917809e-06,2.726974e-05},{1.781239e-05,6.953957e-05}}, //F785LP
   {{7.917809e-06,2.726974e-05},{1.781239e-05,9.479239e-05}}, //F791W
   {{7.917809e-06,2.726974e-05},{1.781239e-05,9.479239e-05}}, //F814W
   {{7.917809e-06,2.726974e-05},{1.781239e-05,9.479239e-05}}, //F850LP
   {{7.917809e-06,2.726974e-05},{1.781239e-05,9.479239e-05}}, //F953N
   {{7.917809e-06,2.726974e-05},{1.781239e-05,9.479239e-05}}  //F1042M
};

//gain ratio corrections;
double gcorr_warm[4]={3.719247e-02,-1.203985e-02,6.861146e-03,3.992702e-03};

//filter zero point corrections;
double fcorr_warm[2][NFILTERS_DEFINED]={
   {0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-2.304644e-02,0.000000e+00,0.000000e+00,-9.948792e-02,0.000000e+00,-7.099982e-02,0.000000e+00,-1.567542e-02,-9.765530e-03,4.599253e-02,0.000000e+00,0.000000e+00,0.000000e+00,-3.731368e-03,-9.215336e-03,7.585554e-03,0.000000e+00,9.387174e-03,-1.170134e-03,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,-6.744189e-03,-1.713857e-02,4.278063e-03,3.849086e-02,1.230014e-02,1.362391e-02,0.000000e+00,1.583694e-02},
   {0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,4.467948e-02,0.000000e+00,0.000000e+00,-5.021531e-02,0.000000e+00,-2.172721e-02,0.000000e+00,3.359719e-02,3.950708e-02,9.526514e-02,0.000000e+00,0.000000e+00,0.000000e+00,1.863842e-02,1.315445e-02,2.995534e-02,0.000000e+00,3.175696e-02,1.809453e-02,0.000000e+00,0.000000e+00,0.000000e+00,0.000000e+00,1.252047e-02,2.126090e-03,-2.218946e-02,1.202334e-02,-1.416738e-02,-1.284361e-02,0.000000e+00,-1.063058e-02}
};
// end of CTE parameters

static void readdata(void) {
   int i,j,n,xct;
   FILE *f;
   char str[1025];
   //double eecorr[2];

   sprintf(str,"%s/wfpc2/data/filters.dat","/Users/patito/Documents/topicos/topicos/ModuloII/dolphot2.0/");
   if ((f=fopen(str,"r"))==NULL) {
      printf("Could not find %s\n",str);
      exit(0);
   }
   fscanf(f,"%d",&WFPC2_NFILTERS);
   if (WFPC2_NFILTERS!=NFILTERS_DEFINED)
   {
      printf("Incorrect number of filters in filters.dat.  Must be %d\n",NFILTERS_DEFINED);
      exit(-1);
   }
   fgets(str,81,f);
   WFPC2filters=(WFPC2filttype*)calloc(sizeof(WFPC2filttype),WFPC2_NFILTERS);
   assert(WFPC2filters!=NULL);
   for (i=0;i<WFPC2_NFILTERS;i++) {
      fgets(WFPC2filters[i].name,11,f); WFPC2filters[i].name[strlen(WFPC2filters[i].name)-1]=0;
      WFPC2filters[i].color=fgetc(f);
      fgets(str,81,f);
      fscanf(f,"%lf",WFPC2filters[i].zp);
      for (j=1;j<4;j++) WFPC2filters[i].zp[j]=WFPC2filters[i].zp[0];
      fgets(str,81,f);
      for (j=0;j<8;j++) WFPC2filters[i].xorder[j]='X';
      fscanf(f,"%d",&n);
      fgets(str,81,f);
      xct=0;
      for (j=0;j<n;j++) {
	 do {WFPC2filters[i].xorder[xct]=fgetc(f);} while (WFPC2filters[i].xorder[xct]==' ');
	 fscanf(f,"%lf %lf %lf %lf %lf",WFPC2filters[i].xform[xct],WFPC2filters[i].xform[xct]+1,WFPC2filters[i].xform[xct]+2,WFPC2filters[i].xform[xct]+3,WFPC2filters[i].xform[xct]+4);
	 WFPC2filters[i].xform[xct][2]-=WFPC2filters[i].zp[0];
	 xct++;
	 fgets(str,81,f);
      }
      //applying encircled energy correction afterwards since .xform is differential;
      //WFPC2filters[i].zp[0]-=eecorr[0];
      //WFPC2filters[i].zp[1]-=eecorr[1];
      //WFPC2filters[i].zp[2]-=eecorr[1];
   }
   fclose(f);
   return;
}

int WFPC2findfilt(char*str) {
   int i;
   void WFPC2initfilters();
   if (WFPC2_NFILTERS<0) WFPC2initfilters();
   for (i=0;i<WFPC2_NFILTERS;i++) if (!strcmp(WFPC2filters[i].name,str)) return i;
   printf("Illegal filter %s\n",str);
   exit(0);
   return -1;
}

/*
void readidcfile(char*base) {
   int i,j,c,d;
   FILE *f;
   char str[1025],tmp[81];

   sprintf(str,"%s/wfpc2/data/%s_idctab.dat",BASEDIR,base);
   if ((f=fopen(str,"r"))==NULL) {
      printf("Could not find %s\n",str);
      exit(0);
   }
   while (fgets(str,1025,f)) {
      c=atoi(str);
      if (c<1 || c>2) {
	 printf("Stupid parse error\n%s",str);
	 exit(-1);
      }
      if (!strcmp(base,"hrc")) c=0;
      if (!strncmp(str+6,"FORWARD",7)) d=0;
      else if (!strncmp(str+6,"INVERSE",7)) d=1;
      else {
	 printf("Stupid parse error\n%s",str);
	 exit(-1);
      }
      if (d==1) {
	 if (strncmp(str+17,"CLEAR",5)) memcpy(tmp,str+17,12);
	 else memcpy(tmp,str+30,12);
	 for (i=12;i>0 && tmp[i-1]==' ';i--);
	 tmp[i]=0;
	 i=findfilt(tmp);
	 for (j=0;j<33;j++) WFPC2filters[i].idc[d][c][j]=atof(str+66+j*14);
      }
      else {
	 if (strncmp(str+15,"CLEAR",5)) memcpy(tmp,str+15,8);
	 else memcpy(tmp,str+24,8);
	 for (i=8;i>0 && tmp[i-1]==' ';i--);
	 tmp[i]=0;
	 i=findfilt(tmp);
	 WFPC2filters[i].idc[d][c][0]=atof(str+50);
	 WFPC2filters[i].idc[d][c][1]=atof(str+61);
	 WFPC2filters[i].idc[d][c][2]=atof(str+72);
	 WFPC2filters[i].idc[d][c][3]=atof(str+85);
	 WFPC2filters[i].idc[d][c][4]=atof(str+98);
	 for (j=5;j<33;j++) WFPC2filters[i].idc[d][c][j]=atof(str+111+(j-5)*21);
      }
   }
   fclose(f);
   return;
}

void readidc(void) {
   readidcfile("hrc");
   readidcfile("wfc");
   return;
}
*/

void WFPC2initfilters(void) {
   readdata();
   //readidc();
   return;
}

static double findcolor(double ci[],char color,int *n) {
   int i;
   if (WFPC2_NFILTERS<0) WFPC2initfilters();
   for (i=0;i<WFPC2_NFILTERS;i++) if (WFPC2filters[i].color==color && ci[i]<99) {
      if (n!=NULL) *n=i;
      return ci[i];
   }
   return 99.999;
}

void WFPC2transform(double vmag[],double dvmag[],double smag[]) {
   int f,i,j,CONT=1,it=0;
   double c,tm,tn;
   static double*smag0=NULL;

   if (WFPC2_NFILTERS<0) WFPC2initfilters();
   if (smag0==NULL) {
      smag0=(double*)calloc(sizeof(double),WFPC2_NFILTERS);
      assert(smag0!=NULL);
   }
   memcpy(smag,vmag,8*WFPC2_NFILTERS);
   while (CONT) {
      CONT=0;
      it++;
      memcpy(smag0,smag,8*WFPC2_NFILTERS);
      for (f=0;f<WFPC2_NFILTERS;f++) smag[f]=99.999;
      for (f=0;f<WFPC2_NFILTERS;f++) if (smag0[f]<99) {
	 tm=tn=0;
	 for (i=0;i<8;i++) if (WFPC2filters[f].xorder[i]!='X' && (c=findcolor(smag0,WFPC2filters[f].xorder[i],&j))<99) {
	    if (j<f) c=c-smag0[f];
	    else c=smag0[f]-c;
	    if ((WFPC2filters[f].xform[i][3]==-99 || c>WFPC2filters[f].xform[i][3]) && (WFPC2filters[f].xform[i][4]==99 || c<=WFPC2filters[f].xform[i][4])) {
	       tm+=(WFPC2filters[f].xform[i][0]*c+WFPC2filters[f].xform[i][1]*c*c+WFPC2filters[f].xform[i][2])/(dvmag[j]+0.01);
	       tn+=1./(dvmag[j]+0.01);
	    }
	 }
	 if (tn>0) {
	    smag[f]=vmag[f]+tm/tn;
	    if (it>5) smag[f]=(smag[f]*5+smag0[f]*(it-5))/it;
	    if (fabs(smag[f]-smag0[f])>0.001) CONT=1;
	 }
      }
   }
}

double WFPC2untransform(int filt,double smag[5]) {
   int i,j,f0=-1;
   double c,tm=0.,tn=0.;
   char UBVRI[6]="UBVRI";

   if (WFPC2_NFILTERS<0) WFPC2initfilters();
   for (i=0;i<5;i++) if (WFPC2filters[filt].color==UBVRI[i]) f0=i;
   if (f0<0 || smag[f0]>90) return 99999.;
   for (i=0;i<8;i++) for (j=1;j<5;j++) if (j!=f0 && WFPC2filters[filt].xorder[i]==UBVRI[j] && smag[j]<99) {
      if (j<f0) c=smag[j]-smag[f0];
      else c=smag[f0]-smag[j];
      if ((WFPC2filters[filt].xform[i][3]==-99 || c>WFPC2filters[filt].xform[i][3]) && (WFPC2filters[filt].xform[i][4]==99 || c<=WFPC2filters[filt].xform[i][4])) {
	 tm-=WFPC2filters[filt].xform[i][0]*c+WFPC2filters[filt].xform[i][1]*c*c+WFPC2filters[filt].xform[i][2];
	 tn++;
      }
   }
   if (!tn) for (i=0;i<8;i++) for (j=0;j<1;j++) if (j!=f0 && WFPC2filters[filt].xorder[i]==UBVRI[j] && smag[j]<99) {
      if (j<f0) c=smag[j]-smag[f0];
      else c=smag[f0]-smag[j];
      if ((WFPC2filters[filt].xform[i][3]==-99 || c>WFPC2filters[filt].xform[i][3]) && (WFPC2filters[filt].xform[i][4]==99 || c<=WFPC2filters[filt].xform[i][4])) {
	 tm-=WFPC2filters[filt].xform[i][0]*c+WFPC2filters[filt].xform[i][1]*c*c+WFPC2filters[filt].xform[i][2];
	 tn++;
      }
   }
   if (tn>0) return tm/tn+smag[f0];
   return 99999.;
}

// correction in magnitudes from Dolphin (2009)
double WFPC2_coldCTE(int chip,float x,float y,float cts,float ctmult,int ign,float sky,float epoch,int temp,int filt) {
   double mult,lct,lbg,yr,dm,cm1,cm2;

   if (cts<=0.) return fcorr[filt][ign][chip];
   //Convert to electrons in original image level;
   mult=1./ctmult;
   if (ign==1) mult*=14;
   else mult*=7;
   //sky-=cts*bgcorr[filt][(chip==0)?0:1];
   sky*=mult;
   cts*=mult;
   //Soften background and signal;
   if (sky<0) sky=1.;
   else sky=sqrt(1.+sky*sky);
   if (cts<0) cts=1.;
   else cts=sqrt(1.+cts*cts);
   lbg=log(sky)-1;
   //sky-=10;
   yr=(epoch-49461.854754)/365.25;
   dm=0.007730*exp(-0.504911*lbg)*(1+0.099582*yr)*x/800.;
   lct=log(cts)-7+0.92103404*dm;
   cm1=1.0-0.200904*lbg+0.039005*lbg*lct+0.002108*lct;
   if (cm1<0.153099) cm1=0.153099;
   cm2=0.957704*(yr-0.025474*yr*yr)*exp(-0.449956*lct);
   if (y<0) y=0.0;
   else if (y>800.0) y=800.0;
   return fcorr[filt][ign][chip]+dm+2.412981*log(exp(0.022389*cm1*y/800.)*(1+cm2)-cm2);
}

// correction in magnitudes from Dolphin (2000)
double WFPC2_CTE(int chip,float x,float y,float cts,float ctmult,int ign,float sky,float epoch,int temp,int filt) {
   double mult,XCTE,YCTE,lct,lbg,yr;
   if (temp==0 && WFPC2filters[filt].color!='X' && WFPC2filters[filt].color!='U') return WFPC2_coldCTE(chip,x,y,cts,ctmult,ign,sky,epoch,temp,filt);
   if (cts<=0.) return fcorr_warm[temp][filt]+gcorr_warm[chip];
   //Convert to electrons in original image level;
   mult=1./ctmult;
   if (ign==1) mult*=14;
   else mult*=7;
   //sky-=cts*bgcorr_warm[filt][temp][(chip==0)?0:1];
   sky*=mult;
   cts*=mult;
   //Soften background and signal;
   if (sky<0) sky=1.;
   else sky=sqrt(1.+sky*sky);
   if (cts<0) cts=1.;
   else cts=sqrt(1.+cts*cts);
   //HSTphot CTE determination;
   lct=log(cts)-7;
   lbg=log(sky)-1;
   if (epoch<49462.5) yr=-2.;
   else yr=(epoch-49462.5)/365.25-2.;
   if (temp==0) YCTE=0.018009+(0.087972+exp(-0.506635*lct))*exp(-0.042368*sky-0.035149*lbg)*(0.096532+0.041442*yr);
   else YCTE=0.102822+0.027988*exp(-0.958904*lct);
   if (temp==0) XCTE=(exp(-0.195966*lct))*exp(-0.126435*lbg)*(0.023881+0.002327*yr);
   else XCTE=0.;
   return fcorr_warm[temp][filt]+gcorr_warm[chip]+YCTE*y/800.+XCTE*x/800.;
}
