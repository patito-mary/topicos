#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

typedef struct {
   char name[11],color;
   double zp[3];
   char xorder[2][8];
   int xformc[2][8];
   double xform[2][8][5];
   double idc[2][3][33];
   double dzp[3];
} ACSfilttype;

int ACS_NFILTERS=-1;
ACSfilttype *ACSfilters;

/*
// New CTE Parameters (cold only);
//background "contamination";
double bgcorr[ACS_NFILTERS][2]={
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
double fcorr[ACS_NFILTERS][2][4]={
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
   {{-0.094952,-0.114463,-0.134266,-0.141137},{-0.091563,-0.124748,-0.126559,-0.129744}}, //F380W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F390N
   {{-0.080175,-0.099687,-0.119490,-0.126360},{-0.076787,-0.109971,-0.111782,-0.114968}}, //F410M
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F437N
   {{-0.016375,-0.035886,-0.055689,-0.062560},{-0.012986,-0.046170,-0.047982,-0.051167}}, //F439W
   {{-0.010060,-0.029571,-0.049374,-0.056245},{-0.006671,-0.039855,-0.041666,-0.044852}}, //F450W
   {{0.042616,0.023105,0.003302,-0.003569},{0.046005,0.012821,0.011009,0.007824}}, //F467M
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F469N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F487N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F502N
   {{0.020200,-0.016446,-0.032966,-0.017032},{-0.000606,-0.017946,-0.032272,-0.029324}}, //F547M
   {{0.012311,-0.024335,-0.040855,-0.024921},{-0.008495,-0.025835,-0.040160,-0.037213}}, //F555W
   {{0.035053,-0.001594,-0.018113,-0.002179},{0.014246,-0.003094,-0.017419,-0.014471}}, //F569W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F588N
   {{0.035776,-0.000870,-0.017390,-0.001456},{0.014970,-0.002370,-0.016696,-0.013748}}, //F606W
   {{-0.022031,-0.036042,-0.033863,-0.030176},{0.002250,-0.042243,-0.045907,-0.024163}}, //F622W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F631N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F656N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F658N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F673N
   {{-0.021277,-0.035288,-0.033109,-0.029422},{0.003003,-0.041489,-0.045153,-0.023409}}, //F675W
   {{-0.022286,-0.036297,-0.034118,-0.030431},{0.001994,-0.042498,-0.046162,-0.024418}}, //F702W
   {{0.021840,-0.013261,-0.022826,-0.010518},{-0.008700,-0.017838,-0.018059,-0.011740}}, //F785LP
   {{0.051464,0.016363,0.006798,0.019106},{0.020924,0.011786,0.011565,0.017884}}, //F791W
   {{0.023574,-0.011527,-0.021092,-0.008784},{-0.006966,-0.016104,-0.016325,-0.010006}}, //F814W
   {{0.021980,-0.013121,-0.022686,-0.010378},{-0.008560,-0.017698,-0.017919,-0.011600}}, //F850LP
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F953N
   {{0.023574,-0.011527,-0.021092,-0.008784},{-0.006966,-0.016104,-0.016325,-0.010006}} //F1042M
};
*/

static void readdata(void) {
   int i,j,n,c,xct[2];
   FILE *f;
   char str[1025];
   double eecorr[2];

   sprintf(str,"%s/acs/data/filters.dat","/Users/patito/Documents/topicos/topicos/ModuloII");
   if ((f=fopen(str,"r"))==NULL) {
      printf("Could not find %s\n",str);
      exit(0);
   }
   fscanf(f,"%d",&ACS_NFILTERS);
   fgets(str,81,f);
   ACSfilters=(ACSfilttype*)calloc(sizeof(ACSfilttype),ACS_NFILTERS);
   assert(ACSfilters!=NULL);
   for (i=0;i<ACS_NFILTERS;i++) {
      fgets(ACSfilters[i].name,11,f); ACSfilters[i].name[strlen(ACSfilters[i].name)-1]=0;
      ACSfilters[i].color=fgetc(f);
      fgets(str,81,f);
      fscanf(f,"%lf %lf",ACSfilters[i].zp,ACSfilters[i].zp+1);
      ACSfilters[i].zp[2]=ACSfilters[i].zp[1];
      fscanf(f,"%lf %lf",eecorr,eecorr+1);
      fscanf(f,"%lf %lf %lf",ACSfilters[i].dzp,ACSfilters[i].dzp+1,ACSfilters[i].dzp+2);
      ACSfilters[i].dzp[0]-=ACSfilters[i].zp[0]; // HRC
      ACSfilters[i].dzp[1]-=ACSfilters[i].zp[1]; // WFC pre 7/4/06
      ACSfilters[i].dzp[2]-=ACSfilters[i].zp[1]; // WFC post 7/4/06
      fgets(str,81,f);
      for (j=0;j<8;j++) ACSfilters[i].xorder[0][j]=ACSfilters[i].xorder[1][j]='X';
      fscanf(f,"%d",&n);
      fgets(str,81,f);
      xct[0]=xct[1]=0;
      for (j=0;j<n;j++) {
	 fscanf(f,"%d",&c);
	 do {ACSfilters[i].xorder[c][xct[c]]=fgetc(f);} while (ACSfilters[i].xorder[c][xct[c]]==' ');
	 fscanf(f,"%lf %lf %lf %lf %lf",ACSfilters[i].xform[c][xct[c]]+2,ACSfilters[i].xform[c][xct[c]],ACSfilters[i].xform[c][xct[c]]+1,ACSfilters[i].xform[c][xct[c]]+3,ACSfilters[i].xform[c][xct[c]]+4);
	 ACSfilters[i].xform[c][xct[c]][2]-=ACSfilters[i].zp[c];
	 xct[c]++;
	 fgets(str,81,f);
      }
      //applying encircled energy correction afterwards since .xform is differential;
      ACSfilters[i].zp[0]-=eecorr[0];
      ACSfilters[i].zp[1]-=eecorr[1];
      ACSfilters[i].zp[2]-=eecorr[1];
   }
   fclose(f);
   return;
}

int ACSfindfilt(char*str) {
   int i;
   void ACSinitfilters();
   if (ACS_NFILTERS<0) ACSinitfilters();
   for (i=0;i<ACS_NFILTERS;i++) if (!strcmp(ACSfilters[i].name,str)) return i;
   printf("Illegal filter %s\n",str);
   exit(0);
   return -1;
}

static void readidcfile(char*base) {
   int i,j,c,d;
   FILE *f;
   char str[1025],tmp[81];

   sprintf(str,"%s/acs/data/%s_idctab.dat","/Users/patito/Documents/topicos/topicos/ModuloII",base);
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
	 i=ACSfindfilt(tmp);
	 for (j=0;j<33;j++) ACSfilters[i].idc[d][c][j]=atof(str+66+j*14);
      }
      else {
	 if (strncmp(str+15,"CLEAR",5)) memcpy(tmp,str+15,8);
	 else memcpy(tmp,str+24,8);
	 for (i=8;i>0 && tmp[i-1]==' ';i--);
	 tmp[i]=0;
	 i=ACSfindfilt(tmp);
	 ACSfilters[i].idc[d][c][0]=atof(str+50) - 0.5;
	 ACSfilters[i].idc[d][c][1]=atof(str+61) - 0.5;
	 ACSfilters[i].idc[d][c][2]=atof(str+72);
	 ACSfilters[i].idc[d][c][3]=atof(str+85);
	 ACSfilters[i].idc[d][c][4]=atof(str+98);
	 for (j=5;j<33;j++) ACSfilters[i].idc[d][c][j]=atof(str+111+(j-5)*21);
      }
   }
   fclose(f);
   return;
}

static void readidc(void) {
   int i,c;
   readidcfile("hrc");
   readidcfile("wfc");
   for (i=0;i<ACS_NFILTERS;i++) {
      ACSfilters[i].idc[0][0][2] = 0;
      ACSfilters[i].idc[0][0][3] = 0;
      double v2_1 = 0.5*(ACSfilters[i].idc[0][1][2]-ACSfilters[i].idc[0][2][2]);
      double v3_1 = 0.5*(ACSfilters[i].idc[0][1][3]-ACSfilters[i].idc[0][2][3]);
      ACSfilters[i].idc[0][1][2] = v2_1; // delta x
      ACSfilters[i].idc[0][1][3] = -v3_1; // delta y
      ACSfilters[i].idc[0][2][2] = -ACSfilters[i].idc[0][1][2];
      ACSfilters[i].idc[0][2][3] = -ACSfilters[i].idc[0][1][3];
      for (c=0;c<3;c++) {
	 ACSfilters[i].idc[1][c][2] = ACSfilters[i].idc[0][c][2];
	 ACSfilters[i].idc[1][c][3] = ACSfilters[i].idc[0][c][3];
      }
   }
   return;
}

void ACSinitfilters(void) {
   readdata();
   readidc();
   return;
}

static double findcolor(double ci[],char color,int *n) {
   int i;
   if (ACS_NFILTERS<0) ACSinitfilters();
   for (i=0;i<ACS_NFILTERS;i++) if (ACSfilters[i].color==color && ci[i]<99) {
      if (n!=NULL) *n=i;
      return ci[i];
   }
   return 99.999;
}

void ACStransform(int cm,double vmag[],double dvmag[],double smag[]) {
   int f,i,j,CONT=1,it=0;
   double c,tm,tn;
   static double*smag0=NULL;

   if (ACS_NFILTERS<0) ACSinitfilters();
   if (smag0==NULL) {
      smag0=(double*)calloc(sizeof(double),ACS_NFILTERS);
      assert(smag0!=NULL);
   }
   memcpy(smag,vmag,8*ACS_NFILTERS);
   while (CONT) {
      CONT=0;
      it++;
      memcpy(smag0,smag,8*ACS_NFILTERS);
      for (f=0;f<ACS_NFILTERS;f++) smag[f]=99.999;
      if (cm<0) return;
      for (f=0;f<ACS_NFILTERS;f++) if (smag0[f]<99) {
	 tm=tn=0;
	 for (i=0;i<8;i++) if (ACSfilters[f].xorder[cm][i]!='X' && (c=findcolor(smag0,ACSfilters[f].xorder[cm][i],&j))<99) {
	    if (j<f) c=c-smag0[f];
	    else c=smag0[f]-c;
	    if ((ACSfilters[f].xform[cm][i][3]==-99 || c>ACSfilters[f].xform[cm][i][3]) && (ACSfilters[f].xform[cm][i][4]==99 || c<=ACSfilters[f].xform[cm][i][4])) {
	       tm+=(ACSfilters[f].xform[cm][i][0]*c+ACSfilters[f].xform[cm][i][1]*c*c+ACSfilters[f].xform[cm][i][2])/(dvmag[j]+0.01);
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

double ACSuntransform(int cm,int filt,double smag[5]) {
   int i,j,f0=-1;
   double c,tm=0.,tn=0.;
   char UBVRI[6]="UBVRI";

   if (ACS_NFILTERS<0) ACSinitfilters();
   for (i=0;i<5;i++) if (ACSfilters[filt].color==UBVRI[i]) f0=i;
   if (f0<0 || smag[f0]>90) return 99999.;
   for (i=0;i<8;i++) for (j=1;j<5;j++) if (j!=f0 && ACSfilters[filt].xorder[cm][i]==UBVRI[j] && smag[j]<99) {
      if (j<f0) c=smag[j]-smag[f0];
      else c=smag[f0]-smag[j];
      if ((ACSfilters[filt].xform[cm][i][3]==-99 || c>ACSfilters[filt].xform[cm][i][3]) && (ACSfilters[filt].xform[cm][i][4]==99 || c<=ACSfilters[filt].xform[cm][i][4])) {
	 tm-=ACSfilters[filt].xform[cm][i][0]*c+ACSfilters[filt].xform[cm][i][1]*c*c+ACSfilters[filt].xform[cm][i][2];
	 tn++;
      }
   }
   if (!tn) for (i=0;i<8;i++) for (j=0;j<1;j++) if (j!=f0 && ACSfilters[filt].xorder[cm][i]==UBVRI[j] && smag[j]<99) {
      if (j<f0) c=smag[j]-smag[f0];
      else c=smag[f0]-smag[j];
      if ((ACSfilters[filt].xform[cm][i][3]==-99 || c>ACSfilters[filt].xform[cm][i][3]) && (ACSfilters[filt].xform[cm][i][4]==99 || c<=ACSfilters[filt].xform[cm][i][4])) {
	 tm-=ACSfilters[filt].xform[cm][i][0]*c+ACSfilters[filt].xform[cm][i][1]*c*c+ACSfilters[filt].xform[cm][i][2];
	 tn++;
      }
   }
   if (tn>0) return tm/tn+smag[f0];
   return 99999.;
}

double ACS_CTE(int chip,float x,float y,float cts,float ctmult,float gn,float sky,float epoch) {
   double YCTE,SKY,FLUX;

   SKY=sky*gn/ctmult;
   if (SKY<0.2) SKY=0.2;
   FLUX=cts*gn/ctmult;
   if (FLUX<0.) FLUX=50.;
   else FLUX=sqrt(FLUX*FLUX+2500.);
   if (chip==0) {
      // new correction from Chiaberge et al.
      YCTE=0.363*pow(SKY,-0.15)*pow(FLUX,-0.36)*(epoch-52333.)/365.;
      if (y<0) YCTE = 0.0;
      else if (y>1024) YCTE *= 1.024;
      else YCTE *= y/1000.;
      // old correction from ISR03-09
      //YCTE=0.87*pow(SKY,-0.27)*pow(FLUX,-0.21)*y/2048.*(epoch-52333.)/381.;
   }
   else {
      // new correction from Chiaberge et al.
      YCTE=0.708*pow(SKY,-0.25)*pow(FLUX,-0.44)*(epoch-52333.)/365.;
      if (chip==1) y = 2048.0-y;
      if (y<0) YCTE = 0.0;
      else if (y>2048) YCTE *= 1.024;
      else YCTE *= y/2000.0;
      // old correction from ISR04-06
      //YCTE=0.0037818747*pow(SKY,-0.31)*pow(FLUX,-0.64)*(epoch-52333.);
      //if (chip==1) YCTE*=1.-y/2048.;
      //else YCTE*=y/2048.;
   }
   //printf("chip=%d, mlt=%3.1f, SKY=%4.2f, FLUX=%4.2f, EPOCH=%3.1f, Y=%3.1f, YCTE=%f\n",chip,ctmult,SKY,FLUX,epoch,y,YCTE);
   return YCTE;
}

double ACS_ZP(int filt,int cm,double MJD) {
   if (ACS_NFILTERS<0) ACSinitfilters();
   if (cm==0) return ACSfilters[filt].zp[cm]+ACSfilters[filt].dzp[0];
   if (MJD<53920) return ACSfilters[filt].zp[cm]+ACSfilters[filt].dzp[1];
   return ACSfilters[filt].zp[cm]+ACSfilters[filt].dzp[2];
}
