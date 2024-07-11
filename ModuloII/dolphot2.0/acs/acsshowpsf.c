#include "../include/dolphot.h"
#include "acspsfdata.h"

void usage(char *exe) {
   printf("Usage: %s <filt> <chip> <Xpos> <Ypos> <dX> <dY> <<flags>>\n",exe);
   printf("  -anderson  to use PSFs from Jay Anderson's library\n");
   exit(-1);
}

int main(int argc,char**argv) {
   int cm,x,y,fx,fy,j,k,l,m,n,psfType=0;
   FILE *fpsf;
   char str[161];
   float **psf;

   if (argc<7) usage(*argv);
   for (j=7;j<argc;j++) {
      if (!strcasecmp(argv[j],"-anderson")) psfType=1;
      else usage(*argv);
   }
   cm=atoi(argv[2]);
   if (cm<0 || cm>2) {
      printf("Chip must be 0-2\n");
      return -1;
   }
   x=atoi(argv[3]);
   if (x<0 || x>=acs_nxpsfpos[cm]) {
      printf("X position must be 0-%d\n",acs_nxpsfpos[cm]-1);
      return -1;
   }
   y=atoi(argv[4]);
   if (y<0 || y>=acs_nypsfpos[cm]) {
      printf("Y position must be 0-%d\n",acs_nypsfpos[cm]-1);
      return -1;
   }
   fx=atoi(argv[5]);
   fy=atoi(argv[6]);
   if (fx<-acs_n2psf[cm] || fx>acs_n2psf[cm] || fy<-acs_n2psf[cm] || fy>acs_n2psf[cm]) {
      printf("dX and dY must be -%d to %d\n",acs_n2psf[cm],acs_n2psf[cm]);
      return -1;
   }
   if (psfType==0) sprintf(str,"%s/acs/data/%s.%s.psf","/Users/patito/Documents/topicos/topicos/ModuloII",argv[1],acs_cn[cm]);
   else sprintf(str,"%s/acs/data/%s_anderson.%s.psf","/Users/patito/Documents/topicos/topicos/ModuloII",argv[1],acs_cn[cm]);
   if ((fpsf=fopen(str,"rb"))==NULL) {
      printf("Cannot open %s\n",str);
      return -1;
   }
   psf=(float**)calloc(sizeof(float*),acs_rpsf[cm]*2+1);
   if (!psf) merr();
   psf+=acs_rpsf[cm];
   for (j=-acs_rpsf[cm];j<=acs_rpsf[cm];j++) {
      psf[j]=(float*)calloc(sizeof(float),acs_rpsf[cm]*2+1);
      if (!psf[j]) merr();
      psf[j]+=acs_rpsf[cm];
   }
   for (j=0;j<y;j++) for (k=0;k<acs_nxpsfpos[cm];k++) for (l=-acs_n2psf[cm];l<=acs_n2psf[cm];l++) for (m=-acs_n2psf[cm];m<=acs_n2psf[cm];m++) for (n=-acs_rpsf[cm];n<=acs_rpsf[cm];n++) fread(psf[n]-acs_rpsf[cm],4,2*acs_rpsf[cm]+1,fpsf);
   for (k=0;k<x;k++) for (l=-acs_n2psf[cm];l<=acs_n2psf[cm];l++) for (m=-acs_n2psf[cm];m<=acs_n2psf[cm];m++) for (n=-acs_rpsf[cm];n<=acs_rpsf[cm];n++) fread(psf[n]-acs_rpsf[cm],4,2*acs_rpsf[cm]+1,fpsf);
   for (l=-acs_n2psf[cm];l<fy;l++) for (m=-acs_n2psf[cm];m<=acs_n2psf[cm];m++) for (n=-acs_rpsf[cm];n<=acs_rpsf[cm];n++) fread(psf[n]-acs_rpsf[cm],4,2*acs_rpsf[cm]+1,fpsf);
   for (m=-acs_n2psf[cm];m<fx;m++) for (n=-acs_rpsf[cm];n<=acs_rpsf[cm];n++) fread(psf[n]-acs_rpsf[cm],4,2*acs_rpsf[cm]+1,fpsf);
   for (n=-acs_rpsf[cm];n<=acs_rpsf[cm];n++) ffread(psf[n]-acs_rpsf[cm],4,2*acs_rpsf[cm]+1,fpsf);
   k=acs_rpsf[cm];
   if (k>6) k=6;
   for (n=k;n>=-k;n--) {
      for (j=-k;j<=k;j++) printf("%5d ",(int)(psf[n][j]*100000+0.5));
      printf("\n");
   }
   fclose(fpsf);
   return 0;
}
