#include "../include/dolphot.h"
#include "wfpc2psfdata.h"

int main(int argc,char**argv) {
   int cm,x,y,fx,fy,j,k,l,m,n;
   FILE *fpsf;
   char str[161];
   float **psf;

   if (argc!=7) {
      printf("Usage: %s <filt> <chip> <Xpos> <Ypos> <dX> <dY>\n",*argv);
      return -1;
   }
   cm=atoi(argv[2])-1;
   if (cm<0 || cm>3) {
      printf("Chip must be 1-4\n");
      return -1;
   }
   x=atoi(argv[3]);
   if (x<0 || x>=wfpc2_nxpsfpos[cm]) {
      printf("X position must be 0-%d\n",wfpc2_nxpsfpos[cm]-1);
      return -1;
   }
   y=atoi(argv[4]);
   if (y<0 || y>=wfpc2_nypsfpos[cm]) {
      printf("Y position must be 0-%d\n",wfpc2_nypsfpos[cm]-1);
      return -1;
   }
   fx=atoi(argv[5]);
   fy=atoi(argv[6]);
   if (fx<-wfpc2_n2psf[cm] || fx>wfpc2_n2psf[cm] || fy<-wfpc2_n2psf[cm] || fy>wfpc2_n2psf[cm]) {
      printf("dX and dY must be -%d to %d\n",wfpc2_n2psf[cm],wfpc2_n2psf[cm]);
      return -1;
   }
   sprintf(str,"%s/wfpc2/data/%s.%s.psf","/Users/patito/Documents/topicos/topicos/ModuloII/dolphot2.0/",argv[1],wfpc2_cn[cm]);
   if ((fpsf=fopen(str,"rb"))==NULL) {
      printf("Cannot open %s\n",str);
      return -1;
   }
   psf=(float**)calloc(sizeof(float*),wfpc2_rpsf[cm]*2+1);
   if (!psf) merr();
   psf+=wfpc2_rpsf[cm];
   for (j=-wfpc2_rpsf[cm];j<=wfpc2_rpsf[cm];j++) {
      psf[j]=(float*)calloc(sizeof(float),wfpc2_rpsf[cm]*2+1);
      if (!psf[j]) merr();
      psf[j]+=wfpc2_rpsf[cm];
   }
   for (j=0;j<y;j++) for (k=0;k<wfpc2_nxpsfpos[cm];k++) for (l=-wfpc2_n2psf[cm];l<=wfpc2_n2psf[cm];l++) for (m=-wfpc2_n2psf[cm];m<=wfpc2_n2psf[cm];m++) for (n=-wfpc2_rpsf[cm];n<=wfpc2_rpsf[cm];n++) fread(psf[n]-wfpc2_rpsf[cm],4,2*wfpc2_rpsf[cm]+1,fpsf);
   for (k=0;k<x;k++) for (l=-wfpc2_n2psf[cm];l<=wfpc2_n2psf[cm];l++) for (m=-wfpc2_n2psf[cm];m<=wfpc2_n2psf[cm];m++) for (n=-wfpc2_rpsf[cm];n<=wfpc2_rpsf[cm];n++) fread(psf[n]-wfpc2_rpsf[cm],4,2*wfpc2_rpsf[cm]+1,fpsf);
   for (l=-wfpc2_n2psf[cm];l<fy;l++) for (m=-wfpc2_n2psf[cm];m<=wfpc2_n2psf[cm];m++) for (n=-wfpc2_rpsf[cm];n<=wfpc2_rpsf[cm];n++) fread(psf[n]-wfpc2_rpsf[cm],4,2*wfpc2_rpsf[cm]+1,fpsf);
   for (m=-wfpc2_n2psf[cm];m<fx;m++) for (n=-wfpc2_rpsf[cm];n<=wfpc2_rpsf[cm];n++) fread(psf[n]-wfpc2_rpsf[cm],4,2*wfpc2_rpsf[cm]+1,fpsf);
   for (n=-wfpc2_rpsf[cm];n<=wfpc2_rpsf[cm];n++) ffread(psf[n]-wfpc2_rpsf[cm],4,2*wfpc2_rpsf[cm]+1,fpsf);
   k=wfpc2_rpsf[cm];
   if (k>6) k=6;
   for (n=k;n>=-k;n--) {
      for (j=-k;j<=k;j++) printf("%5d ",(int)(psf[n][j]*100000+0.5));
      printf("\n");
   }
   fclose(fpsf);
   return 0;
}
