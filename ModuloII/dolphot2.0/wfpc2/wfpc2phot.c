#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/fits.h"
#include "../dolphot_defs.h"
#include "wfpc2psfdata.h"
#include "wfpc2filters.h"
#include "wfpc2distort.h"

extern void shift(int img,double x0,double y0,double *x,double *y,int dir);

static float********wfpc2psflib=NULL;
static int WFPC2_RPSF[4]={29,12,12,12};
static int ANY_WFPC2=0;

void wfpc2initparam(void) {
   char str[161],*ptr;
   int img;
   ANY_WFPC2=0;
   for (img=0;img<Timg && !ANY_WFPC2;img++) {
      strcpy(str,getcardval(dataim+img,"DOLWFPC2",0));
      strtol(str,&ptr,10);
      if (ptr!=str) ANY_WFPC2=1;
   }
   if (ANY_WFPC2==0) return;
   FPSF=4;
   SubPixel=1;
   EPSF=0;
   PSFsol=-1;
   PSFStep=0.;
   Zero=0.;
   MinS=0.2;
   MaxS=15.;
   MaxE=0.1;
   for (img=0;img<Timg;img++) {
      if (RPSF[img]>=wfpc2_rpsf[0]) {RPSF[img]=wfpc2_rpsf[0]-1; printf("Lowering RPSF to %d\n",RPSF[img]);}
      if (RAper[img]>wfpc2_rpsf[0]-1.) {RAper[img]=wfpc2_rpsf[0]-1; printf("Lowering RAper to %d\n",wfpc2_rpsf[0]-1);}
   }
   return;
}

void wfpc2radii(void) {
   char str[161],*ptr;
   int img;
   for (img=0;img<Timg;img++) {
      strcpy(str,getcardval(dataim+img,"DOLWFPC2",0));
      int cm=strtol(str,&ptr,10);
      if (ptr==str) {
	 // AEDDEBUG need warning: printf("**Image %d has not been preprocessed with wfpc2mask; cannot proceed\n",img+1);
      }
      else {
	 if (cm<-1 || cm>3) {
	    printf("**Image %d's chip cannot be identified; please report bug.\n",img+1);
	    exit(-1);
	 }
	 hstmode[img].inst=WFPC2;
	 hstmode[img].cm=cm;
	 if (img==Nimg) {
	    if (hstmode[img].cm==-1) {
	       // use WFC3 PSFs for drizzled frame
	       hstmode[img].cm=2;
	       DRIZZLE_BASE=1;
	    }
	    else DRIZZLE_BASE=0;
	 }
	 if (hstmode[img].cm<0) {
	    printf("**Image %d has been drizzled; cannot run photometry\n",img+1);
	    exit(-1);
	 }
	 if (hstmode[img].cm>0) {
	    RAper[img]/=1.5;
	    RSky0[img]/=1.5;
	    RSky1[img]/=1.5;
	    RSky20[img]/=1.5;
	    RSky21[img]/=1.5;
	    RChi[img]/=1.5;
	    RPSF[img]=(RPSF[img]*2+2)/3;
	    apsize[img]/=1.5;
	    apsky[img][0]/=1.5;
	    apsky[img][1]/=1.5;
	 }
      }
   }
}

void wfpc2tweakshift(void) {
   int img;
   for (img=0;img<Timg;img++) {
      if (hstmode[img].inst==WFPC2 && hstmode[img].cm>0) {
	 switch(hstmode[img].cm) {
	 case 1:
	    dpos[img][0]= 0.457*dpos0[img][1];
	    dpos[img][1]=-0.457*dpos0[img][0];
	    break;
	 case 2:
	    dpos[img][0]=-0.457*dpos0[img][0];
	    dpos[img][1]=-0.457*dpos0[img][1];
	    break;
	 case 3:
	    dpos[img][0]=-0.457*dpos0[img][1];
	    dpos[img][1]= 0.457*dpos0[img][0];
	    break;
	 }
      }
   }
}

void wfpc2initpsf(void) {
   int img,i,j,y1,x1,y2,x2,y3,n2,n3,n3skip;
   char str[161];
   FILE *f;
   float***** ptr1;
   float**** ptr2;
   float*** ptr3;
   float** ptr4;
   float* ptr5;
   float Temp;

   if (WFPC2_NFILTERS<0) WFPC2initfilters();
   if (wfpc2psflib==NULL) {
      wfpc2psflib=(float********)calloc(sizeof(float*******),4);
      if (!wfpc2psflib) merr();
      for (i=0;i<4;i++) {
	 wfpc2psflib[i]=(float*******)calloc(sizeof(float******),WFPC2_NFILTERS);
	 if (!wfpc2psflib[i]) merr();
      }
   }
   for (i=0;i<4;i++) for (j=0;j<WFPC2_NFILTERS;j++) wfpc2psflib[i][j]=NULL;
   for (img=0;img<Timg;img++) if (hstmode[img].inst==WFPC2) {
      strcpy(str,getcardval(datahd+img,"FILTNAM1",0));
      if (!strcmp(str,"")) strcpy(str,getcardval(datahd+img,"FILTNAM2",0));
      hstmode[img].filt=WFPC2findfilt(str);
      //printf("Image %d: cm=%d, filt=%d (%s)\n",img+1,hstmode[img].cm,hstmode[img].filt,WFPC2filters[hstmode[img].filt].name);
      if (RPSF[img]>=wfpc2_rpsf[hstmode[img].cm]) {printf("ERROR: RPSF must be less than %d for WFPC2/WFC data\n",wfpc2_rpsf[hstmode[img].cm]); exit(-1);}
      if (RAper[img]>wfpc2_rpsf[hstmode[img].cm]-1.) {printf("ERROR: RAper must be no more than %d for WFPC2/WFC data\n",wfpc2_rpsf[hstmode[img].cm]-1); exit(-1);}
      // Extract temperature
      Temp=atof(getcardval(datahd+img,"UCH1CJTM",1));
      if (Temp<-82) hstmode[img].i1=0;
      else hstmode[img].i1=1;
      // Extract gain
      if (fabs(iGAIN[img]-7.)<0.001) {
	 hstmode[img].i2=0;
	 iGAIN[img]=7.;
      }
      else if (fabs(iGAIN[img]-14.)<0.001 || fabs(iGAIN[img]-15.)<0.001) {
	 hstmode[img].i2=1;
	 iGAIN[img]=14.;
      }
      else {
	 printf("ERROR: WFPC2 GAIN must be 7 or 14.\n");
      }
   }
   for (i=0;i<4;i++) WFPC2_RPSF[i]=0;
   for (img=0;img<Timg;img++) if (hstmode[img].inst==WFPC2) {
      i=hstmode[img].cm;
      apsf[img][0][0]=1.;
      apsf[img][1][0]=1.;
      apsf[img][2][0]=0.;
      if (i==0) apsize[img]=11.;
      else apsize[img]=5.;
      if (WFPC2_RPSF[i]<RPSF[img]) WFPC2_RPSF[i]=RPSF[img];
      if (WFPC2_RPSF[i]<rphot[img]) WFPC2_RPSF[i]=rphot[img];
      for (i=1;i<5;i++) apsf[img][0][i]=apsf[img][1][i]=apsf[img][2][i]=0.;
   }
   for (img=0;img<Timg;img++) if (hstmode[img].inst==WFPC2 && wfpc2psflib[i=hstmode[img].cm][j=hstmode[img].filt]==NULL) {
      sprintf(str,"%s/wfpc2/data/%s.%s.psf","/Users/patito/Documents/topicos/topicos/ModuloII/dolphot2.0/",WFPC2filters[hstmode[img].filt].name,wfpc2_cn[i]);
      if ((f=fopen(str,"rb"))==NULL) {
	 printf("Cannot open %s\n",str);
	 exit(-1);
      }
      n2=2*wfpc2_n2psf[i]+1;
      n3=2*WFPC2_RPSF[i]+1;
      n3skip=wfpc2_rpsf[i]-WFPC2_RPSF[i];
      wfpc2psflib[i][j]=(float******)calloc(sizeof(float*****),wfpc2_nypsfpos[i]);
      ptr1=(float*****)calloc(sizeof(float****),wfpc2_nypsfpos[i]*wfpc2_nxpsfpos[i]);
      ptr2=(float****)calloc(sizeof(float***),wfpc2_nypsfpos[i]*wfpc2_nxpsfpos[i]*n2);
      ptr3=(float***)calloc(sizeof(float**),wfpc2_nypsfpos[i]*wfpc2_nxpsfpos[i]*n2*n2);
      ptr4=(float**)calloc(sizeof(float*),wfpc2_nypsfpos[i]*wfpc2_nxpsfpos[i]*n2*n2*n3);
      ptr5=(float*)calloc(sizeof(float),wfpc2_nypsfpos[i]*wfpc2_nxpsfpos[i]*n2*n2*n3*n3);
      if (!wfpc2psflib[i][j] || !ptr1) merr();
      for (y1=0;y1<wfpc2_nypsfpos[i];y1++) {
	 wfpc2psflib[i][j][y1]=ptr1;
	 ptr1+=wfpc2_nxpsfpos[i];
	 for (x1=0;x1<wfpc2_nxpsfpos[i];x1++) {
	    wfpc2psflib[i][j][y1][x1]=ptr2+wfpc2_n2psf[i];
	    ptr2+=n2;
	    for (y2=-wfpc2_n2psf[i];y2<=wfpc2_n2psf[i];y2++) {
	       wfpc2psflib[i][j][y1][x1][y2]=ptr3+wfpc2_n2psf[i];
	       ptr3+=n2;
	       for (x2=-wfpc2_n2psf[i];x2<=wfpc2_n2psf[i];x2++) {
		  wfpc2psflib[i][j][y1][x1][y2][x2]=ptr4+WFPC2_RPSF[i];
		  ptr4+=n3;
		  if (n3skip) fseek(f,4*n3skip*(2*wfpc2_rpsf[i]+1),SEEK_CUR);
		  for (y3=-WFPC2_RPSF[i];y3<=WFPC2_RPSF[i];y3++) {
		     wfpc2psflib[i][j][y1][x1][y2][x2][y3]=ptr5+WFPC2_RPSF[i];
		     ptr5+=n3;
		     if (n3skip) fseek(f,4*n3skip,SEEK_CUR);
		     ffread(wfpc2psflib[i][j][y1][x1][y2][x2][y3]-WFPC2_RPSF[i],4,n3,f);
		     if (n3skip) fseek(f,4*n3skip,SEEK_CUR);
		  }
		  if (n3skip) fseek(f,4*n3skip*(2*wfpc2_rpsf[i]+1),SEEK_CUR);
	       }
	    }
	 }
      }
      fclose(f);
   }
   return;
}

void wfpc2freepsf(void) {
   int i,j;
   for (i=0;i<4;i++) for (j=0;j<WFPC2_NFILTERS;j++) if (wfpc2psflib[i][j]) {
      free(wfpc2psflib[i][j][0][0][-wfpc2_n2psf[i]][-wfpc2_n2psf[i]][-WFPC2_RPSF[i]]-WFPC2_RPSF[i]);
      free(wfpc2psflib[i][j][0][0][-wfpc2_n2psf[i]][-wfpc2_n2psf[i]]-WFPC2_RPSF[i]);
      free(wfpc2psflib[i][j][0][0][-wfpc2_n2psf[i]]-wfpc2_n2psf[i]);
      free(wfpc2psflib[i][j][0][0]-wfpc2_n2psf[i]);
      free(wfpc2psflib[i][j][0]);
      free(wfpc2psflib[i][j]);
      wfpc2psflib[i][j]=NULL;
   }
   return;
}

//#define DEBUG_PSF
int calcwfpc2psf(int img,float x,float y,int r,int force) {
   int i,j,y1,x1,y2,x2,yy,xx;
   float mx1=1,my1=1,imx1=0,imy1=0;
   float mx2,my2,imx2,imy2;
   static int first=1,lastr=0,lastimg=0;
   static float lastx=0,lasty=0;

   if (hstmode[img].inst!=WFPC2) {
      printf("Stupid error; called wfpc2psf for non-WFPC2 data\n");
      exit(-1);
   }
   if (!first && lastpsftype==1 && img==lastimg && x==lastx && y==lasty && r<=lastr && !poffreset && !force) return 0;
   first=0;lastpsftype=1;lastimg=img;lastx=x;lasty=y;lastr=r;
   i=hstmode[img].cm;
   j=hstmode[img].filt;
#ifdef DEBUG_PSF
   printf("%d %d %d %f %f %d",img,i,j,x,y,r);
   printf("%s/%s PSF at %f,%f; r=%d:\n",wfpc2_cn[i],WFPC2filters[j].name,x,y,r);
   fflush(stdout);
#endif
   imy2=y-(int)y;
   if (imy2<0) imy2++;
   imy2=(imy2-0.5)*wfpc2_sub[i];
   y2=(int)(imy2+wfpc2_sub[i])-wfpc2_sub[i];
   imy2-=y2; my2=1-imy2;
   imx2=x-(int)x;
   if (imx2<0) imx2++;
   imx2=(imx2-0.5)*wfpc2_sub[i];
   x2=(int)(imx2+wfpc2_sub[i])-wfpc2_sub[i];
   imx2-=x2; mx2=1-imx2;
   if (InterpPSFlib && (img<Nimg || !DRIZZLE_BASE)) {
      imx1 = x/100.0 - 0.5;
      imy1 = y/100.0 - 0.5;
      if (imx1<=0.0) {x1=0; imx1=0.0;}
      else if (imx1>=wfpc2_nxpsfpos[i]-1.0) {x1=wfpc2_nxpsfpos[i]-2; imx1=1.0;}
      else {x1=(int)imx1; imx1-=x1;}
      mx1 = 1.0-imx1;
      if (imy1<=0.0) {y1=0; imy1=0.0;}
      else if (imy1>=wfpc2_nypsfpos[i]-1.0) {y1=wfpc2_nypsfpos[i]-2; imy1=1.0;}
      else {y1=(int)imy1; imy1-=y1;}
      my1 = 1.0-imy1;
      for (yy=-r;yy<=r;yy++) for (xx=-r;xx<=r;xx++) {
	 psf[yy][xx] = 
	    ( wfpc2psflib[i][j][y1][x1][y2][x2][yy][xx]*mx2*my2+wfpc2psflib[i][j][y1][x1][y2][x2+1][yy][xx]*imx2*my2+wfpc2psflib[i][j][y1][x1][y2+1][x2][yy][xx]*mx2*imy2+wfpc2psflib[i][j][y1][x1][y2+1][x2+1][yy][xx]*imx2*imy2 ) * mx1*my1 +
	    ( wfpc2psflib[i][j][y1+1][x1][y2][x2][yy][xx]*mx2*my2+wfpc2psflib[i][j][y1+1][x1][y2][x2+1][yy][xx]*imx2*my2+wfpc2psflib[i][j][y1+1][x1][y2+1][x2][yy][xx]*mx2*imy2+wfpc2psflib[i][j][y1+1][x1][y2+1][x2+1][yy][xx]*imx2*imy2 ) * mx1*imy1 +
	    ( wfpc2psflib[i][j][y1][x1+1][y2][x2][yy][xx]*mx2*my2+wfpc2psflib[i][j][y1][x1+1][y2][x2+1][yy][xx]*imx2*my2+wfpc2psflib[i][j][y1][x1+1][y2+1][x2][yy][xx]*mx2*imy2+wfpc2psflib[i][j][y1][x1+1][y2+1][x2+1][yy][xx]*imx2*imy2 ) * imx1*my1 +
	    ( wfpc2psflib[i][j][y1+1][x1+1][y2][x2][yy][xx]*mx2*my2+wfpc2psflib[i][j][y1+1][x1+1][y2][x2+1][yy][xx]*imx2*my2+wfpc2psflib[i][j][y1+1][x1+1][y2+1][x2][yy][xx]*mx2*imy2+wfpc2psflib[i][j][y1+1][x1+1][y2+1][x2+1][yy][xx]*imx2*imy2 ) * imx1*imy1;
      }
   }
   else {
      if (img==Nimg && DRIZZLE_BASE) {
	 y1=4;
	 x1=7;
      }
      else {
	 y1=(int)y/100; if (y1<0) y1=0; if (y1>=wfpc2_nypsfpos[i]) y1=wfpc2_nypsfpos[i]-1;
	 x1=(int)x/100; if (x1<0) x1=0; if (x1>=wfpc2_nxpsfpos[i]) x1=wfpc2_nxpsfpos[i]-1;
      }
      for (yy=-r;yy<=r;yy++) for (xx=-r;xx<=r;xx++) psf[yy][xx]=wfpc2psflib[i][j][y1][x1][y2][x2][yy][xx]*mx2*my2+wfpc2psflib[i][j][y1][x1][y2][x2+1][yy][xx]*imx2*my2+wfpc2psflib[i][j][y1][x1][y2+1][x2][yy][xx]*mx2*imy2+wfpc2psflib[i][j][y1][x1][y2+1][x2+1][yy][xx]*imx2*imy2;
   }
#ifdef DEBUG_PSF
   printf("y1=%d, my1=%f; x1=%d, mx1=%f\n",y1,my1,x1,mx1);
   printf("y2=%d, my2=%f; x2=%d, mx2=%f\n",y2,my2,x2,mx2);
   for (yy=6;yy>=-6;yy--) if (yy>=-r && yy<=r) {for (xx=-6;xx<=6;xx++) if (xx>=-r && xx<=r) printf("%5d ",(int)(psf[yy][xx]*100000+0.5)); printf("\n");} fflush(stdout);
#endif
   return 1;
}

#ifdef DEBUG_PSF
void WFPC2dumpPSFs(void) {
   int img,xphase,yphase,x,y;
   for (img=0;img<Nimg;img++) if (hstmode[img].inst==WFPC2) {
      for (xphase=-5;xphase<=5;xphase++) calcwfpc2psf(img,450.5+0.1*xphase,450.5,RPSF[img],1);
      for (yphase=-5;yphase<=5;yphase++) calcwfpc2psf(img,450.5,450.5+0.1*yphase,RPSF[img],1);
      for (x=0;x<=800;x+=25) calcwfpc2psf(img,0.5+x,450.5,RPSF[img],1);
      for (y=0;y<=800;y+=25) calcwfpc2psf(img,450.5,0.5+y,RPSF[img],1);
   }
}
#endif

double WFPC2calcmag(int img,float x0,float y0,float ct0,float bg,int useCTE) {
   float cm=1.;
   double x,y,m;

   if (hstmode[img].inst!=WFPC2) {
      printf("Stupid error; called wfpc2calcmag for non-WFPC2 data\n");
      exit(-1);
   }
   if (WFPC2_NFILTERS<0) WFPC2initfilters();
   m=-2.5*log10(ct0*apcor[img]/iEXP[img]*wfpc2_ctmult[hstmode[img].i2][hstmode[img].cm])+WFPC2filters[hstmode[img].filt].zp[hstmode[img].cm];
   if (useCTE) {
      shift(img,x0,y0,&x,&y,1);
      if (iEXP0[img]>0.) cm=iEXP[img]/iEXP0[img];
      m-=WFPC2_CTE(hstmode[img].cm,x,y,ct0,cm,hstmode[img].i2,bg,iEPOCH[img],hstmode[img].i1,hstmode[img].filt);
   }
   return m;
}

void WFPC2outstarinfo(FILE *f,int*ct) {
   int i,j,n;
   for (i=0;i<WFPC2_NFILTERS;i++) {
      n=0;
      for (j=0;j<Nimg;j++) if (hstmode[j].inst==WFPC2 && hstmode[j].filt==i) n++;
      if (n>1) {
	 fprintf(f,"%d. Total counts, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Total sky level, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Normalized count rate, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Normalized count rate uncertainty, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Instrumental VEGAMAG magnitude, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Transformed UBVRI magnitude, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Magnitude uncertainty, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Chi, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Signal-to-noise, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Sharpness, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Roundness, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Crowding, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
	 fprintf(f,"%d. Photometry quality flag, WFPC2_%s\n",++(*ct),WFPC2filters[i].name);
      }
   }
}

static double*smag=NULL,*vmag=NULL,*dvmag=NULL;
void WFPC2outstar(FILE *of,float x0,float y0,photdatatype*pdata) {
   int i,j;
   float m=1.;
   static int *fused=NULL;
   static photdatatype*fphot=NULL;
   double x,y;

   if (WFPC2_NFILTERS<0) WFPC2initfilters();
   if (!fphot) {
      fused=(int*)calloc(sizeof(int),WFPC2_NFILTERS);
      fphot=(photdatatype*)calloc(sizeof(photdatatype),WFPC2_NFILTERS);
      smag=(double*)calloc(sizeof(double),WFPC2_NFILTERS);
      vmag=(double*)calloc(sizeof(double),WFPC2_NFILTERS);
      dvmag=(double*)calloc(sizeof(double),WFPC2_NFILTERS);
      if (!fused || !fphot || !smag || !vmag || !dvmag) merr();
   }
   for (j=0;j<Nimg;j++) if (hstmode[j].inst==WFPC2) {
      double dm = -2.5*log10(wfpc2_ctmult[hstmode[j].i2][hstmode[j].cm]) + WFPC2filters[hstmode[j].filt].zp[hstmode[j].cm]-Zero;
      if (pdata[j].ct>0) {
	 if (WFPC2useCTE) {
	    shift(j,x0,y0,&x,&y,1);
	    if (iEXP0[j]>0.) m=iEXP[j]/iEXP0[j];
	    else m=1.0;
	    dm -= WFPC2_CTE(hstmode[j].cm,x,y,pdata[j].ct,m,hstmode[j].i2,pdata[j].sky,iEPOCH[j],hstmode[j].i1,hstmode[j].filt);
	 }
	 pdata[j].m += dm;
      }
      pdata[j].ctcorr *= pow(10,-0.4*dm);
      pdata[j].dctcorr *= pow(10,-0.4*dm);
   }
   for (i=0;i<WFPC2_NFILTERS;i++) {
      double wt,cwt,swt,twt=0.,tcwt=0.,tswt=0.,is,iss,tcm=0.;
      fphot[i].ct0=fphot[i].ct=fphot[i].chi=fphot[i].sh=fphot[i].sky=fphot[i].ctcorr=fphot[i].dctcorr=fphot[i].rnd=fphot[i].crowd=0.;
      fused[i]=0;
      fphot[i].flag=0;
      for (j=0;j<Nimg;j++) if (hstmode[j].inst==WFPC2 && hstmode[j].filt==i) {
	 fused[i]++;
	 if (pdata[j].flag<8 && !(pdata[j].flag&FlagMask)) {
	    is=pdata[j].ct/iEXP[j];
	    iss=pdata[j].dct/iEXP[j];
	    wt=1./iss/iss;
	    cwt=1./pdata[j].dctcorr/pdata[j].dctcorr; // AEDDEBUG new
	    if (CombineChi && pdata[j].chi>1) {
	       wt/=pdata[j].chi*pdata[j].chi;
	       cwt/=pdata[j].chi*pdata[j].chi; // AEDDEBUG new
	    }
	    if (is>0) swt=wt*is;
	    else swt=0.0;
	    twt+=wt;
	    tcwt+=cwt; // AEDDEBUG new
	    tswt+=swt;
	    fphot[i].ct0+=pdata[j].ct0/iEXP[j]*wt;
	    fphot[i].ct+=is*wt;
	    fphot[i].chi+=pdata[j].chi*pdata[j].chi*cwt; // AEDDEBUG was wt
	    fphot[i].sh+=pdata[j].sh*swt;
	    fphot[i].sky+=pdata[j].sky/iEXP[j]*wt;
	    fphot[i].ctcorr+=pdata[j].ctcorr*cwt; // AEDDEBUG was wt
	    //fphot[i].dctcorr+=pdata[j].dctcorr*pdata[j].dctcorr*wt*wt;
	    fphot[i].rnd+=pdata[j].rnd*swt; // AEDDEBUG was wt
	    fphot[i].crowd+=pdata[j].crowd*swt;
	    fphot[i].flag|=pdata[j].flag;
	 }
	 tcm+=iEXP[j];
      }
      if (twt>0.) {
	 fphot[i].ct0/=twt/tcm;
	 fphot[i].ct/=twt/tcm;
	 fphot[i].dct=tcm/sqrt(twt);
	 //fphot[i].chi=sqrt(fphot[i].chi/twt);
	 fphot[i].sky/=twt/tcm;
	 //fphot[i].ctcorr/=twt;
	 //fphot[i].dctcorr=sqrt(fphot[i].dctcorr)/twt;
	 //if (fphot[i].ctcorr>0) fphot[i].m=-2.5*log10(fphot[i].ctcorr);
	 //else fphot[i].m=99.999;
	 //if (fphot[i].ct>0) fphot[i].dm=1.0857362*fphot[i].dct/fphot[i].ct;
      }
      if (tcwt>0.) { // AEDDEBUG new
	 fphot[i].chi=sqrt(fphot[i].chi/tcwt); // AEDDEBUG new
	 fphot[i].ctcorr/=tcwt;
	 fphot[i].dctcorr=1./sqrt(tcwt);
	 if (fphot[i].ctcorr>0) {
	    fphot[i].m=-2.5*log10(fphot[i].ctcorr);
	    fphot[i].dm=1.0857362*fphot[i].dctcorr/fphot[i].ctcorr;
	 }
	 else {
	    fphot[i].m=99.999;
	    fphot[i].dm=9.999;
	 }
      }
      else {
	 fphot[i].dct=9999;
	 fphot[i].dctcorr=9999;
	 fphot[i].m=99.999;
	 fphot[i].dm=9.999;
      }
      if (tswt>0.) {
	 fphot[i].sh/=tswt;
	 fphot[i].rnd/=tswt;
	 fphot[i].crowd/=tswt; // AEDDEBUG was twt
      }
      vmag[i]=fphot[i].m;
      dvmag[i]=fphot[i].dm;
   }
   WFPC2transform(vmag,dvmag,smag);
   for (i=0;i<WFPC2_NFILTERS;i++) if (fused[i]>1) {
      if (fphot[i].ct<999999.5) fprintf(of,"  %8.1f",fphot[i].ct);
      else fprintf(of,"  %8.2e",fphot[i].ct);
      if (fphot[i].sky<99999.95) fprintf(of," %8.2f",fphot[i].sky);
      else if (fphot[i].sky<999999.5) fprintf(of," %8.1f",fphot[i].sky);
      else fprintf(of," %8.1e",fphot[i].sky);

      if (fphot[i].ctcorr<0.0) fprintf(of," %8.1e",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<0.99) fprintf(of," %8.2e",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<9.999995) fprintf(of," %8.6f",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<99.99995) fprintf(of," %8.5f",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<999.9995) fprintf(of," %8.4f",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<9999.995) fprintf(of," %8.3f",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<99999.95) fprintf(of," %8.2f",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<999999.5) fprintf(of," %8.1f",fphot[i].ctcorr);
      else fprintf(of," %8.1e",fphot[i].ctcorr);

      if (fphot[i].dctcorr<0.0) fprintf(of," %8.1e",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<0.99) fprintf(of," %8.2e",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<9.999995) fprintf(of," %8.6f",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<99.99995) fprintf(of," %8.5f",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<999.9995) fprintf(of," %8.4f",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<9999.995) fprintf(of," %8.3f",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<99999.95) fprintf(of," %8.2f",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<999999.5) fprintf(of," %8.1f",fphot[i].dctcorr);
      else fprintf(of," %8.1e",fphot[i].dctcorr);

      fprintf(of," %6.3f %6.3f",fphot[i].m,smag[i]);
      if (fphot[i].dm>9.999) fprintf(of," 9.999");
      else fprintf(of," %5.3f",fphot[i].dm);
      fprintf(of," %6.2f %7.1f %6.3f %6.3f %5.3f %2d",fphot[i].chi,fphot[i].ctcorr/fphot[i].dctcorr,fphot[i].sh,fphot[i].rnd,fphot[i].crowd,fphot[i].flag);
   }
}

void WFPC2outstarimg(int img,FILE *of,float x0,float y0,photdatatype*pdata) {
   float x;
   if (hstmode[img].inst!=WFPC2) {
      printf("Stupid error; called outstarimg for non-WFPC2 data\n");
      exit(-1);
   }
   if (pdata[img].ctcorr<=0) x=pdata[img].m;
   else if (smag[hstmode[img].filt]>99) x=smag[hstmode[img].filt];
   else x=pdata[img].m+smag[hstmode[img].filt]-vmag[hstmode[img].filt];
   fprintf(of," %6.3f %6.3f",pdata[img].m,x);
   if (pdata[img].dm>9.999) fprintf(of," 9.999");
   else fprintf(of," %5.3f",pdata[img].dm);
   return;
}

float wfpc2_apsize(int img,float x,float y) {
   double area;
   if (hstmode[img].inst!=WFPC2) {
      printf("Stupid error; called wfpc2_apsize for non-WFPC2 data\n");
      exit(-1);
   }
   int cm=hstmode[img].cm;
   if (x<0) x=0;
   else if (x>800) x=800;
   if (y<0) y=0;
   else if (y>800) y=800;
   area = WFPC2pixsize(cm,x,y);
   //printf("%d %f %f %f\n",cm,x,y,area);
   return 10.979/sqrt(area);
}

void WFPC2shift(int img,double*x,double*y) {
   double xx,yy;
   if (hstmode[img].inst!=WFPC2) {
      printf("Stupid error; called shift for non-WFPC2 data\n");
      exit(-1);
   }
   xx=*x; yy=*y;
   if (hstmode[img].cm>=0) WFPC2fwddistort(hstmode[img].cm,&xx,&yy);
   *x=xx; *y=yy;
   return;
}

void WFPC2unshift(int img,double*x,double*y) {
   double xx,yy;
   if (hstmode[img].inst!=WFPC2) {
      printf("Stupid error; called wfpc2unshift for non-WFPC2 data\n");
      exit(-1);
   }
   xx=*x; yy=*y;
   if (hstmode[img].cm>=0) WFPC2revdistort(hstmode[img].cm,&xx,&yy);
   *x=xx; *y=yy;
   return;
}

void writewfpc2info(void) {
   int img;
   fprintf(finfo,"* WFPC2-specific info\n");
   for (img=0;img<Nimg;img++) if (hstmode[img].inst==WFPC2) {
      fprintf(finfo,"* image %d: %s %d %f\n",img+1,WFPC2filters[hstmode[img].filt].name,hstmode[img].cm,iEXP[img]);
   }
   return;
}

static int*ffused=NULL;
static double*fakem0;
void WFPC2readfakemag(FILE*f) {
   int i,img;

   if (WFPC2_NFILTERS<0) WFPC2initfilters();
   if (ffused==NULL) {
      ffused=(int*)calloc(sizeof(int),WFPC2_NFILTERS);
      fakem0=(double*)calloc(sizeof(double),WFPC2_NFILTERS);
      if (!ffused || !fakem0) merr();
      for (img=0;img<Nimg;img++) if (hstmode[img].inst==WFPC2) ffused[hstmode[img].filt]++;
   }
   for (i=0;i<WFPC2_NFILTERS;i++) if (ffused[i]) fscanf(f,"%lf",fakem0+i);
   return;
}

void WFPC2fixfakemag(int img,float x0,float y0,double*ct0,float*bg) {
   int i;
   double dm;

   if (hstmode[img].inst!=WFPC2) {
      printf("Stupid error; called wfpc2fixfakemag for non-WFPC2 data\n");
      exit(-1);
   }
   dm=WFPC2calcmag(img,x0,y0,1.0,bg[img],0)-fakem0[hstmode[img].filt];
   ct0[img]=pow(10,0.4*dm);
   for (i=0;i<5;i++) {
      dm=WFPC2calcmag(img,x0,y0,ct0[img],bg[img],WFPC2useCTE)-fakem0[hstmode[img].filt];
      ct0[img]*=pow(10,0.4*dm);
   }
   return;
}

char *WFPC2imagestring(int img) {
   if (hstmode[img].inst!=WFPC2) {
      printf("Stupid error; called wfpc2imagestring for non-WFPC2 data\n");
      exit(-1);
   }
   return WFPC2filters[hstmode[img].filt].name;
}
