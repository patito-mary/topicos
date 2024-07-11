#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fits.h"
#include "dolphot_defs.h"
#include "acspsfdata.h"
#include "acsfilters.h"
#include "acsdistort.h"

extern void shift(int img,double x0,double y0,double *x,double *y,int dir);

static float********acspsflib=NULL;
static int ACS_RPSF=34;
static int ANY_ACS=0;

void acsinitparam(void) {
   char str[161],*ptr;
   int img;
   ANY_ACS=0;
   for (img=0;img<Timg && !ANY_ACS;img++) {
      strcpy(str,getcardval(dataim+img,"DOL_ACS",0));
      strtol(str,&ptr,10);
      if (ptr!=str) ANY_ACS=1;
   }
   if (ANY_ACS==0) return;
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
      if (RPSF[img]>=acs_rpsf[0]) {RPSF[img]=acs_rpsf[0]-1; printf("Lowering RPSF to %d\n",RPSF[img]);}
      if (RAper[img]>acs_rpsf[0]-1.) {RAper[img]=acs_rpsf[0]-1; printf("Lowering RAper to %d\n",acs_rpsf[0]-1);}
   }
   return;
}

void acsinitpsf(void) {
   int img,i,j,y1,x1,y2,x2,y3,n2,n3,n3skip;
   char str[161],*ptr;
   FILE *f;
   float***** ptr1;
   float**** ptr2;
   float*** ptr3;
   float** ptr4;
   float* ptr5;

   if (ACS_NFILTERS<0) ACSinitfilters();
   if (acspsflib==NULL) {
      acspsflib=(float********)calloc(sizeof(float*******),3);
      if (!acspsflib) merr();
      for (i=0;i<3;i++) {
	 acspsflib[i]=(float*******)calloc(sizeof(float******),ACS_NFILTERS);
	 if (!acspsflib[i]) merr();
      }
   }
   for (i=0;i<3;i++) for (j=0;j<ACS_NFILTERS;j++) acspsflib[i][j]=NULL;
   for (img=0;img<Timg;img++) {
      strcpy(str,getcardval(dataim+img,"DOL_ACS",0));
      int cm=strtol(str,&ptr,10);
      if (ptr==str) {
	 // AEDDEBUG need warning: printf("**Image %d has not been preprocessed with acsmask; cannot proceed\n",img+1);
      }
      else {
	 if (cm<-2 || cm>2) {
	    printf("**Image %d's chip cannot be identified; please report bug.\n",img+1);
	    exit(-1);
	 }
	 hstmode[img].inst=ACS;
	 hstmode[img].cm=cm;
	 if (img==Nimg) {
	    if (hstmode[img].cm==-2) {
	       hstmode[img].cm=1;
	       DRIZZLE_BASE=1;
	    }
	    else if (hstmode[img].cm==-1) {
	       hstmode[img].cm=0;
	       DRIZZLE_BASE=1;
	    }
	    else DRIZZLE_BASE=0;
	 }
	 if (hstmode[img].cm<0) {
	    printf("**Image %d has been drizzled; cannot run photometry\n",img+1);
	    exit(-1);
	 }
	 strcpy(str,getcardval(datahd+img,"FILTER1",0));
	 if (hstmode[img].cm>0 && !strcmp(str,"CLEAR1L")) strcpy(str,getcardval(datahd+img,"FILTER2",0));
	 else if (hstmode[img].cm==0 && !strcmp(str,"CLEAR1S")) strcpy(str,getcardval(datahd+img,"FILTER2",0));
	 hstmode[img].filt=ACSfindfilt(str);
	 //printf("Image %d: cm=%d, filt=%d (%s)\n",img+1,hstmode[img].cm,hstmode[img].filt,ACSfilters[hstmode[img].filt].name);
	 if (RPSF[img]>=acs_rpsf[hstmode[img].cm]) {printf("ERROR: RPSF must be less than %d for ACS/WFC data\n",acs_rpsf[hstmode[img].cm]); exit(-1);}
	 if (RAper[img]>acs_rpsf[hstmode[img].cm]-1.) {printf("ERROR: RAper must be no more than %d for ACS/WFC data\n",acs_rpsf[hstmode[img].cm]-1); exit(-1);}
      }
   }
   ACS_RPSF=0;
   for (img=0;img<Timg;img++) if (hstmode[img].inst==ACS) {
      i=hstmode[img].cm;
      apsf[img][0][0]=1.;
      apsf[img][1][0]=1.;
      apsf[img][2][0]=0.;
      if (i==0) apsize[img]=11.;
      else apsize[img]=10.;
      for (i=1;i<5;i++) apsf[img][0][i]=apsf[img][1][i]=apsf[img][2][i]=0.;
      if (ACS_RPSF<RPSF[img]) ACS_RPSF=RPSF[img];
      if (ACS_RPSF<rphot[img]) ACS_RPSF=rphot[img];
   }
   for (img=0;img<Timg;img++) if (hstmode[img].inst==ACS && acspsflib[i=hstmode[img].cm][j=hstmode[img].filt]==NULL) {
      f=0;
      if (ACSpsfType==1) {
	 sprintf(str,"%s/acs/data/%s_anderson.%s.psf","/Users/patito/Documents/topicos/topicos/ModuloII",ACSfilters[hstmode[img].filt].name,acs_cn[i]);
	 f=fopen(str,"rb");
      }
      if (f==0) {
	 sprintf(str,"%s/acs/data/%s.%s.psf","/Users/patito/Documents/topicos/topicos/ModuloII",ACSfilters[hstmode[img].filt].name,acs_cn[i]);
	 if ((f=fopen(str,"rb"))==NULL) {
	    printf("Cannot open %s\n",str);
	    exit(-1);
	 }
      }
      n2=2*acs_n2psf[i]+1;
      n3=2*ACS_RPSF+1;
      n3skip=acs_rpsf[i]-ACS_RPSF;
      acspsflib[i][j]=(float******)calloc(sizeof(float*****),acs_nypsfpos[i]);
      ptr1=(float*****)calloc(sizeof(float****),acs_nypsfpos[i]*acs_nxpsfpos[i]);
      ptr2=(float****)calloc(sizeof(float***),acs_nypsfpos[i]*acs_nxpsfpos[i]*n2);
      ptr3=(float***)calloc(sizeof(float**),acs_nypsfpos[i]*acs_nxpsfpos[i]*n2*n2);
      ptr4=(float**)calloc(sizeof(float*),acs_nypsfpos[i]*acs_nxpsfpos[i]*n2*n2*n3);
      ptr5=(float*)calloc(sizeof(float),acs_nypsfpos[i]*acs_nxpsfpos[i]*n2*n2*n3*n3);
      if (!acspsflib[i][j] || !ptr1) merr();
      for (y1=0;y1<acs_nypsfpos[i];y1++) {
	 acspsflib[i][j][y1]=ptr1;
	 ptr1+=acs_nxpsfpos[i];
	 for (x1=0;x1<acs_nxpsfpos[i];x1++) {
	    acspsflib[i][j][y1][x1]=ptr2+acs_n2psf[i];
	    ptr2+=n2;
	    for (y2=-acs_n2psf[i];y2<=acs_n2psf[i];y2++) {
	       acspsflib[i][j][y1][x1][y2]=ptr3+acs_n2psf[i];
	       ptr3+=n2;
	       for (x2=-acs_n2psf[i];x2<=acs_n2psf[i];x2++) {
		  acspsflib[i][j][y1][x1][y2][x2]=ptr4+ACS_RPSF;
		  ptr4+=n3;
		  if (n3skip) fseek(f,4*n3skip*(2*acs_rpsf[i]+1),SEEK_CUR);
		  for (y3=-ACS_RPSF;y3<=ACS_RPSF;y3++) {
		     acspsflib[i][j][y1][x1][y2][x2][y3]=ptr5+ACS_RPSF;
		     ptr5+=n3;
		     if (n3skip) fseek(f,4*n3skip,SEEK_CUR);
		     ffread(acspsflib[i][j][y1][x1][y2][x2][y3]-ACS_RPSF,4,n3,f);
		     if (n3skip) fseek(f,4*n3skip,SEEK_CUR);
		  }
		  if (n3skip) fseek(f,4*n3skip*(2*acs_rpsf[i]+1),SEEK_CUR);
	       }
	    }
	 }
      }
      fclose(f);
   }
   return;
}

void acsfreepsf(void) {
   int i,j;
   for (i=0;i<3;i++) for (j=0;j<ACS_NFILTERS;j++) if (acspsflib[i][j]) {
      free(acspsflib[i][j][0][0][-acs_n2psf[i]][-acs_n2psf[i]][-ACS_RPSF]-ACS_RPSF);
      free(acspsflib[i][j][0][0][-acs_n2psf[i]][-acs_n2psf[i]]-ACS_RPSF);
      free(acspsflib[i][j][0][0][-acs_n2psf[i]]-acs_n2psf[i]);
      free(acspsflib[i][j][0][0]-acs_n2psf[i]);
      free(acspsflib[i][j][0]);
      free(acspsflib[i][j]);
      acspsflib[i][j]=NULL;
   }
   return;
}

//#define DEBUG_PSF
int calcacspsf(int img,float x,float y,int r,int force) {
   int i,j,y1,x1,y2,x2,yy,xx;
   float mx1=1,my1=1,imx1=0,imy1=0;
   float mx2,my2,imx2,imy2;
   static int first=1,lastr=0,lastimg=0;
   static float lastx=0,lasty=0;

   if (hstmode[img].inst!=ACS) {
      printf("Stupid error; called acspsf for non-ACS data\n");
      exit(-1);
   }
   if (!first && lastpsftype==1 && img==lastimg && x==lastx && y==lasty && r<=lastr && !poffreset && !force) return 0;
   first=0;lastpsftype=1;lastimg=img;lastx=x;lasty=y;lastr=r;
   i=hstmode[img].cm;
   j=hstmode[img].filt;
#ifdef DEBUG_PSF
   printf("%d %d %d %f %f %d",img,i,j,x,y,r);
   printf("%s/%s PSF at %f,%f; r=%d:\n",acs_cn[i],ACSfilters[j].name,x,y,r);
   fflush(stdout);
#endif
   imy2=y-(int)y;
   if (imy2<0) imy2++;
   imy2=(imy2-0.5)*acs_sub[i];
   y2=(int)(imy2+acs_sub[i])-acs_sub[i];
   imy2-=y2; my2=1-imy2;
   imx2=x-(int)x;
   if (imx2<0) imx2++;
   imx2=(imx2-0.5)*acs_sub[i];
   x2=(int)(imx2+acs_sub[i])-acs_sub[i];
   imx2-=x2; mx2=1-imx2;
   if (InterpPSFlib && (img<Nimg || !DRIZZLE_BASE)) {
      imx1 = x/256.0 - 0.5;
      imy1 = y/256.0 - 0.5;
      if (imx1<=0.0) {x1=0; imx1=0.0;}
      else if (imx1>=acs_nxpsfpos[i]-1.0) {x1=acs_nxpsfpos[i]-2; imx1=1.0;}
      else {x1=(int)imx1; imx1-=x1;}
      mx1 = 1.0-imx1;
      if (imy1<=0.0) {y1=0; imy1=0.0;}
      else if (imy1>=acs_nypsfpos[i]-1.0) {y1=acs_nypsfpos[i]-2; imy1=1.0;}
      else {y1=(int)imy1; imy1-=y1;}
      my1 = 1.0-imy1;
      for (yy=-r;yy<=r;yy++) for (xx=-r;xx<=r;xx++) {
	 psf[yy][xx] = 
	    ( acspsflib[i][j][y1][x1][y2][x2][yy][xx]*mx2*my2+acspsflib[i][j][y1][x1][y2][x2+1][yy][xx]*imx2*my2+acspsflib[i][j][y1][x1][y2+1][x2][yy][xx]*mx2*imy2+acspsflib[i][j][y1][x1][y2+1][x2+1][yy][xx]*imx2*imy2 ) * mx1*my1 +
	    ( acspsflib[i][j][y1+1][x1][y2][x2][yy][xx]*mx2*my2+acspsflib[i][j][y1+1][x1][y2][x2+1][yy][xx]*imx2*my2+acspsflib[i][j][y1+1][x1][y2+1][x2][yy][xx]*mx2*imy2+acspsflib[i][j][y1+1][x1][y2+1][x2+1][yy][xx]*imx2*imy2 ) * mx1*imy1 +
	    ( acspsflib[i][j][y1][x1+1][y2][x2][yy][xx]*mx2*my2+acspsflib[i][j][y1][x1+1][y2][x2+1][yy][xx]*imx2*my2+acspsflib[i][j][y1][x1+1][y2+1][x2][yy][xx]*mx2*imy2+acspsflib[i][j][y1][x1+1][y2+1][x2+1][yy][xx]*imx2*imy2 ) * imx1*my1 +
	    ( acspsflib[i][j][y1+1][x1+1][y2][x2][yy][xx]*mx2*my2+acspsflib[i][j][y1+1][x1+1][y2][x2+1][yy][xx]*imx2*my2+acspsflib[i][j][y1+1][x1+1][y2+1][x2][yy][xx]*mx2*imy2+acspsflib[i][j][y1+1][x1+1][y2+1][x2+1][yy][xx]*imx2*imy2 ) * imx1*imy1;
      }
   }
   else {
      if (img==Nimg && DRIZZLE_BASE) {
	 y1=acs_nypsfpos[i]/2-1;
	 x1=acs_nxpsfpos[i]/2-1;
      }
      else {
	 y1=(int)y/256; if (y1<0) y1=0; if (y1>=acs_nypsfpos[i]) y1=acs_nypsfpos[i]-1;
	 x1=(int)x/256; if (x1<0) x1=0; if (x1>=acs_nxpsfpos[i]) x1=acs_nxpsfpos[i]-1;
      }
      for (yy=-r;yy<=r;yy++) for (xx=-r;xx<=r;xx++) psf[yy][xx]=acspsflib[i][j][y1][x1][y2][x2][yy][xx]*mx2*my2+acspsflib[i][j][y1][x1][y2][x2+1][yy][xx]*imx2*my2+acspsflib[i][j][y1][x1][y2+1][x2][yy][xx]*mx2*imy2+acspsflib[i][j][y1][x1][y2+1][x2+1][yy][xx]*imx2*imy2;
   }
#ifdef DEBUG_PSF
   printf("y1=%d, my1=%f; x1=%d, mx1=%f\n",y1,my1,x1,mx1);
   printf("y2=%d, my2=%f; x2=%d, mx2=%f\n",y2,my2,x2,mx2);
   for (yy=6;yy>=-6;yy--) if (yy>=-r && yy<=r) {for (xx=-6;xx<=6;xx++) if (xx>=-r && xx<=r) printf("%5d ",(int)(psf[yy][xx]*100000+0.5)); printf("\n");}
   fflush(stdout);
#endif
   return 1;
}

#ifdef DEBUG_PSF
void ACSdumpPSFs(void) {
   int img,xphase,yphase,x,y;
   for (img=0;img<Nimg;img++) if (hstmode[img].inst==ACS) {
      if (hstmode[img].cm==0) {
	 for (xphase=-5;xphase<=5;xphase++) calcacspsf(img,448.5+0.1*xphase,448.5,RPSF[img],1);
	 for (yphase=-5;yphase<=5;yphase++) calcacspsf(img,448.5,448.5+0.1*yphase,RPSF[img],1);
	 for (x=0;x<=1024;x+=32) calcacspsf(img,0.5+x,448.5,RPSF[img],1);
	 for (y=0;y<=1024;y+=32) calcacspsf(img,448.5,0.5+y,RPSF[img],1);
      }
      else {
	 for (xphase=-5;xphase<=5;xphase++) calcacspsf(img,1920.5+0.1*xphase,896.5,RPSF[img],1);
	 for (yphase=-5;yphase<=5;yphase++) calcacspsf(img,1920.5,896.5+0.1*yphase,RPSF[img],1);
	 for (x=0;x<=4096;x+=64) calcacspsf(img,0.5+x,896.5,RPSF[img],1);
	 for (y=0;y<=2048;y+=64) calcacspsf(img,1920.5,0.5+y,RPSF[img],1);
      }
   }
}
#endif

double ACScalcmag(int img,float x0,float y0,float ct0,float bg,int useCTE) {
   float cm=1.;
   double x,y,m;

   if (hstmode[img].inst!=ACS) {
      printf("Stupid error; called acscalcmag for non-ACS data\n");
      exit(-1);
   }
   if (ACS_NFILTERS<0) ACSinitfilters();
   m=-2.5*log10(ct0*apcor[img]/iEXP[img]*acs_ctmult[hstmode[img].cm])+ACS_ZP(hstmode[img].filt,hstmode[img].cm,iEPOCH[img]);
   if (useCTE) {
      shift(img,x0,y0,&x,&y,1);
      if (iEXP0[img]>0.) cm=iEXP[img]/iEXP0[img];
      m-=ACS_CTE(hstmode[img].cm,x,y,ct0,cm,iGAIN[img],bg,iEPOCH[img]);
   }
   return m;
}

void ACSoutstarinfo(FILE *f,int*ct) {
   int i,j,n;
   for (i=0;i<ACS_NFILTERS;i++) {
      n=0;
      for (j=0;j<Nimg;j++) if (hstmode[j].inst==ACS && hstmode[j].filt==i) n++;
      if (n>1) {
	 fprintf(f,"%d. Total counts, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Total sky level, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Normalized count rate, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Normalized count rate uncertainty, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Instrumental VEGAMAG magnitude, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Transformed UBVRI magnitude, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Magnitude uncertainty, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Chi, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Signal-to-noise, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Sharpness, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Roundness, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Crowding, ACS_%s\n",++(*ct),ACSfilters[i].name);
	 fprintf(f,"%d. Photometry quality flag, ACS_%s\n",++(*ct),ACSfilters[i].name);
      }
   }
}

static double*smag=NULL,*vmag=NULL,*dvmag=NULL;
void ACSoutstar(FILE *of,float x0,float y0,photdatatype*pdata) {
   int i,j;
   static int *fused=NULL;
   static photdatatype*fphot=NULL;
   int cm=-1;
   double x,y,m=1.0;

   if (ACS_NFILTERS<0) ACSinitfilters();
   if (!fphot) {
      fused=(int*)calloc(sizeof(int),ACS_NFILTERS);
      fphot=(photdatatype*)calloc(sizeof(photdatatype),ACS_NFILTERS);
      smag=(double*)calloc(sizeof(double),ACS_NFILTERS);
      vmag=(double*)calloc(sizeof(double),ACS_NFILTERS);
      dvmag=(double*)calloc(sizeof(double),ACS_NFILTERS);
      if (!fused || !fphot || !smag || !vmag || !dvmag) merr();
   }
   cm=-1;
   for (j=0;j<Nimg;j++) if (hstmode[j].inst==ACS) {
      if (cm==-1) cm=(1+hstmode[j].cm)/2;
      else if (cm>=0 && cm!=(1+hstmode[j].cm)/2) {
	 printf("Multi-mode transformations not yet supported; providing only VEGAMAG\n");
	 cm=-2;
      }

      double dm = -2.5*log10(acs_ctmult[hstmode[j].cm]) + ACS_ZP(hstmode[j].filt,hstmode[j].cm,iEPOCH[j])-Zero;
      if (pdata[j].ct>0) {
	 if (ACSuseCTE) {
	    shift(j,x0,y0,&x,&y,1);
	    if (iEXP0[j]>0.) m=iEXP[j]/iEXP0[j];
	    else m=1.0;
	    dm -= ACS_CTE(hstmode[j].cm,x,y,pdata[j].ct,m,iGAIN[j],pdata[j].sky,iEPOCH[j]);
	 }
	 pdata[j].m += dm;
      }
      pdata[j].ctcorr *= pow(10,-0.4*dm);
      pdata[j].dctcorr *= pow(10,-0.4*dm);
   }
   for (i=0;i<ACS_NFILTERS;i++) {
      double wt,cwt,swt,twt=0.,tcwt=0.,tswt=0.,is,iss,tcm=0.;
      fphot[i].ct0=fphot[i].ct=fphot[i].chi=fphot[i].sh=fphot[i].sky=fphot[i].ctcorr=fphot[i].dctcorr=fphot[i].rnd=fphot[i].crowd=0.;
      fused[i]=0;
      fphot[i].flag=0;
      for (j=0;j<Nimg;j++) if (hstmode[j].inst==ACS && hstmode[j].filt==i) {
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
   ACStransform(cm,vmag,dvmag,smag);
   for (i=0;i<ACS_NFILTERS;i++) if (fused[i]>1) {
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

void ACSoutstarimg(int img,FILE *of,float x0,float y0,photdatatype*pdata) {
   float x;
   if (hstmode[img].inst!=ACS) {
      printf("Stupid error; called outstarimg for non-ACS data\n");
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

double ACSidcfwd0[3][25]={
   {531,512,205.1358,470.1861,-0.07,0.025,0.000030,0.028415,-3.998E-08,2.373E-07,-1.047E-07,9.507E-12,-2.084E-11,6.676E-12,-2.627E-11,2.469E-02,2.804E-03,2.771E-07,-5.868E-08,3.675E-08,6.095E-11,-1.247E-11,3.535E-11,2.432E-11,1.123113},
   {2072,1024,256.4158,237.0618,-2.38,0.050,0.002032,0.049160,1.029E-07,-3.527E-07,4.240E-07,4.461E-12,-2.258E-11,-5.327E-12,-2.187E-11,4.895E-02,2.231E-03,-4.904E-07,2.955E-07,-1.304E-07,-1.571E-11,-1.996E-12,-1.846E-11,3.686E-12,0.9809135},
   {2072,1024,254.2550,301.1133,-1.94,0.050,0.001718,0.049779,9.929E-08,-2.541E-07,4.316E-07,-1.055E-13,-2.606E-11,-2.602E-12,-2.286E-11,5.068E-02,1.583E-03,-3.677E-07,3.068E-07,-8.124E-08,-1.960E-11,-4.799E-12,-2.736E-11,3.714E-12,0.9809135}
};

float acs_apsize(int img,float x,float y) {
   float xc,yc,p2x,p2y,p2,p1x,p1y,p1,th,area;
   int cm=hstmode[img].cm;
   if (hstmode[img].inst!=ACS) {
      printf("Stupid error; called acs_apsize for non-ACS data\n");
      exit(-1);
   }
   if (x<0) x=0;
   if (y<0) y=0;
   if (cm==0) {
      if (x>1024) x=1024;
      if (y>1024) y=1024;
   }
   else {
      if (x>4096) x=4096;
      if (y>2048) y=2048;
   }
   //compute pixel areas;
   xc=x-ACSidcfwd0[cm][0];
   yc=y-ACSidcfwd0[cm][1];
   p2x=(ACSidcfwd0[cm][7]+ACSidcfwd0[cm][9]*yc+2*ACSidcfwd0[cm][10]*xc+ACSidcfwd0[cm][12]*yc*yc+2*ACSidcfwd0[cm][13]*xc*yc+3*ACSidcfwd0[cm][14]*xc*xc)/ACSidcfwd0[cm][5]; //dX/dx
   p2y=(ACSidcfwd0[cm][6]+2*ACSidcfwd0[cm][8]*yc+ACSidcfwd0[cm][9]*xc+3*ACSidcfwd0[cm][11]*yc*yc+2*ACSidcfwd0[cm][12]*xc*yc+ACSidcfwd0[cm][13]*xc*xc)/ACSidcfwd0[cm][5]; //dY/dx
   p2=sqrt(p2x*p2x+p2y*p2y);
   p1x=(ACSidcfwd0[cm][16]+ACSidcfwd0[cm][18]*yc+2*ACSidcfwd0[cm][19]*xc+ACSidcfwd0[cm][21]*yc*yc+2*ACSidcfwd0[cm][22]*xc*yc+3*ACSidcfwd0[cm][23]*xc*xc)/ACSidcfwd0[cm][5]; //dX/dy
   p1y=(ACSidcfwd0[cm][15]+2*ACSidcfwd0[cm][17]*yc+ACSidcfwd0[cm][18]*xc+3*ACSidcfwd0[cm][20]*yc*yc+2*ACSidcfwd0[cm][21]*xc*yc+ACSidcfwd0[cm][22]*xc*xc)/ACSidcfwd0[cm][5]; //dY/dy
   p1=sqrt(p1x*p1x+p1y*p1y);
   th=acos((p1x*p2x+p1y*p2y)/p1/p2);
   area=p1*p2*sin(th);
   if (cm==0) return 12./sqrt(area);
   return 10./sqrt(area);
}

void ACSshift(int img,double*x,double*y) {
   double xx,yy;
   if (hstmode[img].inst!=ACS) {
      printf("Stupid error; called acsshift for non-ACS data\n");
      exit(-1);
   }
   xx=*x; yy=*y;
   if (hstmode[img].cm>=0) ACSfwddistort(hstmode[img].cm,hstmode[img].filt,&xx,&yy);
   *x=xx; *y=yy;
   return;
}

void ACSunshift(int img,double*x,double*y) {
   double xx,yy;
   if (hstmode[img].inst!=ACS) {
      printf("Stupid error; called acsunshift for non-ACS data\n");
      exit(-1);
   }
   xx=*x; yy=*y;
   if (hstmode[img].cm>=0) ACSrevdistort(hstmode[img].cm,hstmode[img].filt,&xx,&yy);
   *x=xx; *y=yy;
   return;
}

void writeacsinfo(void) {
   int img;
   fprintf(finfo,"* ACS-specific info\n");
   for (img=0;img<Nimg;img++) if (hstmode[img].inst==ACS) {
      fprintf(finfo,"* image %d: %s %d %f\n",img+1,ACSfilters[hstmode[img].filt].name,hstmode[img].cm,iEXP[img]);
      if (fabs(iGAIN[img]-1.)>0.001) printf("ERROR: All ACS data should have gain of 1\n");
   }
   return;
}

static int*ffused=NULL;
static double*fakem0;
void ACSreadfakemag(FILE*f) {
   int i,img;

   if (ACS_NFILTERS<0) ACSinitfilters();
   if (ffused==NULL) {
      ffused=(int*)calloc(sizeof(int),ACS_NFILTERS);
      fakem0=(double*)calloc(sizeof(double),ACS_NFILTERS);
      if (!ffused || !fakem0) merr();
      for (img=0;img<Nimg;img++) if (hstmode[img].inst==ACS) ffused[hstmode[img].filt]++;
   }
   for (i=0;i<ACS_NFILTERS;i++) if (ffused[i]) fscanf(f,"%lf",fakem0+i);
   return;
}

void ACSfixfakemag(int img,float x0,float y0,double*ct0,float*bg) {
   int i;
   double dm;

   if (hstmode[img].inst!=ACS) {
      printf("Stupid error; called acsfixfakemag for non-ACS data\n");
      exit(-1);
   }
   dm=ACScalcmag(img,x0,y0,1.0,bg[img],0)-fakem0[hstmode[img].filt];
   ct0[img]=pow(10,0.4*dm);
   for (i=0;i<5;i++) {
      dm=ACScalcmag(img,x0,y0,ct0[img],bg[img],ACSuseCTE)-fakem0[hstmode[img].filt];
      ct0[img]*=pow(10,0.4*dm);
   }
   return;
}

char *ACSimagestring(int img) {
   if (hstmode[img].inst!=ACS) {
      printf("Stupid error; called acsimagestring for non-ACS data\n");
      exit(-1);
   }
   return ACSfilters[hstmode[img].filt].name;
}
