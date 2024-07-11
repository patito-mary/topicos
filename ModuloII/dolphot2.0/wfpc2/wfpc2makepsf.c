#include <fits.h>
#include "wfpc2psfdata.h"
#include "wfpc2distort.h"

ftype fits;
double dQE=0.00; // using tiny tim kernel instead
double kernel[3][3]=
{
   {0.0125,0.05,0.0125},
   {0.05  ,0.75,0.05  },
   {0.0125,0.05,0.0125}
}; // from tiny tim users' guide

/*
void getpsf(char*filt,int cm) {
   int x,y,z,cx,cy;
   int x1,y1,x2,y2,dx,dy,ct,kx,ky,xx,yy;
   FILE *f;
   char fn[321];
   float **psf,**ttpsf,*xoff;
   double norm,m,xpsf,ypsf,zz;

   sprintf(fn,"%s/wfpc2/data/%s.%s.psf",BASEDIR,filt,wfpc2_cn[cm]);
   if ((f=fopen(fn,"wb"))==NULL) {
      printf("Cannot write %s\n",fn);
      exit(1);
   }
   psf=(float**)calloc(2*wfpc2_rpsf[cm]+1,sizeof(float*));
   if (!psf) merr();
   psf+=wfpc2_rpsf[cm];
   for (z=-wfpc2_rpsf[cm];z<=wfpc2_rpsf[cm];z++) {
      psf[z]=(float*)calloc(2*wfpc2_rpsf[cm]+1,sizeof(float));
      if (!psf[z]) merr();
      psf[z]+=wfpc2_rpsf[cm];
   }
   ttpsf=(float**)calloc(2*wfpc2_rpsf[cm]+1,sizeof(float*));
   if (!ttpsf) merr();
   ttpsf+=wfpc2_rpsf[cm];
   for (z=-wfpc2_rpsf[cm];z<=wfpc2_rpsf[cm];z++) {
      ttpsf[z]=(float*)calloc(2*wfpc2_rpsf[cm]+1,sizeof(float));
      if (!ttpsf[z]) merr();
      ttpsf[z]+=wfpc2_rpsf[cm];
   }
   xoff=(float*)calloc(2*wfpc2_rpsf[cm]+1,sizeof(float));
   if (!xoff) merr();
   xoff+=wfpc2_rpsf[cm];
   // loop over all 64 chip positions
   for (ct=0;ct<wfpc2_npsfpos[cm];ct++) {
      sprintf(fn,"%s/tmp/%s_WFPC2/%s%d%d.fits",BASEDIR,filt,wfpc2_cn[cm],ct/10,ct%10);
      readfits(fn,&fits,0);
      if (fits.Next!=0 || fits.img.Z!=1) {printf("Illegal format of PSF file\n"); exit(-1);}
      xpsf = atof(getcardval(&fits.img,"X_PSF",1));
      ypsf = atof(getcardval(&fits.img,"Y_PSF",1));
      // locate maximum at cx,cy
      cx=cy=0;
      for (y=0;y<fits.img.Y;y++) for (x=0;x<fits.img.X;x++) if (fits.img.data[0][y][x]>fits.img.data[0][cy][cx]) {cx=x; cy=y;}
      if (wfpc2_rpsf[cm]>(cx+1)/wfpc2_sub[cm]-2) {printf("Insufficient room on left; max wfpc2_rpsf[%d]=%d\n",cm,(cx+1)/wfpc2_sub[cm]-2); exit(-1);}
      if (wfpc2_rpsf[cm]>(fits.img.X-cx+1)/wfpc2_sub[cm]-2) {printf("Insufficient room on right; max wfpc2_rpsf[%d]=%d\n",cm,(fits.img.X-cx+1)/wfpc2_sub[cm]-2); exit(-1);}
      if (wfpc2_rpsf[cm]>(cy+1)/wfpc2_sub[cm]-2) {printf("Insufficient room on bottom; max wfpc2_rpsf[%d]=%d\n",cm,(cy+1)/wfpc2_sub[cm]-2); exit(-1);}
      if (wfpc2_rpsf[cm]>(fits.img.Y-cy+1)/wfpc2_sub[cm]-2) {printf("Insufficient room on top; max wfpc2_rpsf[%d]=%d\n",cm,(fits.img.Y-cy+1)/wfpc2_sub[cm]-2); exit(-1);}
      norm=0;
      for (y=0;y<fits.img.Y;y++) for (x=0;x<fits.img.X;x++) norm+=fits.img.data[0][y][x];
      printf("%3d: Center=%d,%d; total PSF=%f\n",ct,cx,cy,norm);
      // loop over sampled PSF phases (x,y)
      for (y=-wfpc2_n2psf[cm];y<=wfpc2_n2psf[cm];y++) for (x=-wfpc2_n2psf[cm];x<=wfpc2_n2psf[cm];x++) {
	 for (y1=-wfpc2_rpsf[cm];y1<=wfpc2_rpsf[cm];y1++) for (x1=-wfpc2_rpsf[cm];x1<=wfpc2_rpsf[cm];x1++) psf[y1][x1]=0;
	 // convert tiny tim PSF into sampled PSF by looping over tiny tim PSF
	 for (y1=0;y1<fits.img.Y;y1++) for (x1=0;x1<fits.img.X;x1++) {
	    dx=x1-cx+x; // distance from center plus phasing
	    dy=y1-cy+y;
	    x2=(dx+wfpc2_sub[cm]/2+1000*wfpc2_sub[cm])/wfpc2_sub[cm]-1000; // pixel where this Tiny Tim pixel falls = round(dx/wfpc2_sub[cm])
	    y2=(dy+wfpc2_sub[cm]/2+1000*wfpc2_sub[cm])/wfpc2_sub[cm]-1000;
	    dx-=x2*wfpc2_sub[cm]; // dx is modulo x2
	    dy-=y2*wfpc2_sub[cm];
	    m=1+dQE*(0.5-3.*(float)(dx*dx+dy*dy)/(float)(wfpc2_sub[cm]*wfpc2_sub[cm])); // apply subpixel QE variations
	    m*=fits.img.data[0][y1][x1];
	    //ignore normalization since TT PSFs are already normalized;
	    //m/=norm;
	    for (kx=-1;kx<=1;kx++) for (ky=-1;ky<=1;ky++) if (abs(x2+kx)<=wfpc2_rpsf[cm] && abs(y2+ky)<=wfpc2_rpsf[cm]) psf[y2+ky][x2+kx]+=m*kernel[ky+1][kx+1];
	 }
	 // apply distortion corrections
	 zz=x/(double)wfpc2_sub[cm];
	 m = WFPC2xsize(cm,xpsf,ypsf);
	 xoff[1]=(1-m)*(0.5-zz);
	 xoff[0]=-(1-m)*(0.5+zz);
	 for (xx=2;xx<=wfpc2_rpsf[cm]+1;xx++) xoff[xx]=xoff[xx-1]+(1-m);
	 for (xx=-1;xx>=-wfpc2_rpsf[cm];xx--) xoff[xx]=xoff[xx+1]-(1-m);
	 for (xx=-wfpc2_rpsf[cm];xx<=wfpc2_rpsf[cm];xx++) for (yy=-wfpc2_rpsf[cm];yy<=wfpc2_rpsf[cm];yy++) ttpsf[xx][yy]=0;
	 for (xx=-wfpc2_rpsf[cm];xx<=wfpc2_rpsf[cm];xx++) for (yy=-wfpc2_rpsf[cm];yy<=wfpc2_rpsf[cm];yy++) {
	    float q;
	    ttpsf[xx][yy]+=(q=psf[xx][yy]);
	    if (xoff[yy]<0) {
	       if (yy>-wfpc2_rpsf[cm]) {
		  ttpsf[xx][yy-1]+=0.5*(q+psf[xx][yy-1])*(-xoff[yy]);
		  ttpsf[xx][yy]-=0.5*(q+psf[xx][yy-1])*(-xoff[yy]);
	       }
	       else ttpsf[xx][yy]-=0.75*q*(-xoff[yy]);
	    }
	    if (xoff[yy+1]>0) {
	       if (yy<wfpc2_rpsf[cm]) {
		  ttpsf[xx][yy+1]+=0.5*(q+psf[xx][yy+1])*xoff[yy+1];
		  ttpsf[xx][yy]-=0.5*(q+psf[xx][yy+1])*xoff[yy+1];
	       }
	       else ttpsf[xx][yy]-=0.75*q*xoff[yy+1];
	    }
	    if (xoff[yy]>0 && yy==-wfpc2_rpsf[cm]) ttpsf[xx][yy]+=0.75*q*xoff[yy];
	    if (xoff[yy+1]<0 && yy==wfpc2_rpsf[cm]) ttpsf[xx][yy]+=0.75*q*(-xoff[yy+1]);
	 }
	 zz=y/(double)wfpc2_sub[cm];
	 m = WFPC2ysize(cm,xpsf,ypsf);
	 xoff[1]=(1-m)*(0.5-zz);
	 xoff[0]=-(1-m)*(0.5+zz);
	 for (yy=2;yy<=wfpc2_rpsf[cm]+1;yy++) xoff[yy]=xoff[yy-1]+(1-m);
	 for (yy=-1;yy>=-wfpc2_rpsf[cm];yy--) xoff[yy]=xoff[yy+1]-(1-m);
	 for (xx=-wfpc2_rpsf[cm];xx<=wfpc2_rpsf[cm];xx++) for (yy=-wfpc2_rpsf[cm];yy<=wfpc2_rpsf[cm];yy++) psf[xx][yy]=0;
	 for (yy=-wfpc2_rpsf[cm];yy<=wfpc2_rpsf[cm];yy++) for (xx=-wfpc2_rpsf[cm];xx<=wfpc2_rpsf[cm];xx++) {
	    float q;
	    psf[xx][yy]+=(q=ttpsf[xx][yy]);
	    if (xoff[xx]<0) {
	       if (xx>-wfpc2_rpsf[cm]) {
		  psf[xx-1][yy]+=0.5*(q+ttpsf[xx-1][yy])*(-xoff[xx]);
		  psf[xx][yy]-=0.5*(q+ttpsf[xx-1][yy])*(-xoff[xx]);
	       }
	       else psf[xx][yy]-=0.75*q*(-xoff[xx]);
	    }
	    if (xoff[xx+1]>0) {
	       if (xx<wfpc2_rpsf[cm]) {
		  psf[xx+1][yy]+=0.5*(q+ttpsf[xx+1][yy])*xoff[xx+1];
		  psf[xx][yy]-=0.5*(q+ttpsf[xx+1][yy])*xoff[xx+1];
	       }
	       else psf[xx][yy]-=0.75*q*xoff[xx+1];
	    }
	    if (xoff[xx]>0 && xx==-wfpc2_rpsf[cm]) psf[xx][yy]+=0.75*q*xoff[xx];
	    if (xoff[xx+1]<0 && xx==wfpc2_rpsf[cm]) psf[xx][yy]+=0.75*q*(-xoff[xx+1]);
	 }
	 // output
	 for (y1=-wfpc2_rpsf[cm];y1<=wfpc2_rpsf[cm];y1++) ffwrite(psf[y1]-wfpc2_rpsf[cm],2*wfpc2_rpsf[cm]+1,sizeof(float),f);
      }
   }
   for (z=-wfpc2_rpsf[cm];z<=wfpc2_rpsf[cm];z++) free(psf[z]-wfpc2_rpsf[cm]);
   free(psf-wfpc2_rpsf[cm]);
   for (z=-wfpc2_rpsf[cm];z<=wfpc2_rpsf[cm];z++) free(ttpsf[z]-wfpc2_rpsf[cm]);
   free(ttpsf-wfpc2_rpsf[cm]);
   free(xoff-wfpc2_rpsf[cm]);
   fclose(f);
   return;
}
*/

void getpsf(char*filt,int cm) {
   FILE *f;
   char fn[321];
   float **psf;
   int x,y,cx,cy,z,x1,y1,ct,kx,ky;
   int ix1,ix2,iy1,iy2,ix,iy;
   double norm,m,xpsf,ypsf,xmag,ymag,dx,dy,dx1,dx2,dy1,dy2;

   sprintf(fn,"%s/wfpc2/data/%s.%s.psf",BASEDIR,filt,wfpc2_cn[cm]);
   if ((f=fopen(fn,"wb"))==NULL) {
      printf("Cannot write %s\n",fn);
      exit(1);
   }
   psf=(float**)calloc(2*wfpc2_rpsf[cm]+1,sizeof(float*));
   if (!psf) merr();
   psf+=wfpc2_rpsf[cm];
   for (z=-wfpc2_rpsf[cm];z<=wfpc2_rpsf[cm];z++) {
      psf[z]=(float*)calloc(2*wfpc2_rpsf[cm]+1,sizeof(float));
      if (!psf[z]) merr();
      psf[z]+=wfpc2_rpsf[cm];
   }
   // loop over all 64 chip positions
   for (ct=0;ct<wfpc2_npsfpos[cm];ct++) {
      sprintf(fn,"%s/tmp/%s_WFPC2/%s%d%d.fits",BASEDIR,filt,wfpc2_cn[cm],ct/10,ct%10);
      readfits(fn,&fits,0);
      if (fits.Next!=0 || fits.img.Z!=1) {printf("Illegal format of PSF file\n"); exit(-1);}
      xpsf = atof(getcardval(&fits.img,"X_PSF",1));
      ypsf = atof(getcardval(&fits.img,"Y_PSF",1));
      xmag = WFPC2xsize(cm,xpsf,ypsf);
      ymag = WFPC2ysize(cm,xpsf,ypsf);
      // locate maximum at cx,cy
      cx=cy=0;
      for (y=0;y<fits.img.Y;y++) for (x=0;x<fits.img.X;x++) if (fits.img.data[0][y][x]>fits.img.data[0][cy][cx]) {cx=x; cy=y;}
      // lines below OK as long as mag <= 1.0 (i.e., 1/mag >= 1.0)
      if (wfpc2_rpsf[cm]>(cx+1)/wfpc2_sub[cm]-2) {printf("Insufficient room on left; max wfpc2_rpsf[%d]=%d\n",cm,(cx+1)/wfpc2_sub[cm]-2); exit(-1);}
      if (wfpc2_rpsf[cm]>(fits.img.X-cx+1)/wfpc2_sub[cm]-2) {printf("Insufficient room on right; max wfpc2_rpsf[%d]=%d\n",cm,(fits.img.X-cx+1)/wfpc2_sub[cm]-2); exit(-1);}
      if (wfpc2_rpsf[cm]>(cy+1)/wfpc2_sub[cm]-2) {printf("Insufficient room on bottom; max wfpc2_rpsf[%d]=%d\n",cm,(cy+1)/wfpc2_sub[cm]-2); exit(-1);}
      if (wfpc2_rpsf[cm]>(fits.img.Y-cy+1)/wfpc2_sub[cm]-2) {printf("Insufficient room on top; max wfpc2_rpsf[%d]=%d\n",cm,(fits.img.Y-cy+1)/wfpc2_sub[cm]-2); exit(-1);}
      norm=0;
      for (y=0;y<fits.img.Y;y++) for (x=0;x<fits.img.X;x++) norm+=fits.img.data[0][y][x];
      printf("%3d (%3d,%3d): Center=%d,%d; total PSF=%f; mag=%5.3f,%5.3f\n",ct,(int)xpsf,(int)ypsf,cx,cy,norm,xmag,ymag);
      // loop over sampled PSF phases (x,y)
      for (y=-wfpc2_n2psf[cm];y<=wfpc2_n2psf[cm];y++) for (x=-wfpc2_n2psf[cm];x<=wfpc2_n2psf[cm];x++) {
	 for (y1=-wfpc2_rpsf[cm];y1<=wfpc2_rpsf[cm];y1++) for (x1=-wfpc2_rpsf[cm];x1<=wfpc2_rpsf[cm];x1++) psf[y1][x1]=0;
	 // convert tiny tim PSF into sampled PSF by looping over tiny tim PSF
	 for (y1=0;y1<fits.img.Y;y1++) for (x1=0;x1<fits.img.X;x1++) {
	    dx1=((x1-cx-0.5)/xmag+x)/wfpc2_sub[cm]; // pixel where this light should be added
	    dx2=((x1-cx+0.5)/xmag+x)/wfpc2_sub[cm];
	    dy1=((y1-cy-0.5)/ymag+y)/wfpc2_sub[cm];
	    dy2=((y1-cy+0.5)/ymag+y)/wfpc2_sub[cm];
	    ix1=round(dx1);
	    ix2=round(dx2);
	    iy1=round(dy1);
	    iy2=round(dy2);
	    for (ix=ix1;ix<=ix2;ix++) for (iy=iy1;iy<=iy2;iy++) {
	       m = 1.0/((dx2-dx1)*(dy2-dy1));
	       if (ix1==ix2) {
		  dx = 0.5*(dx1+dx2);
		  m *= dx2-dx1;
	       }
	       else if (ix==ix1) {
		  m *= ix1+0.5-dx1;
		  dx = 0.5*(dx1+ix1+0.5);
	       }
	       else if (ix==ix2) {
		  m *= dx2-(ix2-0.5);
		  dx = 0.5*(dx2+ix2-0.5);
	       }
	       else {
		  printf("Error: ix1=%d, ix2=%d\n",ix1,ix2);
		  exit(-1);
	       }
	       if (iy1==iy2) {
		  dy = 0.5*(dy1+dy2);
		  m *= dy2-dy1;
	       }
	       else if (iy==iy1) {
		  m *= iy1+0.5-dy1;
		  dy = 0.5*(dy1+iy1+0.5);
	       }
	       else if (iy==iy2) {
		  m *= dy2-(iy2-0.5);
		  dy = 0.5*(dy2+iy2-0.5);
	       }
	       else {
		  printf("Error: iy1=%d, iy2=%d\n",iy1,iy2);
		  exit(-1);
	       }
	       m*=1+dQE*(0.5-3.*(float)(dx*dx+dy*dy)/(float)(wfpc2_sub[cm]*wfpc2_sub[cm])); // apply subpixel QE variations
	       m*=fits.img.data[0][y1][x1];
	       //ignore normalization since TT PSFs are already normalized;
	       //m/=norm;
	       for (kx=-1;kx<=1;kx++) for (ky=-1;ky<=1;ky++) if (abs(ix+kx)<=wfpc2_rpsf[cm] && abs(iy+ky)<=wfpc2_rpsf[cm]) psf[iy+ky][ix+kx]+=m*kernel[ky+1][kx+1];
	    }
	 }
	 // output
	 for (y1=-wfpc2_rpsf[cm];y1<=wfpc2_rpsf[cm];y1++) ffwrite(psf[y1]-wfpc2_rpsf[cm],2*wfpc2_rpsf[cm]+1,sizeof(float),f);
      }
   }
   for (z=-wfpc2_rpsf[cm];z<=wfpc2_rpsf[cm];z++) free(psf[z]-wfpc2_rpsf[cm]);
   free(psf-wfpc2_rpsf[cm]);
   fclose(f);
   return;
}

int main(int argc,char**argv) {
   int usecm[4]={0,0,0,0},nfilt=0,i,j;
   char filt[100][20];

   if (argc<3) {
      printf("Usage: %s <filters> <chips (1-4)>\n",*argv);
      return -1;
   }
   for (i=1;i<argc;i++) {
      if (argv[i][0]=='F') strcpy(filt[nfilt++],argv[i]);
      else if (argv[i][0]>='1' && argv[i][0]<='4' && argv[i][1]==0) usecm[argv[i][0]-'1']=1;
      else {
	 printf("Unknown option: \"%s\"\n",argv[i]);
	 return -1;
      }
   }
   for (i=0;i<nfilt;i++) for (j=0;j<4;j++) if (usecm[j]) {
      getpsf(filt[i],j);
   }
   return 0;
}
