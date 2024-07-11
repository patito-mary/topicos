#include "fits.h"
#include <assert.h>
#include "acspsfdata.h"

int fiducialX[9]={0,512,1024,1536,2168,2800,3192,3584,4096};
int fiducialY[10]={0,512,1024,1536,2048,2049,2560,3072,3584,4096};
ftype fits;
double apsf[101][101];
double ipsf[23][23];

void interpAndersonPSF(int x,int y)
{
   int ix,iy,xx,yy;
   double mx,my;
   for (ix=0;ix<7 && x>fiducialX[ix+1];ix++);
   mx=(double)(x-fiducialX[ix])/(double)(fiducialX[ix+1]-fiducialX[ix]);
   for (iy=0;iy<8 && y>fiducialY[iy+1];iy++);
   my=(double)(y-fiducialY[iy])/(double)(fiducialY[iy+1]-fiducialY[iy]);
   for (xx=0;xx<101;xx++) for (yy=0;yy<101;yy++) {
      apsf[yy][xx] = (1-mx)*(1-my)*fits.img.data[0][iy*100+yy][ix*100+xx]
	 + mx*(1-my)*fits.img.data[0][iy*100+yy][ix*100+100+xx]
	 + (1-mx)*my*fits.img.data[0][iy*100+100+yy][ix*100+xx]
	 + mx*my*fits.img.data[0][iy*100+100+yy][ix*100+100+xx];
   }
   //printf("%d %d -> %d %f, %d %f (%f)\n",x,y,ix,mx,iy,my,apsf[50][50]);
}

void cubicSplineWeights(double t,double w[4])
{
   w[0] = -0.5*t + t*t - 0.5*t*t*t;
   w[1] = 1.0 - 2.5*t*t + 1.5*t*t*t;
   w[2] = 0.5*t + 2.0*t*t - 1.5*t*t*t;
   w[3] = -0.5*t*t + 0.5*t*t*t;
}

double cubicSpline(double x[4],double w[4])
{
   return x[0]*w[0] + x[1]*w[1] + x[2]*w[2] + x[3]*w[3];
}

// extrema will be +/- 0.6 in both dx and dy, or +/- 3 in 0.25-sampled space.  Bicubic scheme will require +/-4 shift from center, so this works out to 11 pixels from center
void calcAndersonPSF(double dx,double dy)
{
   int i,x,y,ix,iy;
   double yy,xx;
   double xinterp[4],wy[4],wx[4];

   for (y=-11;y<=11;y++) {
      yy = dy*4.0+4*y+50;
      iy = (int)(yy);
      assert(iy>=1 && iy<=98);
      cubicSplineWeights(yy-iy,wy);
      for (x=-11;x<=11;x++) {
	 xx = dx*4.0+4*x+50;
	 ix = (int)(xx);
	 assert(ix>=1 && ix<=98);
	 cubicSplineWeights(xx-ix,wx);

	 // perform cubic interpolations along X axis first
	 for (i=-1;i<=2;i++) {
	    xinterp[i+1]=cubicSpline(&(apsf[iy+i][ix-1]),wx);
	 }
	 ipsf[y+11][x+11]=cubicSpline(xinterp,wy);
      }
   }
}

void mergepsf(char*filt)
{
   int cm,x0,y0,x,y,x1,y1,t;
   FILE *f,*fout;
   char fn[321];
   float **psf;
   double m,ttt,tand;

   sprintf(fn,"PSFEFF.%s.fits",filt);
   readfits(fn,&fits,0);
   if (fits.Next!=0 || fits.img.X!=901 || fits.img.Y!=1001 || fits.img.Z!=1) {
      printf("Illegal format of %s\n",fn);
      exit(-1);
   }
   for (cm=1;cm<3;cm++)
   {
      // open file, allocate storage
      printf(fn,"%s/acs/data/%s.%s.psf","/Users/patito/Documents/topicos/topicos/ModuloII/",filt,acs_cn[cm]);
      if ((f=fopen(fn,"rb"))==NULL) {
	 printf("Cannot open %s\n",fn);
	 exit(1);
      }
      sprintf(fn,"%s_anderson.%s.psf",filt,acs_cn[cm]);
      if ((fout=fopen(fn,"wb"))==NULL) {
	 printf("Cannot open %s\n",fn);
	 exit(1);
      }
      psf=(float**)calloc(2*acs_rpsf[cm]+1,sizeof(float*));
      if (!psf) merr();
      psf+=acs_rpsf[cm];
      for (y=-acs_rpsf[cm];y<=acs_rpsf[cm];y++) {
	 psf[y]=(float*)calloc(2*acs_rpsf[cm]+1,sizeof(float));
	 if (!psf[y]) merr();
	 psf[y]+=acs_rpsf[cm];
      }

      // loop over positions
      for (y0=0;y0<acs_nypsfpos[cm];y0++) for (x0=0;x0<acs_nxpsfpos[cm];x0++) {
	 if (cm==2) interpAndersonPSF(128+256*x0,128+256*y0);
	 else interpAndersonPSF(128+256*x0,2048+128+256*y0);

	 // loop over phases
	 for (y=-acs_n2psf[cm];y<=acs_n2psf[cm];y++) for (x=-acs_n2psf[cm];x<=acs_n2psf[cm];x++) {
	    for (y1=-acs_rpsf[cm];y1<=acs_rpsf[cm];y1++) ffread(psf[y1]-acs_rpsf[cm],2*acs_rpsf[cm]+1,sizeof(float),f);
	    calcAndersonPSF(-x/(double)acs_sub[cm],-y/(double)acs_sub[cm]);
	    /*
	    for (y1=6;y1>=-6;y1--) {
	       for (x1=-6;x1<=6;x1++) printf("%5d ",(int)(psf[y1][x1]*100000+0.5));
	       printf("\n");
	    }
	    printf("\n");
	    for (y1=17;y1>=5;y1--) {
	       for (x1=5;x1<=17;x1++) printf("%5d ",(int)(ipsf[y1][x1]*100000+0.5));
	       printf("\n");
	    }
	    printf("\n");
	    */
	    ttt=tand=0.0;
	    for (y1=-11;y1<=11;y1++) for (x1=-11;x1<=11;x1++) {
	       t = abs(y1);
	       if (abs(x1)>t) t=abs(x1);
	       if (t<=8) m=1.0;
	       else m=3.0-0.25*t;
	       ttt += m*psf[y1][x1];
	       tand += m*ipsf[11+y1][11+x1];
	    }
	    for (y1=-11;y1<=11;y1++) for (x1=-11;x1<=11;x1++) {
	       t = abs(y1);
	       if (abs(x1)>t) t=abs(x1);
	       if (t<=8) m=1.0;
	       else m=3.0-0.25*t;
	       psf[y1][x1] = (1-m)*psf[y1][x1] + m*ttt/tand*ipsf[11+y1][11+x1];
	    }
	    for (y1=-acs_rpsf[cm];y1<=acs_rpsf[cm];y1++) ffwrite(psf[y1]-acs_rpsf[cm],2*acs_rpsf[cm]+1,sizeof(float),fout);
	    /*
	    printf("%f\n",ttt/tand);
	    for (y1=6;y1>=-6;y1--) {
	       for (x1=-6;x1<=6;x1++) printf("%5d ",(int)(psf[y1][x1]*100000+0.5));
	       printf("\n");
	    }
	    exit(0);
	    */
	  }

      }

      // deallocate storage, close file
      for (y=-acs_rpsf[cm];y<=acs_rpsf[cm];y++) free(psf[y]-acs_rpsf[cm]);
      free(psf-acs_rpsf[cm]);
      fclose(fout);
      fclose(f);
   }
   return;
}

void testmath(void)
{
   int ix;
   double w[4];
   double y0[4]={-1.0,0.0,1.0,2.0};
   double y1[4]={1.0,0.0,1.0,4.0};
   double y2[4]={-6.0,0.0,0.0,6.0};
   for (ix=0;ix<=10;ix++) {
      cubicSplineWeights(0.1*ix,w);
      printf("%2d %f %f %f\n",ix,cubicSpline(y0,w),cubicSpline(y1,w),cubicSpline(y2,w));
   }
}

int main(int argc,char**argv)
{
   int i;

   if (argc<2) {
      printf("Usage: %s <filters>\n",*argv);
      return -1;
   }
   //testmath();
   for (i=1;i<argc;i++) {
      mergepsf(argv[i]);
   }
   return 0;
}
