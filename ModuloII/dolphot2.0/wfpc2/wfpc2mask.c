#include "../include/fits.h"
#include "wfpc2psfdata.h"
#include "wfpc2distort.h"

ftype fits,dq;
double GAIN,RN,EXP,EXP0,EPOCH,FIXED_ET=-1.;
float DMIN,DMAX;
int MASKCR=1;

int WFPC2type(ftype *f) {
   int i;
   if (f->img.X==800 && f->img.Y==800 && f->img.Z==4) return 1;
   if (f->Next>=4) {
      for (i=0;i<4 && f->ext[i].X==800 && f->ext[i].Y==800 && f->ext[i].Z==1;i++);
      if (i==4) {
	 int ext;
	 for (ext=4;ext<fits.Next;ext++) freeim(fits.ext+ext);
	 fits.Next=4;
	 return 2;
      }
   }
   if (f->Next>=3 && f->ext[0].Z==1) {
      for (i=1;i<3 && f->ext[i].X==f->ext[0].X && f->ext[i].Y==f->ext[0].Y && (f->ext[i].Z==1 || (i==2 && f->ext[i].Z>1));i++);
      if (i==3) {
	 printf("Irregular size; assuming drizzled\n");
	 return -1;
      }
   }
   return 0;
}

void WFPC2getcards(int drz) {
   if (strcmp(getcardval(&(fits.img),"FILETYPE",1),"SCI") || strcmp(getcardval(&(fits.img),"TELESCOP",1),"HST") || strcmp(getcardval(&(fits.img),"INSTRUME",1),"WFPC2")) {printf("**Format error (filetype,telescop,instrume)\n"); exit(-1);}
   GAIN=atof(getcardval(&(fits.img),"ATODGAIN",1));
   if (GAIN==7) RN=5.0;
   else if (GAIN==14) RN=7.5;
   else if (GAIN==15) {
      RN=7.5;
      GAIN=14;
   }
   else {printf("**Format error (atodgain)\n"); exit(-1);}
   if (FIXED_ET>0) EXP=FIXED_ET;
   else EXP=atof(getcardval(&(fits.img),"EXPTIME",1));
   DMIN=safedown(atof(getcardval(&(fits.img),"RSDPFILL",1)));
   DMAX=safeup(atof(getcardval(&(fits.img),"SATURATE",1)));
   EPOCH=0.5*(atof(getcardval(&(fits.img),"EXPSTART",1))+atof(getcardval(&(fits.img),"EXPEND",1)));
   EXP0=EXP;
   return;
}

void WFPC2toextensions(void) {
   int ext,i;
   if (fits.Next>0) {
      for (ext=0;ext<fits.Next;ext++) freeim(fits.ext+ext);
      free(fits.ext);
   }
   fits.Next=4;
   fits.ext=(imtype*)calloc(sizeof(imtype),fits.Next);
   if (!fits.ext) merr();
   for (ext=0;ext<4;ext++) {
      fits.ext[ext].Ncards=0;
      fits.ext[ext].Nmax=0;
      fits.ext[ext].X=fits.ext[ext].Y=800;
      fits.ext[ext].Z=1;
      fits.ext[ext].bits=-32;
      fits.ext[ext].pcount=0;
      fits.ext[ext].bzero=0.;
      fits.ext[ext].bscale=1.;
      strcpy(fits.ext[ext].xtension,"IMAGE");
      fits.ext[ext].cards=NULL;
      fits.ext[ext].data=allocimg(800,800,1);
      memcpy(fits.ext[ext].data[0][0],fits.img.data[ext][0],sizeof(float)*800*800);
      if (!strcmp(getcardval(&(fits.img),"CTYPE1",0),"RA---TAN") && !strcmp(getcardval(&(fits.img),"CTYPE2",0),"DEC--TAN")) {
	 double xref0,yref0,xref,yref,x,y,x0,y0,m0[2][2],m1[2][2],m[2][2];
	 xref0 = xref = atof(getcardval(&(fits.img),"CRPIX1",1))-0.5;
	 yref0 = yref = atof(getcardval(&(fits.img),"CRPIX2",1))-0.5;
	 WFPC2fwddistort(0,&xref,&yref);
	 WFPC2revdistort(ext,&xref,&yref);
	 x0=xref; y0=yref; WFPC2fwddistort(ext,&x0,&y0); WFPC2revdistort(0,&x0,&y0);
	 for (i=0;i<5;i++) {
	    switch(ext) {
	    case 0:
	       xref+=xref0-x0;
	       yref+=yref0-y0;
	       break;
	    case 1:
	       xref+=0.457*(yref0-y0);
	       yref+=-0.457*(xref0-x0);
	       break;
	    case 2:
	       xref+=-0.457*(xref0-x0);
	       yref+=-0.457*(yref0-y0);
	       break;
	    case 3:
	       xref+=-0.457*(yref0-y0);
	       yref+=0.457*(xref0-x0);
	       break;
	    }
	    x0=xref; y0=yref; WFPC2fwddistort(ext,&x0,&y0); WFPC2revdistort(0,&x0,&y0);
	 }
	 x=xref+1; y=yref; WFPC2fwddistort(ext,&x,&y); WFPC2revdistort(0,&x,&y);
	 m1[0][0] = x-x0; // dxPC/dx
	 m1[1][0] = y-y0; // dyPC/dx
	 x=xref; y=yref+1; WFPC2fwddistort(ext,&x,&y); WFPC2revdistort(0,&x,&y);
	 m1[0][1] = x-x0; // dxPC/dy
	 m1[1][1] = y-y0; // dyPC/dy
	 m0[0][0] = atof(getcardval(&(fits.img),"CD1_1",1)); // dRA/dxPC
	 m0[0][1] = atof(getcardval(&(fits.img),"CD1_2",1)); // dRA/dyPC
	 m0[1][0] = atof(getcardval(&(fits.img),"CD2_1",1)); // dDec/dxPC
	 m0[1][1] = atof(getcardval(&(fits.img),"CD2_2",1)); // dDec/dyPC
	 m[0][0] = m0[0][0]*m1[0][0] + m0[0][1]*m1[1][0]; // dRA/dx
	 m[0][1] = m0[0][0]*m1[0][1] + m0[0][1]*m1[1][1]; // dRA/dy
	 m[1][0] = m0[1][0]*m1[0][0] + m0[1][1]*m1[1][0]; // dDec/dx
	 m[1][1] = m0[1][0]*m1[0][1] + m0[1][1]*m1[1][1]; // dDec/dy
	 char card[81];
	 if (ext==0) printf("Converting WCS\n");
	 addcard(fits.ext+ext,"CTYPE1  =             RA---TAN /                                                ");
	 addcard(fits.ext+ext,"CTYPE2  =             DEC--TAN /                                                ");
	 sprintf(card,"CRPIX1  = %20.4f /                                                ",xref+0.5);
	 addcard(fits.ext+ext,card);
	 sprintf(card,"CRPIX2  = %20.4f /                                                ",yref+0.5);
	 addcard(fits.ext+ext,card);
	 sprintf(card,"CRVAL1  = %20s /                                                ",getcardval(&(fits.img),"CRVAL1",1));
	 addcard(fits.ext+ext,card);
	 sprintf(card,"CRVAL2  = %20s /                                                ",getcardval(&(fits.img),"CRVAL2",1));
	 addcard(fits.ext+ext,card);
	 sprintf(card,"CD1_1   = %20.10e /                                                ",m[0][0]);
	 addcard(fits.ext+ext,card);
	 sprintf(card,"CD1_2   = %20.10e /                                                ",m[0][1]);
	 addcard(fits.ext+ext,card);
	 sprintf(card,"CD2_1   = %20.10e /                                                ",m[1][0]);
	 addcard(fits.ext+ext,card);
	 sprintf(card,"CD2_2   = %20.10e /                                                ",m[1][1]);
	 addcard(fits.ext+ext,card);
      }
   }
   freeimg(fits.img.data,fits.img.X,fits.img.Y,fits.img.Z);
   fits.img.X=fits.img.Y=fits.img.Z=0;
   fits.img.bits=16;
   fits.img.data=allocimg(0,0,0);
   return;
}

void WFPC2mask(int tp) {
   int x,y,z,idq;
   for (z=0;z<fits.Next;z++) for (y=0;y<fits.ext[z].Y;y++) for (x=0;x<fits.ext[z].X;x++) {
      if (tp==1) idq = dq.img.data[z][y][x];
      else idq = dq.ext[z].data[0][y][x];
      if (idq&8 || fits.ext[z].data[0][y][x]>3500) fits.ext[z].data[0][y][x]=safeup(DMAX);
      else if (idq&5047 || y==fits.ext[z].Y-1 || x==fits.ext[z].X-1) fits.ext[z].data[0][y][x]=safedown(DMIN); // flags 1, 2, 4, 16, 32, 128, 256, 512, 4096=CR
      if (idq&2112) printf("Illegal data quality value %d %d\n",idq&64,idq&2048); // flags 64, 2048(no 1024 check??)
   }
   return;
}

void WFPC2pixcorr(void) {
   int x,y,z;
   float DMIN1,DMAX1;

   if (DMAX>0.) DMAX1=1.5*DMAX;
   else DMAX1=0.;
   if (DMIN<0.) DMIN1=1.5*DMIN;
   else DMIN1=0.;
   for (z=0;z<fits.Next;z++) for (y=0;y<fits.ext[z].Y;y++) for (x=0;x<fits.ext[z].X;x++) if (fits.ext[z].data[0][y][x]>DMIN && fits.ext[z].data[0][y][x]<DMAX) {
      fits.ext[z].data[0][y][x]*=WFPC2xsize(z,x,y)*WFPC2ysize(z,x,y);
   }
   else if (fits.ext[z].data[0][y][x]>=DMAX) fits.ext[z].data[0][y][x]=safeup(DMAX1);
   else fits.ext[z].data[0][y][x]=safedown(DMIN1);
   DMIN=DMIN1;
   DMAX=DMAX1;
   return;
}

void WFPC2drzpixcorr(double mult) {
   int x,y;

   for (y=0;y<fits.ext[0].Y;y++) for (x=0;x<fits.ext[0].X;x++) {
      if (fits.ext[1].data[0][y][x]==0.) fits.ext[0].data[0][y][x]=safedown(DMIN); // no weight
   }
   // mult*=GAIN;  // assuming drizzled image is in cps
   DMAX*=EXP/mult;
   DMIN*=EXP/mult;
   for (y=0;y<fits.ext[0].Y;y++) for (x=0;x<fits.ext[0].X;x++) {
      fits.ext[0].data[0][y][x]*=EXP/mult;
   }
   for (y=1;y<fits.Next;y++) freeim(fits.ext+y);
   fits.Next=1;
   return;
}

void WFPC2setcards(int tp) {
   int ext;
   if (tp==-1) {
      insertcards(fits.ext,GAIN,RN,EXP,DMIN,DMAX,EPOCH,0.0,EXP0);
      return;
   }
   for (ext=0;ext<4;ext++) {
      insertcards(fits.ext+ext,GAIN,RN,EXP,DMIN,DMAX,EPOCH,0.0,EXP0);
   }
   return;
}

int main(int argc,char**argv) {
   int i,tp;
   //int fcorr=0;

   if (argc<2) {
      printf("Usage: %s <<-flags>> <<fits file, data quality file> for each file>\n",*argv);
      printf(" -keepcr      leaves fixed cosmic ray pixels in image\n");
      //printf(" -flashdir=X  sets directory prefix for flash correction files\n");
      return 1;
   }
   for (i=1;i<argc;i++) if (!strcasecmp(argv[i],"-KEEPCR")) MASKCR=0;
   else if (!strncasecmp(argv[i],"-EXPTIME=",9)) FIXED_ET=atof(argv[i]+9);
   //else if (!strncasecmp(argv[i],"-FLASHDIR=",10)) strcpy(FLASHDIR,argv[i]+10);
   else {
      readfits(argv[i],&fits,1);
      tp=WFPC2type(&fits);
      if (strcmp(getcardval(&fits.img,"DOLWFPC2",0),"")) printf("%s already run through wfpc2mask\n",argv[i]);
      else if (tp==0) printf("%s is not a WFPC2 fits file\n",argv[i]);
      else {
	 if (tp>0) {
	    if (i==argc-1) {
	       printf("No data quality file listed for %s\n",argv[i]);
	       return 1;
	    }
	    readfits(argv[i+1],&dq,1);
	    if (tp!=WFPC2type(&dq)) printf("%s and %s have different image types\n",argv[i],argv[i+1]);
	 }
	 WFPC2getcards(0);
	 if (tp==1) {
	    WFPC2toextensions();
	 }
	 if (tp>0) {
	    WFPC2mask(tp);
	    WFPC2pixcorr();
	 }
	 else WFPC2drzpixcorr(1.0);
	 WFPC2setcards(tp);
	 if (tp>0) {
	    addcard(fits.ext,"DOLWFPC2=                    0 / DOLPHOT WFPC2 tag                              ");
	    addcard(fits.ext+1,"DOLWFPC2=                    1 / DOLPHOT WFPC2 tag                              ");
	    addcard(fits.ext+2,"DOLWFPC2=                    2 / DOLPHOT WFPC2 tag                              ");
	    addcard(fits.ext+3,"DOLWFPC2=                    3 / DOLPHOT WFPC2 tag                              ");
	 }
	 else {
	    addcard(fits.ext,"DOLWFPC2=                   -1 / DOLPHOT WFPC2 tag                              ");
	 }
	 writefits(argv[i],&fits,1);
      }
      //if (fcorr) freefits(&flash);
      freefits(&fits);
      if (tp>0) {
	 freefits(&dq);
	 i++;  // skip DQ image
      }
   }
   return 0;
}
