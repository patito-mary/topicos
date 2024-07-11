#include "../include/fits.h"
#include "wfpc2psfdata.h"
#include "wfpc2distort.h"

ftype fits,fits2;

void WFPC2unmask(void) {
   int x,y,z;
   double GAIN,RN,EXP,EPOCH,AIR,EXP0;
   float DMIN,DMAX,DMIN1,DMAX1;

   fitscopy(&fits,&fits2);
   for (z=0;z<4;z++) {
      fits2.ext[z].bits=16;
      // parse cards
      parsecards(&fits.img,&GAIN,&RN,&EXP,&DMIN,&DMAX,&EPOCH,&AIR,&EXP0,0,1);
      parsecards(fits.ext+z,&GAIN,&RN,&EXP,&DMIN,&DMAX,&EPOCH,&AIR,&EXP0,1,0);

      
      // undo pixel corrections
      if (DMAX>0.) DMAX1=1.5*DMAX;
      else DMAX1=0.;
      if (DMAX1<4000) DMAX1=4000;
      if (DMIN<0.) DMIN1=1.5*DMIN;
      else DMIN1=0.;
      for (y=0;y<fits.ext[z].Y;y++) for (x=0;x<fits.ext[z].X;x++) {
	 if (fits.ext[z].data[0][y][x]>DMIN && fits.ext[z].data[0][y][x]<DMAX) {
	    fits.ext[z].data[0][y][x]/=WFPC2xsize(z,x,y)*WFPC2ysize(z,x,y);
	    fits2.ext[z].data[0][y][x]=0;
	 }
	 else if (fits.ext[z].data[0][y][x]>=DMAX) {
	    fits.ext[z].data[0][y][x]=safeup(DMAX1);
	    fits2.ext[z].data[0][y][x]=8;
	 }
	 else {
	    fits.ext[z].data[0][y][x]=safedown(DMIN1);
	    fits2.ext[z].data[0][y][x]=951;
	 }
      }
   }
}

void clearcards(imtype *img) {
   int i,j;
   for (i=0;i<img->Ncards;i++) if (!strncasecmp(img->cards[i],"DOLWFPC2=",9)) {
      img->Ncards--;
      for (j=i;j<img->Ncards;j++) strcpy(img->cards[j],img->cards[j+1]);
      i--;
   }
}

int main(int argc,char**argv) {
   int i,i0,j;

   if (argc<2) {
      printf("Usage: %s <fits files>\n",*argv);
      return 1;
   }
   for (i=1;i<argc;i++) {
      readfits(argv[i],&fits,1);
      i0=i;
      // With splitgroups
      if (fits.Next==0 && !strcmp(getcardval(&fits.img,"DOLWFPC2",0),"0")) {
	 if (i+3>=argc) {
	    fprintf(stderr,"PC1 provided without full set of input files\n");
	    return -1;
	 }
	 if (i+4>=argc) {
	    fprintf(stderr,"No DQ file specified for %s-%s\n",argv[i],argv[i+3]);
	    return -1;
	 }
	 fits.Next=4;
	 fits.ext=(imtype*)realloc(fits.ext,sizeof(imtype)*4);
	 if (!fits.ext) {fprintf(stderr,"Memory allocation error\n"); return -1;}
	 imcopy(&fits.img,fits.ext);
	 freeimg(fits.img.data,fits.img.X,fits.img.Y,fits.img.Z);
	 fits.img.X=fits.img.Y=fits.img.Z=0;
	 fits.img.bits=16;
	 fits.img.data=allocimg(0,0,0);
	 for (j=1;j<4;j++) {
	    readfits(argv[i+j],&fits2,1);
	    if (fits2.Next!=0 || atoi(getcardval(&fits2.img,"DOLWFPC2",0))!=j) {
	       fprintf(stderr,"%s is not chip WF%d\n",argv[i+j],j+1);
	       return -1;
	    }
	    imcopy(&fits2.img,fits.ext+j);
	    freefits(&fits2);
	 }
	 i+=4;
      }
      // WFC without splitgroups
      else if (fits.Next==4 && !strcmp(getcardval(fits.ext,"DOLWFPC2",0),"0") && !strcmp(getcardval(fits.ext+1,"DOLWFPC2",0),"1") && !strcmp(getcardval(fits.ext+2,"DOLWFPC2",0),"2") && !strcmp(getcardval(fits.ext+3,"DOLWFPC2",0),"3")) {
	 if (i+1>=argc) {
	    fprintf(stderr,"No DQ specified for %s\n",argv[i]);
	    return -1;
	 }
	 i++;
      }
      else {
	 fprintf(stderr,"%s is not an ACS file run through acsmask\n",argv[i]);
	 return -1;
      }
      clearcards(&fits.img);
      for (j=0;j<fits.Next;j++) clearcards(fits.ext+j);
      WFPC2unmask();
      writefits(argv[i0],&fits,1);
      writefits(argv[i],&fits2,1);
      i++;
      freefits(&fits);
      freefits(&fits2);
   }
   return 0;
}
