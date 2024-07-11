#include "fits.h"
#include "acspsfdata.h"

ftype fits,fits2,pam;

void ACSunmask(int ext,char*pamfn) {
   int x,y,cm;
   double GAIN,RN,EXP,EPOCH,AIR,EXP0;
   float DMIN,DMAX,DMIN1,DMAX1;
   char str[321],*ptr;

   // parse cards
   strcpy(str,getcardval(fits.ext+ext,"DOL_ACS",0));
   cm=strtol(str,&ptr,10);
   parsecards(&fits.img,&GAIN,&RN,&EXP,&DMIN,&DMAX,&EPOCH,&AIR,&EXP0,0,1);
   parsecards(fits.ext+ext,&GAIN,&RN,&EXP,&DMIN,&DMAX,&EPOCH,&AIR,&EXP0,1,0);

   // undo pixel corrections
   if (DMAX>0.) DMAX1=1.5*DMAX;
   else DMAX1=0.;
   if (DMIN<0.) DMIN1=1.5*DMIN;
   else DMIN1=0.;
   sprintf(str,"%s/acs/data/%s","/Users/patito/Documents/topicos/topicos/ModuloII",pamfn);
   readfits(str,&pam,0);
   if (pam.img.X!=fits.ext[ext].X || pam.img.Y!=fits.ext[ext].Y || pam.img.Z!=fits.ext[ext].Z) {
      printf("Error in PAM file %s\n",pamfn);
      return;
   }
   for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) if (fits.ext[ext].data[0][y][x]>DMIN && fits.ext[ext].data[0][y][x]<DMAX) {
      fits.ext[ext+1].data[0][y][x]=RN;
      if (fits.ext[ext].data[0][y][x]>0) fits.ext[ext+1].data[0][y][x]+=fits.ext[ext].data[0][y][x];
      fits.ext[ext].data[0][y][x]/=pam.img.data[0][y][x]/acs_ctmult[cm];
      fits.ext[ext+1].data[0][y][x]=sqrt(fits.ext[ext+1].data[0][y][x])/(pam.img.data[0][y][x]/acs_ctmult[cm]);
      fits.ext[ext+2].data[0][y][x]=0;
   }
   else if (fits.ext[ext].data[0][y][x]>=DMAX) {
      fits.ext[ext].data[0][y][x]=safeup(DMAX1);
      fits.ext[ext+1].data[0][y][x]=0.0;
      fits.ext[ext+2].data[0][y][x]=2304;
   }
   else {
      fits.ext[ext].data[0][y][x]=safedown(DMIN1);
      fits.ext[ext+1].data[0][y][x]=0.0;
      fits.ext[ext+2].data[0][y][x]=5823;
   }
   freefits(&pam);
   DMIN=DMIN1;
   DMAX=DMAX1;
}

void clearcards(imtype *img) {
   int i,j;
   for (i=0;i<img->Ncards;i++) if (!strncasecmp(img->cards[i],"DOL_ACS =",9)) {
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
      // HRC with splitgroups
      if (fits.Next==0 && !strcmp(getcardval(&fits.img,"DOL_ACS",0),"0")) {
	 fits.Next=3;
	 fits.ext=(imtype*)realloc(fits.ext,sizeof(imtype)*3);
	 if (!fits.ext) {fprintf(stderr,"Memory allocation error\n"); return -1;}
	 imcopy(&fits.img,fits.ext);
	 imcopy(&fits.img,fits.ext+1);
	 imcopy(&fits.img,fits.ext+2); fits.ext[2].bits=16;
	 freeimg(fits.img.data,fits.img.X,fits.img.Y,fits.img.Z);
	 fits.img.X=fits.img.Y=fits.img.Z=0;
	 fits.img.bits=16;
	 fits.img.data=allocimg(0,0,0);
	 ACSunmask(0,"hrc_pam.fits");
      }
      // HRC without splitgroups
      else if (fits.Next==1 && !strcmp(getcardval(fits.ext,"DOL_ACS",0),"0")) {
	 fits.Next=3;
	 fits.ext=(imtype*)realloc(fits.ext,sizeof(imtype)*3);
	 if (!fits.ext) {fprintf(stderr,"Memory allocation error\n"); return -1;}
	 imcopy(fits.ext,fits.ext+1);
	 imcopy(fits.ext,fits.ext+2); fits.ext[2].bits=16;
	 ACSunmask(0,"hrc_pam.fits");
      }
      // WFC with splitgroups
      else if (fits.Next==0 && !strcmp(getcardval(&fits.img,"DOL_ACS",0),"2")) {
	 if (i+1>=argc) {
	    fprintf(stderr,"WFC chip1 provided without chip2\n");
	    return -1;
	 }
	 readfits(argv[i+1],&fits2,1);
	 if (fits2.Next!=0 || strcmp(getcardval(&fits2.img,"DOL_ACS",0),"1")) {
	    fprintf(stderr,"WFC chip1 provided without chip2\n");
	    return -1;
	 }
	 fits.Next=6;
	 fits.ext=(imtype*)realloc(fits.ext,sizeof(imtype)*6);
	 if (!fits.ext) {fprintf(stderr,"Memory allocation error\n"); return -1;}
	 imcopy(&fits.img,fits.ext);
	 imcopy(&fits.img,fits.ext+1);
	 imcopy(&fits.img,fits.ext+2); fits.ext[2].bits=16;
	 imcopy(&fits2.img,fits.ext+3);
	 imcopy(&fits2.img,fits.ext+4);
	 imcopy(&fits2.img,fits.ext+5); fits.ext[5].bits=16;
	 freefits(&fits2);
	 freeimg(fits.img.data,fits.img.X,fits.img.Y,fits.img.Z);
	 fits.img.X=fits.img.Y=fits.img.Z=0;
	 fits.img.bits=16;
	 fits.img.data=allocimg(0,0,0);
	 ACSunmask(0,"wfc2_pam.fits");
	 ACSunmask(3,"wfc1_pam.fits");
	 i++;
      }
      // WFC without splitgroups
      else if (fits.Next==2 && !strcmp(getcardval(fits.ext,"DOL_ACS",0),"2") && !strcmp(getcardval(fits.ext+1,"DOL_ACS",0),"1")) {
	 fits.Next=6;
	 fits.ext=(imtype*)realloc(fits.ext,sizeof(imtype)*6);
	 if (!fits.ext) {fprintf(stderr,"Memory allocation error\n"); return -1;}
	 imcopy(fits.ext+1,fits.ext+3);
	 imcopy(fits.ext+1,fits.ext+4);
	 imcopy(fits.ext+1,fits.ext+5); fits.ext[5].bits=16;
	 imcopy(fits.ext,fits.ext+1);
	 imcopy(fits.ext,fits.ext+2); fits.ext[2].bits=16;
	 ACSunmask(0,"wfc2_pam.fits");
	 ACSunmask(3,"wfc1_pam.fits");
      }
      else {
	 fprintf(stderr,"%s is not an ACS file run through acsmask\n",argv[i]);
	 return -1;
      }
      clearcards(&fits.img);
      for (j=0;j<fits.Next;j++) clearcards(fits.ext+j);
      writefits(argv[i0],&fits,1);
      freefits(&fits);
   }
   return 0;
}
