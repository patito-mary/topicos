#include "fits.h"
#include "acspsfdata.h"


ftype fits,pam;
//ftype flash;
double RN,EXP,EXP0,EPOCH,FIXED_ET=-1.;
float DMIN,DMAX;
int MASKCR=1,MASKSAT=1,NCOMBINE=1,FIXED_NC=-1,USE_WHT=0;
//char FLASHDIR[161]="";
int offsetx,offsety;
int Next0;
int FORCE_SUBARRAY=0,FORCE_OFFSETX=0,FORCE_OFFSETY=0; // 3=HRC, 5=WFC1, 6=WFC2

char* eitherstring(char*val,int err) {
   // astrodrizzle
   if (Next0>=4 && fits.ext[3].tfields>0) {
      char *rval = gettablevalstring(fits.ext+3,val,0,err);
      if (*rval!=0) return rval;
   }
   else if (Next0>=3 && fits.ext[2].tfields>0) {
      char *rval = gettablevalstring(fits.ext+2,val,0,err);
      if (*rval!=0) return rval;
   }
   // FLT or multidrizzle
   return getcardval(&(fits.img),val,err);
}

double eitherdouble(char*val,int err) {
   // astrodrizzle
   if (Next0==4 && fits.ext[3].tfields>0) {
      double rval = gettablevaldouble(fits.ext+3,val,0,0,err);
      if (rval!=0) return rval;
   }
   else if (Next0>=3 && fits.ext[2].tfields>0) {
      double rval = gettablevaldouble(fits.ext+2,val,0,0,err);
      if (rval!=0) return rval;
   }
   // FLT or multidrizzle
   return atof(getcardval(&(fits.img),val,err));
}

int ACStype(ftype *f) {
   int i;
   char detector[81],aperture[81];

   if (strcmp(eitherstring("FILETYPE",1),"SCI") || strcmp(getcardval(&(fits.img),"TELESCOP",1),"HST") || strcmp(eitherstring("INSTRUME",1),"ACS")) {printf("**Format error (filetype,telescop,instrume)\n"); exit(-1);}

   strcpy(detector,eitherstring("DETECTOR",1));
   strcpy(aperture,eitherstring("APERTURE",1));

   if (FORCE_SUBARRAY!=0) {
      int XMAX=4096,YMAX=2048;
      if (FORCE_SUBARRAY>3) {
	 if (strcmp(detector,"WFC")) printf("Warning: detector is %s; WFC was commanded\n",detector);
      }
      else {
	 XMAX=YMAX=1024;
	 if (strcmp(detector,"HRC")) printf("Warning: detector is %s; HRC was commanded\n",detector);
      }
      if (Next0<3) {
	 printf("Error: at least three extensions required\n");
	 exit(-1);
      }
      else if (f->ext[0].Z!=1) {
	 printf("Error: primary image has third dimension\n");
	 exit(-1);
      }
      else if (FORCE_OFFSETX<0 || FORCE_OFFSETY<0 || f->ext[0].X+FORCE_OFFSETX>XMAX || f->ext[0].Y+FORCE_OFFSETY>YMAX) {
	 printf("Error: subarray would fall outside %dx%d image\n",XMAX,YMAX);
	 exit(-1);
      }
      else if (f->ext[0].X!=f->ext[1].X || f->ext[0].Y!=f->ext[1].Y || f->ext[1].Z!=1 || f->ext[0].X!=f->ext[2].X || f->ext[0].Y!=f->ext[2].Y || f->ext[2].Z!=1) {
	 printf("Error: data quality image does not match primary image size\n");
	 exit(-1);
      }
      else if (FORCE_SUBARRAY>3 && Next0!=3 && Next0<7) printf("Warning: expect three-extension image for subarray\n");
      else if (FORCE_SUBARRAY<=3 && Next0!=3 && Next0<6) printf("Warning: expect three-extension image for subarray\n");

      offsetx = FORCE_OFFSETX;
      offsety = FORCE_OFFSETY;
      return FORCE_SUBARRAY;
   }

   offsetx = 0;
   offsety = 0;

   // WFC image types
   if (!strcmp(detector,"WFC")) {
      // undrizzled WFC
      if (Next0==6 || Next0>=13) {
	 for (i=0;i<6 && f->ext[i].X==4096 && f->ext[i].Y==2048 && f->ext[i].Z==1;i++);
	 if (i==6) return 1;
	 return 0;
      }
      if (Next0==3 || Next0>=7) {
	 if (!strcmp(aperture,"WFC1-512")) {
            for (i=0;i<3 && f->ext[i].X==512 && f->ext[i].Y==512 && f->ext[i].Z==1;i++);
            if (i==3) {offsetx=3584; offsety=1536; return 5;}
	 }
	 if (!strcmp(aperture,"WFC1-1K")) {
            for (i=0;i<3 && f->ext[i].X==1024 && f->ext[i].Y==1024 && f->ext[i].Z==1;i++);
            if (i==3) {offsetx=3072; offsety=1024; return 5;}
	 }
	 if (!strcmp(aperture,"WFC1-2K")) {
            for (i=0;i<3 && f->ext[i].X==2048 && f->ext[i].Y==2046 && f->ext[i].Z==1;i++);
            if (i==3) {offsetx=2048; offsety=1; return 5;}
	 }
	 if (!strcmp(aperture,"WFC2-2K")) {
            for (i=0;i<3 && f->ext[i].X==2048 && f->ext[i].Y==2046 && f->ext[i].Z==1;i++);
            if (i==3) {offsetx=0; offsety=1; return 6;}
	 }
      }
      // drizzled WFC
      if (Next0==3 || Next0==4 || Next0==5) {
         if (f->ext[0].Z==1 && f->ext[1].X==f->ext[0].X && f->ext[1].Y==f->ext[0].Y && f->ext[1].Z==f->ext[0].Z && (USE_WHT!=0 || (f->ext[2].X==f->ext[0].X && f->ext[2].Y==f->ext[0].Y))) {
	    printf("Irregular size; assuming drizzled\n");
            return 2;
         }
      }
   }

   // HRC image types
   if (!strcmp(detector,"HRC")) {
      //undrizzled HRC
      if (Next0==3 || Next0>=6) {
	 for (i=0;i<3 && f->ext[i].X==1024 && f->ext[i].Y==1024 && f->ext[i].Z==1;i++);
	 if (i==3) return 3;
	 if (!strcmp(aperture,"HRC-512")) {
            for (i=0;i<3 && f->ext[i].X==513 && f->ext[i].Y==512 && f->ext[i].Z==1;i++);
            if (i==3) {offsetx=0; offsety=0; return 3;}
	 }
      }
      //drizzled
      if (Next0==3 || Next0==4) {
	 if (f->ext[0].Z==1 && f->ext[1].X==f->ext[0].X && f->ext[1].Y==f->ext[0].Y && f->ext[1].Z==f->ext[0].Z && (USE_WHT!=0 || (f->ext[2].X==f->ext[0].X && f->ext[2].Y==f->ext[0].Y))) {
	    printf("Irregular size; assuming drizzled\n");
	    return 4;
	 }
      }
   }
   return 0;
}

void ACSgetcards(int ext,int MODE,int drz) {
   int x,y;
   char str[81];

   if (!strcmp(eitherstring("CCDAMP",1),"ABCD")) {
      if (MODE==1) {
	 RN=0.5*(eitherdouble("READNSEC",1)+eitherdouble("READNSED",1));
      }
      else if (MODE==2) {
	 RN=0.5*(eitherdouble("READNSEA",1)+eitherdouble("READNSEB",1));
      }
      else RN=0.25*(eitherdouble("READNSEA",1)+eitherdouble("READNSEB",1)+eitherdouble("READNSEC",1)+eitherdouble("READNSED",1));
   }
   else if (!strcmp(eitherstring("CCDAMP",1),"A")) RN=eitherdouble("READNSEA",1);
   else if (!strcmp(eitherstring("CCDAMP",1),"B")) RN=eitherdouble("READNSEB",1);
   else if (!strcmp(eitherstring("CCDAMP",1),"C")) RN=eitherdouble("READNSEC",1);
   else if (!strcmp(eitherstring("CCDAMP",1),"D")) RN=eitherdouble("READNSED",1);
   else {printf("**Format error (ccdamp)\n"); exit(-1);}

   if (FIXED_ET>0) EXP=FIXED_ET;
   else EXP=atof(getcardval(&(fits.img),"EXPTIME",1));
   if (drz) {
      DMIN=DMAX=fits.ext[ext].data[0][0][0];
      for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) {
	 if (fits.ext[ext].data[0][y][x]<DMIN) DMIN=fits.ext[ext].data[0][y][x];
	 else if (fits.ext[ext].data[0][y][x]>DMAX) DMAX=fits.ext[ext].data[0][y][x];
      }
   }
   else {
      DMIN=safedown(atof(getcardval(fits.ext+ext,"GOODMIN",1)));
      DMAX=safeup(atof(getcardval(fits.ext+ext,"GOODMAX",1)));
   }
   EPOCH=0.5*(atof(getcardval(&(fits.img),"EXPSTART",1))+atof(getcardval(&(fits.img),"EXPEND",1)));
   if (FIXED_NC>0) {
      NCOMBINE=FIXED_NC;
      printf("Setting number of images to %d\n",NCOMBINE);
   }
   else {
      strcpy(str,getcardval(&(fits.img),"NRPTEXP",0));
      if (*str==0) {
	 printf("Error: NRPTEXP not found; assuming 1\n");
	 NCOMBINE=1;
      }
      else NCOMBINE=atoi(str);
      strcpy(str,getcardval(&(fits.img),"CRSPLIT",0));
      if (*str==0) printf("Error: CRSPLIT not found; assuming 1\n");
      else NCOMBINE*=atoi(str);
   }
   EXP0=EXP/NCOMBINE;
   return;
}

//ignoring type 64 (warm pixel) on assumption it is fixed OK
void ACSmask(ftype*f,int ext,int drz) {
   int x,y,dq,i;
   if (!drz) {
      float DMAX1=DMAX, DMIN1=DMIN;
      if (DMAX>0.) DMAX1*=2.;
      else DMAX1*=0.5;
      if (DMIN<0.) DMIN1*=2.;
      else DMIN1*=0.5;
      for (y=0;y<f->ext[ext].Y;y++) for (x=0;x<f->ext[ext].X;x++) {
	 dq=(int)(f->ext[ext+2].data[0][y][x]+0.5);
	 if (dq&2048) f->ext[ext].data[0][y][x]=safeup(DMAX1); // A/D saturated
	 else if (dq&256 && MASKSAT) f->ext[ext].data[0][y][x]=safeup(DMAX1); // full well, not A/D
	 else if (dq&1727) f->ext[ext].data[0][y][x]=safedown(DMIN1); // others
	 else if ((dq&12288) && MASKCR) f->ext[ext].data[0][y][x]=safedown(DMIN1); // 4096=found in drizzle; 8192=found in crsplit
	 else if (f->ext[ext].data[0][y][x]>=DMAX && MASKSAT) f->ext[ext].data[0][y][x]=safeup(DMAX1); // exceed GOODMAX
	 else if (f->ext[ext].data[0][y][x]<=DMIN) f->ext[ext].data[0][y][x]=safedown(DMIN1); // below GOODMIN
      }
      DMIN=DMIN1;
      DMAX=DMAX1;
   }
   else {
      for (y=0;y<f->ext[ext].Y;y++) for (x=0;x<f->ext[ext].X;x++) {
	 dq=0;
	 if (USE_WHT) {
	    if (f->ext[ext+1].data[0][y][x]>0.0) dq=1;
	 }
	 else {
	    for (i=0;i<f->ext[ext+2].Z && !dq;i++) if ((int)(f->ext[ext+2].data[i][y][x]+0.5)!=0) dq=1;
	 }
	 if (dq==0) f->ext[ext].data[0][y][x]=safedown(DMIN);
      }
   }
   return;
}

void ACSpixcorr(int ext,char*pamfn,double mult) {
   int x,y;
   float DMIN1,DMAX1;
   char fn[321];

   if (DMAX>0.) DMAX1=1.5*DMAX;
   else DMAX1=0.;
   if (DMIN<0.) DMIN1=1.5*DMIN;
   else DMIN1=0.;
   sprintf(fn,"%s/acs/data/%s","/Users/patito/Documents/topicos/topicos/ModuloII",pamfn);
   readfits(fn,&pam,0);
   if (pam.img.X<fits.ext[ext].X+offsetx || pam.img.Y<fits.ext[ext].Y+offsety || pam.img.Z!=fits.ext[ext].Z) {
      printf("Error in PAM file %s\n",pamfn);
      return;
   }
   for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) if (fits.ext[ext].data[0][y][x]>DMIN && (fits.ext[ext].data[0][y][x]<DMAX || !MASKSAT)) {
      fits.ext[ext].data[0][y][x]*=pam.img.data[0][y+offsety][x+offsetx]/mult;
      fits.ext[ext+1].data[0][y][x]*=pam.img.data[0][y+offsety][x+offsetx]/mult;
   }
   else if (fits.ext[ext].data[0][y][x]>=DMAX) fits.ext[ext].data[0][y][x]=safeup(DMAX1);
   else fits.ext[ext].data[0][y][x]=safedown(DMIN1);
   freefits(&pam);
   DMIN=DMIN1;
   DMAX=DMAX1;
   return;
}

void ACSdrzpixcorr(int ext,double mult) {
   int x,y;

   DMAX*=EXP/mult;
   DMIN*=EXP/mult;
   for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) {
      fits.ext[ext].data[0][y][x]*=EXP/mult;
   }
   return;
}

/*
void ACStrim(int ext,int xf,int yf) {
   int dx=0,dy=0,y,xsize;
   imgtype tmp;
   char str[81],*ptr;

   strcpy(str,getcardval(flash.ext+ext,"LTV1",0));
   if (str[0]) {
      dx=strtol(str,&ptr,10);
      if (ptr==str || (*ptr && strcmp(ptr,".0")) || dx<0 || dx+xf>flash.ext[ext].X) {printf("**Illegal X shift %s\n",str); exit(-1);}
   }
   strcpy(str,getcardval(flash.ext+ext,"LTV2",0));
   if (str[0]) {
      dy=strtol(str,&ptr,10);
      if (ptr==str || (*ptr && strcmp(ptr,".0")) || dy<0 || dy+yf>flash.ext[ext].Y) {printf("**Illegal Y shift %s\n",str); exit(-1);}
   }
   tmp=allocimg(xf,yf,1);
   xsize=sizeof(float)*xf;
   for (y=0;y<yf;y++) memcpy(tmp[0][y],flash.ext[ext].data[0][y+dy]+dx,xsize);
   freeimg(flash.ext[ext].data,flash.ext[ext].X,flash.ext[ext].Y,flash.ext[ext].Z);
   flash.ext[ext].data=tmp;
   flash.ext[ext].X=xf;
   flash.ext[ext].Y=yf;
   return;
}

int loadflash(void) {
   char str[161],*ptr,fcur[81],fn[321];
   int ext,x,y,dq;
   double tflash;

   strcpy(str,getcardval(&(fits.img),"FLASHDUR",1));
   if (!strcmp(str,"0.0")) return 0;
   tflash=strtod(str,&ptr);
   if (ptr==str || *ptr) {printf("**Format error: FLASHDUR\n"); exit(-1);}
   strcpy(str,getcardval(&(fits.img),"FLSHCORR",1));
   if (!strcmp(str,"OMIT")) return 0;
   if (strcmp(str,"COMPLETE")) {printf("**Format error: FLSHCORR\n"); exit(-1);}
   strcpy(fcur,getcardval(&(fits.img),"FLASHCUR",1));
   if (strcmp(fcur,"LOW") && strcmp(fcur,"MED") && strcmp(fcur,"HIGH")) {printf("**Format error: FLASHCUR\n"); exit(-1);}
   strcpy(str,getcardval(&(fits.img),"FLASHSTA",1));
   if (strcmp(str,"SUCCESSFUL")) {printf("**Format error: FLASHSTA\n"); exit(-1);}
   strcpy(str,getcardval(&(fits.img),"FLSHFILE",1));
   if (strncmp(str,"jref$",5)) {printf("**Format error: FLSHFILE\n"); exit(-1);}
   sprintf(fn,"%s%s",FLASHDIR,str+5);
   readfits(fn,&flash,0);
   if (strcmp(getcardval(&(flash.img),"FLASHCUR",1),fcur)) {printf("**Flash/Data mismatch: FLASHCUR\n"); exit(-1);}
   strcpy(str,getcardval(&(flash.img),"FLSHCORR",1));
   if (strcmp(str,"OMIT")) {printf("**Flash format error: FLSHCORR\n"); exit(-1);}
   strcpy(str,getcardval(&(flash.img),"FLASHSTA",1));
   if (strcmp(str,"SUCCESSFUL")) {printf("**Flash format error FLASHSTA\n"); exit(-1);}
   if (fits.Next!=flash.Next) {printf("**Flash format error: Next\n"); exit(-1);}
   for (ext=0;ext<fits.Next;ext+=3) {
      for (y=0;y<flash.ext[ext].Y;y++) for (x=0;x<flash.ext[ext].X;x++) {
	 dq=(int)(flash.ext[ext+2].data[0][y][x]+0.5);
	 if (dq) {printf("**Flash DQ error: nonzero value (%d)\n",dq); exit(-1);}
	 flash.ext[ext].data[0][y][x]*=tflash;
      }
      if (flash.ext[ext].X==4144 && flash.ext[ext].Y==2068) ACStrim(ext,4096,2048);
      else if (flash.ext[ext].X==1062 && flash.ext[ext].Y==1044) ACStrim(ext,1024,1024);
      else {printf("**Flash size error: %dx%d\n",flash.ext[ext].X,flash.ext[ext].Y); exit(-1);}
      if (ext) memcpy(flash.ext+ext/3,flash.ext+ext,sizeof(imtype));
      freeim(flash.ext+ext+1);
      freeim(flash.ext+ext+2);
   }
   flash.Next/=3;
   return 1;
}

void ACSflashcorr(int ext) {
   double DMIN1=0.,DMAX1=0.;
   int x,y;

   if (flash.ext[ext].X!=fits.ext[ext].X || flash.ext[ext].Y!=fits.ext[ext].Y || flash.ext[ext].Z!=fits.ext[ext].Z) {printf("**Error in FLASH file (size mismatch)\n"); exit(-1);}
   for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) {
      if (flash.ext[ext].data[0][y][x]>DMAX1) DMAX1=flash.ext[ext].data[0][y][x];
      else if (flash.ext[ext].data[0][y][x]<DMIN1) DMIN1=flash.ext[ext].data[0][y][x];
   }
   DMAX1+=DMAX+1.;
   DMIN1+=DMIN-1.;
   for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) if (fits.ext[ext].data[0][y][x]>DMIN && fits.ext[ext].data[0][y][x]<DMAX) fits.ext[ext].data[0][y][x]+=flash.ext[ext].data[0][y][x];
   else if (fits.ext[ext].data[0][y][x]>=DMAX) fits.ext[ext].data[0][y][x]=safeup(DMAX1);
   else fits.ext[ext].data[0][y][x]=safedown(DMIN1);
   DMAX=DMAX1;
   DMIN=DMIN1;
   return;
}
*/

void ACSsetcards(int ext) {
   int x,y;
   chiptype data,noise;
   double I=0,X=0,Y=0,XX=0,XY=0,tGN,tRN;

   data=fits.ext[ext].data[0];
   noise=fits.ext[ext+1].data[0];
   for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) if (data[y][x]>DMIN && data[y][x]<DMAX) {
      I++;
      Y+=noise[y][x]*noise[y][x];
      if (data[y][x]>0.) {
	 X+=data[y][x];
	 XX+=data[y][x]*data[y][x];
	 XY+=data[y][x]*noise[y][x]*noise[y][x];
      }
   }
   tGN=(XX*I-X*X)/(XY*I-X*Y);
   tRN=sqrt((XX*Y-X*XY)/(XX*I-X*X))*tGN;
   printf("%f %f\n",tGN,tRN);
   insertcards(fits.ext+ext,1.,RN*sqrt(NCOMBINE),EXP,DMIN,DMAX,EPOCH,0.0,EXP0);
   freeim(fits.ext+ext+1);
   freeim(fits.ext+ext+2);
   fits.Next-=2;
   for (x=ext+1;x<fits.Next;x++) memcpy(fits.ext+x,fits.ext+x+2,sizeof(imtype));
   return;
}

int main(int argc,char**argv) {
   int i,ext,tp;
   //int fcorr=0;
   char card[81];
   char *ptr;

   if (argc<2) {
      printf("Usage: %s <<-flags>> <fits files>\n",*argv);
      printf(" -keepcr      leaves fixed cosmic ray pixels in image\n");
      printf(" -keepsat     leaves saturated pixels in image\n");
      printf(" -exptime=#   overrides exposure time; needed for _drz images\n");
      printf(" -ncombine=#  overrides NCOMBINE keyword\n");
      printf(" -usewht      to use weight extension instead of context for _drz\n");
      printf(" -wfc1sub=#,# to force WFC1 subarray shifted by #,# from full image\n");
      printf(" -wfc2sub=#,# to force WFC2 subarray shifted by #,# from full image\n");
      printf(" -hrcsub=#,#  to force HRC subarray shifted by #,# from full image\n");
      //printf(" -flashdir=X  sets directory prefix for flash correction files\n");
      return 1;
   }
   for (i=1;i<argc;i++) if (!strcasecmp(argv[i],"-KEEPCR")) MASKCR=0;
   else if (!strcasecmp(argv[i],"-KEEPSAT")) MASKSAT=0;
   else if (!strncasecmp(argv[i],"-EXPTIME=",9)) FIXED_ET=atof(argv[i]+9);
   else if (!strncasecmp(argv[i],"-NCOMBINE=",10)) FIXED_NC=atof(argv[i]+10);
   else if (!strcasecmp(argv[i],"-USEWHT")) USE_WHT=1;
   else if (!strncasecmp(argv[i],"-hrcsub=",8)) {
      FORCE_SUBARRAY = 3;
      FORCE_OFFSETX = strtol(argv[i]+8,&ptr,10);
      if (*ptr!=',') {
	 printf("Cannot parse subarray definition \"%s\"\n",argv[i]);
	 return 1;
      }
      FORCE_OFFSETY = atoi(ptr+1);
   }
   else if (!strncasecmp(argv[i],"-wfc1sub=",9)) {
      FORCE_SUBARRAY = 5;
      FORCE_OFFSETX = strtol(argv[i]+9,&ptr,10);
      if (*ptr!=',') {
	 printf("Cannot parse subarray definition \"%s\"\n",argv[i]);
	 return 1;
      }
      FORCE_OFFSETY = atoi(ptr+1);
   }
   else if (!strncasecmp(argv[i],"-wfc2sub=",9)) {
      FORCE_SUBARRAY = 6;
      FORCE_OFFSETX = strtol(argv[i]+9,&ptr,10);
      if (*ptr!=',') {
	 printf("Cannot parse subarray definition \"%s\"\n",argv[i]);
	 return 1;
      }
      FORCE_OFFSETY = atoi(ptr+1);
   }
   //else if (!strncasecmp(argv[i],"-FLASHDIR=",10)) strcpy(FLASHDIR,argv[i]+10);
   else {
      readfits(argv[i],&fits,1);
      Next0=fits.Next;
      tp=ACStype(&fits);
      if (Next0>0 && strcmp(getcardval(fits.ext,"DOL_ACS",0),"")) printf("%s already run through acsmask\n",argv[i]);
      else if (tp<1 || tp>6) printf("%s is not an ACS fits file\n",argv[i]);
      else {
	 //fcorr=loadflash();
	 if (tp==1) { // FLT or CRJ WFC image
	    ACSgetcards(0,2,0);
	    ACSmask(&fits,0,0);
	    ACSpixcorr(0,"wfc2_pam.fits",acs_ctmult[2]);
	    //if (fcorr) ACSflashcorr(0);
	    ACSsetcards(0);
	    ACSgetcards(1,1,0);
	    ACSmask(&fits,1,0);
	    ACSpixcorr(1,"wfc1_pam.fits",acs_ctmult[1]);
	    //if (fcorr) ACSflashcorr(1);
	    ACSsetcards(1);
	    addcard(fits.ext,"DOL_ACS =                    2 / DOLPHOT ACS tag                                ");
	    addcard(fits.ext+1,"DOL_ACS =                    1 / DOLPHOT ACS tag                                ");
	    for (ext=2;ext<fits.Next;ext++) freeim(fits.ext+ext); // new pipeline output
	    fits.Next=2;
	 }
	 else if (tp==2) { // DRZ WFC image
	    //if (fcorr) printf("%s is drizzled and flashed\n",argv[i]);
	    //else {
	    ACSgetcards(0,-1,1);
	    ACSmask(&fits,0,1);
	    ACSdrzpixcorr(0,acs_ctmult[1]);
	    ACSsetcards(0);
	    for (ext=1;ext<fits.Next;ext++) freeim(fits.ext+ext);
	    fits.Next=1;
	    addcard(fits.ext,"DOL_ACS =                   -2 / DOLPHOT ACS tag                                ");
	    //}
	 }
	 else if (tp==3) { // FLT or CRJ HRC image
	    ACSgetcards(0,0,0);
	    ACSmask(&fits,0,0);
	    ACSpixcorr(0,"hrc_pam.fits",acs_ctmult[0]);
	    //if (fcorr) ACSflashcorr(0);
	    ACSsetcards(0);
	    for (ext=1;ext<fits.Next;ext++) freeim(fits.ext+ext);
	    fits.Next=1;
	    addcard(fits.ext,"DOL_ACS =                    0 / DOLPHOT ACS tag                                ");
	 }
	 else if (tp==4) { // DRZ HRC image
	    //if (fcorr) printf("%s is drizzled and flashed\n",argv[i]);
	    //else {
	    ACSgetcards(0,-2,1);
	    ACSmask(&fits,0,1);
	    ACSdrzpixcorr(0,acs_ctmult[0]);
	    ACSsetcards(0);
	    for (ext=1;ext<fits.Next;ext++) freeim(fits.ext+ext);
	    fits.Next=1;
	    addcard(fits.ext,"DOL_ACS =                   -1 / DOLPHOT ACS tag                                ");
	    //}
	 }
	 else if (tp==5) { // FLT or CRJ WFC1 only
	    ACSgetcards(0,1,0);
	    ACSmask(&fits,0,0);
	    ACSpixcorr(0,"wfc1_pam.fits",acs_ctmult[1]);
	    //if (fcorr) ACSflashcorr(0);
	    ACSsetcards(0);
	    for (ext=1;ext<fits.Next;ext++) freeim(fits.ext+ext);
	    fits.Next=1;
	    addcard(fits.ext,"DOL_ACS =                    1 / DOLPHOT ACS tag                                ");
	    if (offsetx!=0) {
	       sprintf(card,"DOL_OFFX=                 %4d / Origin of subarray relative to full chip       ",offsetx);
	       addcard(fits.ext,card);
	    }
	    if (offsety!=0) {
	       sprintf(card,"DOL_OFFY=                 %4d / Origin of subarray relative to full chip       ",offsety);
	       addcard(fits.ext,card);
	    }
	 }
	 else if (tp==6) { // FLT or CRJ WFC2 only
	    ACSgetcards(0,2,0);
	    ACSmask(&fits,0,0);
	    ACSpixcorr(0,"wfc2_pam.fits",acs_ctmult[2]);
	    //if (fcorr) ACSflashcorr(0);
	    ACSsetcards(0);
	    for (ext=1;ext<fits.Next;ext++) freeim(fits.ext+ext);
	    fits.Next=1;
	    addcard(fits.ext,"DOL_ACS =                    2 / DOLPHOT ACS tag                                ");
	    if (offsetx!=0) {
	       sprintf(card,"DOL_OFFX=                 %4d / Origin of subarray relative to full chip       ",offsetx);
	       addcard(fits.ext,card);
	    }
	    if (offsety!=0) {
	       sprintf(card,"DOL_OFFY=                 %4d / Origin of subarray relative to full chip       ",offsety);
	       addcard(fits.ext,card);
	    }
	 }
	 writefits(argv[i],&fits,1);
      }
      //if (fcorr) freefits(&flash);
      freefits(&fits);
   }
   return 0;
}
