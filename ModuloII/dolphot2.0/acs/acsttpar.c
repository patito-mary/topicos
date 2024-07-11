#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc,char**argv) {
   char fn[161],pstr[321],wmag[21]=" -wmag=2";
   FILE *f;
   int x,y,c,ALL=1;

   if (argc<2) {
      printf("Usage: %s <filter> <<options>>\n",*argv);
      printf("Supported options:\n");
      printf("   -hrc  compute HRC psf file only\n");
      return -1;
   }
   for (x=2;x<argc;x++) {
      if (!strcasecmp(argv[x],"-HRC")) ALL=0;
      else {
	 printf("Unknown option \"%s\"\n",argv[x]);
	 return -1;
      }
   }
   if (argv[1][strlen(argv[1])-1]=='N') wmag[0]=0;
   sprintf(fn,"%s/tmp",BASEDIR);
   mkdir(fn,00755);
   sprintf(fn,"%s/tmp/%s",BASEDIR,argv[1]);
   mkdir(fn,00755);
   sprintf(fn,"%s/tmp/%s/wfc.pos",BASEDIR,argv[1]);
   if ((f=fopen(fn,"w"))==NULL) {
      printf("Cannot write %s\n",fn);
      return -1;
   }
   for (y=128;y<2048;y+=256) for (x=128;x<4096;x+=256) fprintf(f,"%d %d\n",x,y);
   fclose(f);
   sprintf(fn,"%s/tmp/%s/hrc.pos",BASEDIR,argv[1]);
   if ((f=fopen(fn,"w"))==NULL) {
      printf("Cannot write %s\n",fn);
      return -1;
   }
   for (y=128;y<1024;y+=256) for (x=128;x<1024;x+=256) fprintf(f,"%d %d\n",x,y);
   fclose(f);
   for (c=0;c<3;c++) if (c==0 || ALL) {
      sprintf(fn,"%s/tmp/%s/tiny1.resp",BASEDIR,argv[1]);
      if ((f=fopen(fn,"w"))==NULL) {
	 printf("Cannot write %s\n",fn);
	 return -1;
      }
      if (c) fprintf(f,"15\n%d\n@%s/tmp/%s/wfc.pos\n",c,BASEDIR,argv[1]);
      else fprintf(f,"16\n@%s/tmp/%s/hrc.pos\n",BASEDIR,argv[1]);
      if (wmag[0]) fprintf(f,"%s\n1\n11\n3.0\n",argv[1]);
      else fprintf(f,"%s\n3.0\n",argv[1]);
      if (!c) fprintf(f,"hrc\n");
      else fprintf(f,"wfc%d\n",c);
      fclose(f);
      if (!c) sprintf(pstr,"cat %s/tmp/%s/tiny1.resp | %s/tiny1 %s/tmp/%s/hrc.par%s > %s/tmp/%s/tiny1.log",BASEDIR,argv[1],TTDIR,BASEDIR,argv[1],wmag,BASEDIR,argv[1]);
      else sprintf(pstr,"cat %s/tmp/%s/tiny1.resp | %s/tiny1 %s/tmp/%s/wfc%d.par%s >> %s/tmp/%s/tiny1.log",BASEDIR,argv[1],TTDIR,BASEDIR,argv[1],c,wmag,BASEDIR,argv[1]);
      if ((f=popen(pstr,"w"))!=NULL) pclose(f);
      else printf("Cannot run script %d\n",c);
   }
   sprintf(fn,"%s/tmp/%s/wfc.pos",BASEDIR,argv[1]);
   unlink(fn);
   sprintf(fn,"%s/tmp/%s/hrc.pos",BASEDIR,argv[1]);
   unlink(fn);
   sprintf(fn,"%s/tmp/%s/tiny1.resp",BASEDIR,argv[1]);
   unlink(fn);
   sprintf(fn,"%s/tmp/%s/runtiny",BASEDIR,argv[1]);
   if ((f=fopen(fn,"w"))==NULL) {
      printf("Cannot write %s\n",fn);
      return -1;
   }
   fprintf(f,"#! /bin/csh\n");
   fprintf(f,"rm tiny.log; touch tiny.log\n");
   fprintf(f,"tiny2 hrc.par >> tiny.log\n");
   if (ALL) {
      fprintf(f,"tiny2 wfc1.par >> tiny.log\n");
      fprintf(f,"tiny2 wfc2.par >> tiny.log\n");
   }
   for (x=0;x<16;x++) fprintf(f,"tiny3 hrc.par sub=3 pos=%d >> tiny.log\n",x);
   if (ALL) {
      for (x=0;x<128;x++) fprintf(f,"tiny3 wfc1.par sub=5 pos=%d >> tiny.log\n",x);
      for (x=0;x<128;x++) fprintf(f,"tiny3 wfc2.par sub=5 pos=%d >> tiny.log\n",x);
   }
   fclose(f);
   chmod(fn,00755);
   return 0;
}
