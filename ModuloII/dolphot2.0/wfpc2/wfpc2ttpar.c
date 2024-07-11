#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc,char**argv) {
   char fn[161],pstr[321],wmag[21]=" -wmag=2";
   FILE *f;
   int x,y,c;

   if (argc<2) {
      printf("Usage: %s <filter> <<options>>\n",*argv);
      printf("Supported options:\n");
      return -1;
   }
   for (x=2;x<argc;x++) {
      printf("Unknown option \"%s\"\n",argv[x]);
      return -1;
   }
   if (argv[1][strlen(argv[1])-1]=='N') wmag[0]=0;
   sprintf(fn,"%s/tmp",BASEDIR);
   mkdir(fn,00755);
   sprintf(fn,"%s/tmp/%s_WFPC2",BASEDIR,argv[1]);
   mkdir(fn,00755);
   sprintf(fn,"%s/tmp/%s_WFPC2/pos.xy",BASEDIR,argv[1]);
   if ((f=fopen(fn,"w"))==NULL) {
      printf("Cannot write %s\n",fn);
      return -1;
   }
   for (y=50;y<800;y+=100) for (x=50;x<800;x+=100) fprintf(f,"%d %d\n",x,y);
   fclose(f);
   for (c=0;c<4;c++) {
      sprintf(fn,"%s/tmp/%s_WFPC2/tiny1.resp",BASEDIR,argv[1]);
      if ((f=fopen(fn,"w"))==NULL) {
	 printf("Cannot write %s\n",fn);
	 return -1;
      }
      if (c) fprintf(f,"5\n%d\n@%s/tmp/%s_WFPC2/pos.xy\n",c+1,BASEDIR,argv[1]);
      else fprintf(f,"6\n@%s/tmp/%s_WFPC2/pos.xy\n",BASEDIR,argv[1]);
      if (wmag[0]) fprintf(f,"%s\n1\n11\n3.0\n",argv[1]);
      else fprintf(f,"%s\n3.0\n",argv[1]);
      if (!c) fprintf(f,"y\n5\npc\n");
      else fprintf(f,"y\n9\nwfc%d\n",c+1);
      fclose(f);
      if (!c) sprintf(pstr,"cat %s/tmp/%s_WFPC2/tiny1.resp | %s/tiny1 %s/tmp/%s_WFPC2/pc.par%s > %s/tmp/%s_WFPC2/tiny1.log",BASEDIR,argv[1],TTDIR,BASEDIR,argv[1],wmag,BASEDIR,argv[1]);
      else sprintf(pstr,"cat %s/tmp/%s_WFPC2/tiny1.resp | %s/tiny1 %s/tmp/%s_WFPC2/wfc%d.par%s >> %s/tmp/%s_WFPC2/tiny1.log",BASEDIR,argv[1],TTDIR,BASEDIR,argv[1],c+1,wmag,BASEDIR,argv[1]);
      if ((f=popen(pstr,"w"))!=NULL) pclose(f);
      else printf("Cannot run script %d\n",c);
   }
   sprintf(fn,"%s/tmp/%s_WFPC2/pos.xy",BASEDIR,argv[1]);
   unlink(fn);
   sprintf(fn,"%s/tmp/%s_WFPC2/tiny1.resp",BASEDIR,argv[1]);
   unlink(fn);
   sprintf(fn,"%s/tmp/%s_WFPC2/runtiny",BASEDIR,argv[1]);
   if ((f=fopen(fn,"w"))==NULL) {
      printf("Cannot write %s\n",fn);
      return -1;
   }
   fprintf(f,"#! /bin/csh\n");
   fprintf(f,"rm tiny.log; touch tiny.log\n");
   fprintf(f,"tiny2 pc.par >> tiny.log\n");
   fprintf(f,"tiny2 wfc2.par >> tiny.log\n");
   fprintf(f,"tiny2 wfc3.par >> tiny.log\n");
   fprintf(f,"tiny2 wfc4.par >> tiny.log\n");
   fclose(f);
   chmod(fn,00755);
   return 0;
}
