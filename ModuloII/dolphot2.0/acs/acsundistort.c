#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/dolphot.h"
#include <assert.h>
#include "acsfilters.h"
#include "acsdistort.h"

void usage(char*exe) {
   fprintf(stderr,"Usage: %s <<flags>>\n",exe);
   fprintf(stderr,"Data given via stdin, in column format:\n");
   fprintf(stderr,"  default format: <camera> <filter> <X> <Y> <optional data>\n");
   fprintf(stderr,"  -file=#       specify filename\n");
   fprintf(stderr,"  -filt=#       specify filter\n");
   fprintf(stderr,"  -hrc          camera is HRC\n");
   fprintf(stderr,"  -wfc1         camera is WFC1\n");
   fprintf(stderr,"  -wfc2         camera is WFC2\n");
   fprintf(stderr,"  -wfcext       camera is WFC; format is raw dolphot output\n");
   fprintf(stderr,"  -hrcext       camera is WFC; format is raw dolphot output\n");
   exit(-1);
}

int main(int argc,char**argv) {
   int filtspec=0,mode=0;
   char*ptr,*ptr2,str[8001];
   int i,cm=-1,filt=-1;
   double x,y;
   FILE *fin=stdin;

   ACSinitfilters();
   for (i=1;i<argc;i++) {
      if (!strncasecmp(argv[i],"-file=",6)) {
	 fin=fopen(argv[i]+6,"r");
	 if (!fin) usage(*argv);
      }
      else if (!strncasecmp(argv[i],"-filt=",6)) {
	 filtspec=1;
	 filt=ACSfindfilt(argv[i]+6);
      }
      else if (!strcasecmp(argv[i],"-hrc")) {cm=0; mode=1;}
      else if (!strcasecmp(argv[i],"-hrcext")) {cm=0; mode=2;}
      else if (!strcasecmp(argv[i],"-wfc1")) {cm=1; mode=1;}
      else if (!strcasecmp(argv[i],"-wfc2")) {cm=1; mode=1;}
      else if (!strcasecmp(argv[i],"-wfcext")) {mode=3;}
      else usage(*argv);
   }
   if ((mode==2 || mode==3) && filtspec==0) {
      fprintf(stderr,"Error: if using raw dolphot output format, must specify filter\n");
      usage(*argv);
   }
   while (fgets(str,8001,fin)) {
      if (mode==2 || mode==3) {
	 i=strtol(str,&ptr,10);  // extension
	 if (ptr==str) {
	    fprintf(stderr,"No extension found:\n%s",str);
	    return -1;
	 }
	 if (mode==3) {
	    assert(i==1 || i==2);
	    cm=3-i;
	 }
	 i=strtol(ptr,&ptr2,10);  // chip
	 if (ptr==ptr2) {
	    fprintf(stderr,"No chip found:\n%s",str);
	    return -1;
	 }
	 assert(i==1);
	 ptr=ptr2;
      }
      else {
	 if (mode==0) {
	    cm=strtol(str,&ptr,10);
	    if (ptr==str) {
	       fprintf(stderr,"No camera found:\n%s",str);
	       return -1;
	    }
	    if (cm<0 || cm>2) {
	       fprintf(stderr,"Illegal camera (0, 1, 2 OK):\n%s",str);
	       return -1;
	    }
	 }
	 else ptr=str;
	 while (*ptr==' ') ptr++;
	 if (filtspec==0) {
	    for (ptr2=ptr;*ptr2 && *ptr2!=' ';ptr2++);
	    if (!*ptr2) {
	       fprintf(stderr,"No filter found:\n%s",str);
	       return -1;
	    }
	    *ptr2=0;
	    filt=ACSfindfilt(ptr);
	    *ptr2=' ';
	 }
	 ptr=ptr2;
      }
      x=strtod(ptr,&ptr2);
      if (ptr2==ptr) {
	 fprintf(stderr,"No X value found:\n%s",str);
	 return -1;
      }
      if (x<0 || (cm==0 && x>1024) || (cm>0 && x>4096)) {
	 fprintf(stderr,"Illegal X value (0-1024 HRC, 0-4096 WFC):\n%s",str);
	 return -1;
      }
      ptr=ptr2;
      y=strtod(ptr,&ptr2);
      if (ptr2==ptr) {
	 fprintf(stderr,"No Y value found:\n%s",str);
	 return -1;
      }
      if (y<0 || (cm==0 && y>1024) || (cm>0 && y>2048)) {
	 fprintf(stderr,"Illegal Y value (0-1024 HRC, 0-2048 WFC):\n%s",str);
	 return -1;
      }
      ACSfwddistort(cm,filt,&x,&y);
      if (cm==0) {
	 x=-x*ACSfilters[filt].idc[0][cm][4]+ACSfilters[filt].idc[0][cm][2];
	 y=y*ACSfilters[filt].idc[0][cm][4]+ACSfilters[filt].idc[0][cm][3];
      }
      else {
	 x=x*ACSfilters[filt].idc[0][cm][4]+ACSfilters[filt].idc[0][cm][2];
	 y=-y*ACSfilters[filt].idc[0][cm][4]+ACSfilters[filt].idc[0][cm][3];
      }
      printf("%8.3f %8.3f %s",x,y,str);
   }
   if (fin!=stdin) fclose(fin);
   return 0;
}
