#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/dolphot.h"
#include "wfpc2distort.h"

void usage(char*exe) {
   printf("Usage: %s <CM> <X> <Y> <fwd/rev>\n",exe);
   printf("  CM: 0=PC, 1-3=WFC2-4\n");
   printf("  fwd = forward (raw image to distortion corrected\n");
   printf("  rev = reverse (distortion corrected to raw image\n");
   exit(-1);
}

int main(int argc,char**argv) {
   int cm;
   double X,Y;

   if (argc!=5) usage(*argv);
   cm=atoi(argv[1]);
   X=atof(argv[2]);
   Y=atof(argv[3]);
   if (!strcasecmp(argv[4],"fwd")) WFPC2fwddistort(cm,&X,&Y);
   else if (!strcasecmp(argv[4],"rev")) WFPC2revdistort(cm,&X,&Y);
   else usage(*argv);
   printf("%f %f\n",X,Y);
   return 0;
}
