#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/dolphot.h"
#include "acsfilters.h"
#include "acsdistort.h"

void usage(char*exe) {
   printf("Usage: %s <CM> <filt> <X> <Y> <fwd/rev>\n",exe);
   printf("  CM: 0=HRC, 1/2=WFC\n");
   printf("  fwd = forward (raw image to distortion corrected\n");
   printf("  rev = reverse (distortion corrected to raw image\n");
   exit(-1);
}

int main(int argc,char**argv) {
   int cm,filt;
   double X,Y;

   if (argc!=6) usage(*argv);
   ACSinitfilters();
   cm=atoi(argv[1]);
   filt=ACSfindfilt(argv[2]);
   X=atof(argv[3]);
   Y=atof(argv[4]);
   if (!strcasecmp(argv[5],"fwd")) ACSfwddistort(cm,filt,&X,&Y);
   else if (!strcasecmp(argv[5],"rev")) ACSrevdistort(cm,filt,&X,&Y);
   else usage(*argv);
   printf("%f %f\n",X,Y);
   return 0;
}
