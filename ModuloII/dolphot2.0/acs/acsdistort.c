#include <stdio.h>
#include "acsfilters.h"

void ACSfwddistort(int cm,int filt,double*x,double*y) {
   double x0,y0;
   x0=*x-ACSfilters[filt].idc[0][cm][0];
   y0=*y-ACSfilters[filt].idc[0][cm][1];
   *x=(ACSfilters[filt].idc[0][cm][5]*y0+ACSfilters[filt].idc[0][cm][6]*x0+ACSfilters[filt].idc[0][cm][7]*y0*y0+ACSfilters[filt].idc[0][cm][8]*x0*y0+ACSfilters[filt].idc[0][cm][9]*x0*x0+ACSfilters[filt].idc[0][cm][10]*y0*y0*y0+ACSfilters[filt].idc[0][cm][11]*x0*y0*y0+ACSfilters[filt].idc[0][cm][12]*x0*x0*y0+ACSfilters[filt].idc[0][cm][13]*x0*x0*x0+ACSfilters[filt].idc[0][cm][14]*y0*y0*y0*y0+ACSfilters[filt].idc[0][cm][15]*x0*y0*y0*y0+ACSfilters[filt].idc[0][cm][16]*x0*x0*y0*y0+ACSfilters[filt].idc[0][cm][17]*x0*x0*y0*x0+ACSfilters[filt].idc[0][cm][18]*x0*x0*x0*x0+ACSfilters[filt].idc[0][cm][2])/ACSfilters[filt].idc[0][cm][4];
   *y=(ACSfilters[filt].idc[0][cm][19]*y0+ACSfilters[filt].idc[0][cm][20]*x0+ACSfilters[filt].idc[0][cm][21]*y0*y0+ACSfilters[filt].idc[0][cm][22]*x0*y0+ACSfilters[filt].idc[0][cm][23]*x0*x0+ACSfilters[filt].idc[0][cm][24]*y0*y0*y0+ACSfilters[filt].idc[0][cm][25]*x0*y0*y0+ACSfilters[filt].idc[0][cm][26]*x0*x0*y0+ACSfilters[filt].idc[0][cm][27]*x0*x0*x0+ACSfilters[filt].idc[0][cm][28]*y0*y0*y0*y0+ACSfilters[filt].idc[0][cm][29]*x0*y0*y0*y0+ACSfilters[filt].idc[0][cm][30]*x0*x0*y0*y0+ACSfilters[filt].idc[0][cm][31]*x0*x0*y0*x0+ACSfilters[filt].idc[0][cm][32]*x0*x0*x0*x0+ACSfilters[filt].idc[0][cm][3])/ACSfilters[filt].idc[0][cm][4];
   return;
}

void ACSrevdistort(int cm,int filt,double*x,double*y) {
   double x0,y0;
   x0=(*x)*ACSfilters[filt].idc[0][cm][4]-ACSfilters[filt].idc[0][cm][2];
   y0=(*y)*ACSfilters[filt].idc[0][cm][4]-ACSfilters[filt].idc[0][cm][3];
   *x=ACSfilters[filt].idc[0][cm][0]+ACSfilters[filt].idc[1][cm][5]*y0+ACSfilters[filt].idc[1][cm][6]*x0+ACSfilters[filt].idc[1][cm][7]*y0*y0+ACSfilters[filt].idc[1][cm][8]*x0*y0+ACSfilters[filt].idc[1][cm][9]*x0*x0+ACSfilters[filt].idc[1][cm][10]*y0*y0*y0+ACSfilters[filt].idc[1][cm][11]*x0*y0*y0+ACSfilters[filt].idc[1][cm][12]*x0*x0*y0+ACSfilters[filt].idc[1][cm][13]*x0*x0*x0+ACSfilters[filt].idc[1][cm][14]*y0*y0*y0*y0+ACSfilters[filt].idc[1][cm][15]*x0*y0*y0*y0+ACSfilters[filt].idc[1][cm][16]*x0*x0*y0*y0+ACSfilters[filt].idc[1][cm][17]*x0*x0*y0*x0+ACSfilters[filt].idc[1][cm][18]*x0*x0*x0*x0;
   *y=ACSfilters[filt].idc[0][cm][1]+ACSfilters[filt].idc[1][cm][19]*y0+ACSfilters[filt].idc[1][cm][20]*x0+ACSfilters[filt].idc[1][cm][21]*y0*y0+ACSfilters[filt].idc[1][cm][22]*x0*y0+ACSfilters[filt].idc[1][cm][23]*x0*x0+ACSfilters[filt].idc[1][cm][24]*y0*y0*y0+ACSfilters[filt].idc[1][cm][25]*x0*y0*y0+ACSfilters[filt].idc[1][cm][26]*x0*x0*y0+ACSfilters[filt].idc[1][cm][27]*x0*x0*x0+ACSfilters[filt].idc[1][cm][28]*y0*y0*y0*y0+ACSfilters[filt].idc[1][cm][29]*x0*y0*y0*y0+ACSfilters[filt].idc[1][cm][30]*x0*x0*y0*y0+ACSfilters[filt].idc[1][cm][31]*x0*x0*y0*x0+ACSfilters[filt].idc[1][cm][32]*x0*x0*x0*x0;
   return;
}
