#include "../include/fits.h"

ftype fitsin,fitsout;

int main(int argc,char**argv) {
   int C,X0,X1,Y0,Y1,x,y;

   if (argc!=4 && argc!=8) {
      printf("Usage: %s <input> <output> <chip> <<X0> <X1> <Y0> <Y1>>\n",*argv);
      return 1;
   }
   readfits(argv[1],&fitsin,1);

   // set up output FITS
   C=atoi(argv[3])-1;
   if (C<0 || C>=fitsin.img.Z) {
      fprintf(stderr,"Chip must be 1-%d\n",fitsin.img.Z);
      exit(-1);
   }
   if (argc==4) {
      X0 = Y0 = 0;
      X1 = fitsin.img.X;
      Y1 = fitsin.img.Y;
   }
   else {
      X0 = atoi(argv[4]);
      X1 = atoi(argv[5]);
      Y0 = atoi(argv[6]);
      Y1 = atoi(argv[7]);
      if (X0<0 || X1<=X0 || X1>fitsin.img.X || Y0<0 || Y1<=Y0 || Y1>fitsin.img.Y) {
	 fprintf(stderr,"Error in X,Y limits\n");
	 exit(-1);
      }
   }
   imcopy(&(fitsin.img),&(fitsout.img));
   freeimg(fitsout.img.data,fitsout.img.X,fitsout.img.Y,fitsout.img.Z);
   fitsout.img.X = X1-X0;
   fitsout.img.Y = Y1-Y0;
   fitsout.img.Z = 1;
   fitsout.img.data = allocimg(fitsout.img.X,fitsout.img.Y,fitsout.img.Z);
   for (y=Y0;y<Y1;y++) {
      for (x=X0;x<X1;x++) {
	 fitsout.img.data[0][y-Y0][x-X0] = fitsin.img.data[C][y][x];
      }
   }
   fitsout.Next = 0;

   // output
   writefits(argv[2],&fitsout,1);

   return 0;
}
