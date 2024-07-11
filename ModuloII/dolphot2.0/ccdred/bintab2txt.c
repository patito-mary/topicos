#include "../include/fits.h"
#include <unistd.h>

ftype fits;

int main(int argc,char**argv) {
   int ext,i,j,k;
   imtype *im;

   if (argc<3 || argc>4) {
      printf("Usage: %s <FITS file> <output> <<header>>\n",*argv);
      return -1;
   }
   readfits(argv[1],&fits,1);
   for (ext=0;ext<fits.Next && fits.ext[ext].bintabdata==0;ext++);
   if (ext>=fits.Next) {
      printf("ERROR: could not locate binary table\n");
      return -1;
   }
   im = fits.ext+ext;
   FILE *f=fopen(argv[2],"w");
   if (!f) {
      printf("ERROR: could not write \"%s\"\n",argv[2]);
      return -1;
   }
   for (i=0;i<im->Y;i++) {
      for (j=0;j<im->tfields;j++) {
	 fprintf(f," ");
	 for (k=0;k<im->bintabdata[j].length;k++) {
	    void *ptr = (void*) ( (char*)im->data[0][0] + i*im->X + im->bintabdata[j].offset + k*im->bintabdata[j].isize );
	    switch(im->bintabdata[j].type) {
	    case 'A':
	       if (*((char*)ptr)==0) fprintf(f," ");
	       else fprintf(f,im->bintabdata[j].format,*((char*)ptr));
	       break;
	    case 'B':
	       fprintf(f,im->bintabdata[j].format,*((int8_t*)ptr));
	       break;
	    case 'D':
	       fprintf(f,im->bintabdata[j].format,*((double*)ptr));
	       break;
	    case 'E':
	       fprintf(f,im->bintabdata[j].format,*((float*)ptr));
	       break;
	    case 'I':
	       fprintf(f,im->bintabdata[j].format,*((int16_t*)ptr));
	       break;
	    case 'J':
	       fprintf(f,im->bintabdata[j].format,*((int32_t*)ptr));
	       break;
	    case 'K':
	       fprintf(f,im->bintabdata[j].format,*((int64_t*)ptr));
	       break;
	    }
	    if (k<im->bintabdata[j].length-1 && im->bintabdata[j].type!='A') fprintf(f,",");
	 }
      }
      fprintf(f,"\n");
   }
   fclose(f);
   if (argc==4) {
      f=fopen(argv[3],"w");
      if (!f) {
	 printf("ERROR: could not write \"%s\"\n",argv[3]);
	 return -1;
      }
      for (i=0;i<im->tfields;i++) {
	 fprintf(f,"%d. %s\n",i+1,im->bintabdata[i].name);
      }
      fclose(f);
   }
   return 0;
}
