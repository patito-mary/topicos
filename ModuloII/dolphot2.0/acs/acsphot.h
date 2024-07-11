extern int ACS_DRIZZLE_BASE;

//cm: HRC=0, WFC1=1, WFC2=2;
//vmag: VEGA magnitudes;
//smag: UBVRI magnitudes;

//correct distortion;
extern void ACSshift(int img,double*x,double*y);

//add distortion;
extern void ACSunshift(int img,double*x,double*y);

//transforms, returns values in smag[];
extern void transform(int cm,double vmag[],double dvmag[],double smag[]);

//returns VEGAMAG;
extern double untransform(int cm,int filt,double smag[]);

//set parameters needed for ACS;
extern void acsinitparam(void);

//initialize PSF libraries;
extern void acsinitpsf(void);

//dispose PSF libraries;
extern void acsfreepsf(void);

//calculate PSF from ACS library;
extern int calcacspsf(int img,float x,float y,int r,int force);

//do CTE corrections, photometric transformations, and output mag only;
extern double ACScalcmag(int img,float x0,float y0,float ct0,float bg,int useCTE);

//do CTE corrections, photometric transformations, and output;
extern void ACSoutstar(FILE*,float x,float y,photdatatype*);
extern void ACSoutstarimg(int img,FILE*,float x,float y,photdatatype*);

//output column headers and image descriptions
extern void ACSoutstarinfo(FILE*,int*ct);
extern char *ACSimagestring(int img);

//aperture correction radius (0.5 arcsec WFC/0.3 arcsec HRC)
extern float acs_apsize(int img,float x,float y);

//put lines in .info file;
extern void writeacsinfo(void);

//correct counts in fake stars;
extern void ACSreadfakemag(FILE*);
extern void ACSfixfakemag(int img,float x0,float y0,double*ct0,float*bg);
