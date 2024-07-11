extern int WFPC2_DRIZZLE_BASE;

//cm: PC=0, WFC2=1, WFC3=2, WFC4=3;
//vmag: VEGA magnitudes;
//smag: UBVRI magnitudes;

//correct distortion;
extern void WFPC2shift(int img,double*x,double*y);

//add distortion;
extern void WFPC2unshift(int img,double*x,double*y);

//transforms, returns values in smag[];
extern void transform(int cm,double vmag[],double dvmag[],double smag[]);

//returns VEGAMAG;
extern double untransform(int cm,int filt,double smag[]);

//set parameters needed for WFPC2;
extern void wfpc2initparam(void);

// update photometry parameters for WFC plate scale;
extern void wfpc2radii(void);

// update alignment parameters for WFC plate scale;
extern void wfpc2tweakshift(void);

//initialize PSF libraries;
extern void wfpc2initpsf(void);

//dispose PSF libraries;
extern void wfpc2freepsf(void);

//calculate PSF from WFPC2 library;
extern int calcwfpc2psf(int img,float x,float y,int r,int force);

//do CTE corrections, photometric transformations, and output mag only;
extern double WFPC2calcmag(int img,float x0,float y0,float ct0,float bg,int useCTE);

//do CTE corrections, photometric transformations, and output;
extern void WFPC2outstar(FILE*,float x,float y,photdatatype*);
extern void WFPC2outstarimg(int img,FILE*,float x,float y,photdatatype*);

//output column headers and image descriptions
extern void WFPC2outstarinfo(FILE*,int*ct);
extern char *WFPC2imagestring(int img);

//aperture correction radius (0.5 arcsec WFC/0.3 arcsec HRC)
extern float wfpc2_apsize(int img,float x,float y);

//put lines in .info file;
extern void writewfpc2info(void);

//correct counts in fake stars;
extern void WFPC2readfakemag(FILE*);
extern void WFPC2fixfakemag(int img,float x0,float y0,double*ct0,float*bg);
