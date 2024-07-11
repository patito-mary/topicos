typedef enum {
   NONE,
   WFPC2,
   ACS,
   WFC3,
   ROMAN,
   NIRCAM,
   NIRISS,
   MIRI
} instrument;

typedef struct {
   instrument inst;
   int cm;
   int filt;
   int i1; // temperature (WFPC2)
   int i2; // gain (WFPC2)
   double data1; // free field - used by NIRCAM for conversion from DN/sec to Jy
} hstmodetype;

typedef struct{float ct0,ct,dct,ctcorr,dctcorr,sky,m,dm,chi,sh,rnd,r,e,a,b,c,crowd;int flag;} photdatatype;

typedef double apsftype[3][6];
typedef float apskytype[2];
typedef float dpostype[40];
typedef float wcstype[80];
typedef double wcsreftype[4];
typedef int hstoffsettype[2];

#ifdef DOLPHOT_MAIN
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN hstmodetype *hstmode;
EXTERN apsftype *apsf;
EXTERN apskytype *apsky;
EXTERN dpostype *dpos,*dpos0,*ref2img;
EXTERN wcstype *wcs;
EXTERN wcsreftype *wcsref,*wcsminmax;
EXTERN int wcsminmax_set;
EXTERN hstoffsettype *hstoffset;
EXTERN imtype *datahd,*dataim;
EXTERN int Timg,Nimg;
EXTERN int *rphot,*RPSF,SubPixel,FPSF,EPSF;
EXTERN double *RAper,*RChi,MinS,MaxS,MaxE,PSFStep,Zero,*apsize;
EXTERN int PSFsol;
EXTERN int FlagMask,CombineChi,WFPC2useCTE,ACSuseCTE,WFC3useCTE,InterpPSFlib;
EXTERN double *iGAIN,*iEXP,*iEXP0,*iEPOCH,*apcor;
EXTERN int lastpsftype;
EXTERN float **psf;
EXTERN FILE *finfo;
EXTERN int poffreset;
EXTERN double *RSky0,*RSky1,*RSky20,*RSky21;
EXTERN int DEBUG;
EXTERN int DRIZZLE_BASE;
EXTERN int ACSpsfType;
EXTERN int WFC3psfType[3];
EXTERN int NIRCAMvega;
EXTERN int NIRISSvega;
EXTERN int MIRIvega;
