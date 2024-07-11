typedef struct {
   char name[11],color;
   double zp[3];
   char xorder[2][8];
   int xformc[2][8];
   double xform[2][8][5];
   double idc[2][3][33];
   double dzp[3];
} ACSfilttype;

extern int ACS_NFILTERS;
extern ACSfilttype *ACSfilters;

extern int ACSfindfilt(char*);
extern void ACSinitfilters(void);

extern void ACStransform(int cm,double vmag[],double dvmag[],double smag[]);
extern double ACSuntransform(int cm,int filt,double smag[5]);
extern double ACS_CTE(int chip,float x,float y,float cts,float ctmult,float gn,float sky,float epoch);
extern double ACS_ZP(int filt,int cm,double MJD);
