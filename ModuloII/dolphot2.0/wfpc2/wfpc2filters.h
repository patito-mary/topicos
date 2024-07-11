typedef struct {
   char name[11],color;
   double zp[4];
   char xorder[8];
   int xformc[8];
   double xform[8][5];
   //double idc[2][3][33];
} WFPC2filttype;

extern int WFPC2_NFILTERS;
extern WFPC2filttype *WFPC2filters;

extern int WFPC2findfilt(char*);
extern void WFPC2initfilters(void);

extern void WFPC2transform(double vmag[],double dvmag[],double smag[]);
extern double WFPC2untransform(int filt,double smag[5]);
extern double WFPC2_CTE(int chip,float x,float y,float cts,float ctmult,int ign,float sky,float epoch,int temp,int filt);
