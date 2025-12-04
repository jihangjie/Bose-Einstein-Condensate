/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
/* swiped from numerical recipes (2nd) bandec.c, banbks.c */
#define swap(a,b) {double temp=(a);(a)=(b);(b)=temp;}
#define TINY 1.0e-20
#define max(a,b) ((a)>(b)? (a):(b))
#define min(a,b) ((b)>(a)? (a):(b))
void band2lu(double **a,int n,int m1,int m2,double **al,int *indx);
void bandsolve(double **a,int n,int m1,int m2,double **al,int *indx,double *b);
void bandmult(double **a,int n,int m1,int m2,double *x,double *b);
void border2bandsolve(double **a, double *b1,double *b2, double *c1,double *c2, 
					double **d, int n, int m1,int m2, double **al,int *indx,double *f, double *g);
void border1bandsolve(double **a, double *b,double *c, 
					double d, int n, int m1,int m2, double **al,int *indx,double *f, double *g);
void border3bandsolve(double **a, double *b1,double *b2,double *b3, double *c1,double *c2, double *c3, 
					double **d, int n, int m1,int m2, double **al,int *indx,double *f, double *g);


/*-------------------------------(Dynamic Memory allocation stuff)------*/
