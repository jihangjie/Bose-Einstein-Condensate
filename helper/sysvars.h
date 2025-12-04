/*=============================================================*/
int N;
double L,dx,dt;
double err_tol;
double tend;

/*---------------------------*/

double *U;
double *U0;
double *U0_record;
double *U_PDE;

double *rhs;
double **J;
double **J1;
int *indx;



/*-------------------------*/


double hmax0,hmin0;
double mass0;
double mass0_PDE;


int flag_plot;

double dt,tt;
double tnext,dtnext;
int BAD;


/*=====================================================*/
