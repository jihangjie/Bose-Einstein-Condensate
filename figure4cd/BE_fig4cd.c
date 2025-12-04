/*=====================================*/
/* PDE simulation for a Bose-Einstein thin film-type equation 
by Roman Taranets, Marina Chugunova, and Hangjie Ji
x^beta*u_t + (x^alpha*(f(u)*u_xx-df(u)*u_x^2))_xx=0, 0<x<L

f(u) = u^4, df(u)=4*u^3, ddf(u) = 12*u^2
     
BCs: 
x^alpha*(f(u)*u_xx-df(u)*u_x^2) = 0 at x = 0 and x = L
(x^alpha*(f(u)*u_xx-df(u)*u_x^2))_x = 0 at x = 0 and x = L
*/

/*COMPILE WITH 

cc BE_fig4cd.c -lm -O3
*/

/* 2025.11.24 by H.Ji*/
/*=====================================*/
/*=====================================*/
#include<stdio.h>
#include<math.h>
#include"../helper/myalloc.h"
/*--------------------------------------------------------*/
#define sqr(x) ((x)*(x))
#define max(a,b) ((a)>(b)? (a):(b))
#define sgn(x) ((x)>=0.0 ? 1.0 : -1.0)
#define pi M_PI
/*------------------------------------------------------------------------*/
#define swap(a,b) {double temp=(a);(a)=(b);(b)=temp;}
#define TINY 1.0e-20
#include "../helper/borderband.h"
#include "../helper/borderband.c"
/*------------------------------------------------------------------------*/

#include "../helper/sysvars.h"

/*-------------------------*/

FILE *out;
FILE *out1;


int memory_setup();
int run();
int timestep(double tt);
int build_rhs_jac(double *u, double *u_old);

double ff(double u);
double dff(double u);
double ddff(double u);

double param_beta;
double param_alpha;


double prewet_thickness=0.01;

double prewet_eps = 1e-6; 

int plot_idx = 0;

/*---------------*/

int main(int argc,char *argv[])
{
	int i,j;
	

	double M0;
	double h0;
	
	char text[100];
	double xx, xxp;
	
	

/*-------------------------*/
 	out=fopen("out_pde.dat","w");
	out1=fopen("out1_pde.dat","w");
	

	flag_plot = 1;
	

	/* -----System parameters -----*/
	param_alpha = 8.0;  // physical: alpha = 13.0/2.0
	param_beta = 1.0/2.0; 
	/*========================================================================*/

  	
  	/* domain from 0 to L */
  	
  	N = 10000;
	L = 10.0; 
	
	dx=L/N;
	memory_setup();
		
	tt=0.0;

     
    printf("\n alpha=%g, beta=%g\n",
    param_alpha, param_beta);
    /*====================== set ICs ================*/
    /*******************************************************************/ 	

   
 	
 	
 	for (i = 0; i<N+1; i++)
 	{
 	    xx = i*dx;
 	   
 	   U0[i] = prewet_thickness+0.4*exp(-10.0*pow(xx-0.3*L,2.0)) + 0.4*exp(-10.0*pow(xx-0.7*L,2.0));  // two-hump initial condition

 	}

	dt=0.001;
	tnext = 0.0;
	dtnext=10.0; // 5
	tend=300;
	
	plot_idx += 1;
	
    for(i=0;i<N+1;++i)
    {			
        fprintf(out,"%g\t %g\t %g\n",
        i*dx,U0[i],tt);
    }
    fprintf(out,"\n\n");
    fflush(out);
	
	
	/* mass of the initial condition*/
	mass0 = 0.0;
	
	for(i=0;i<N;++i)
	{
	    xx = i*dx;
	    xxp = (i+1)*dx;
		mass0+=0.5*dx*(U0[i]*pow(xx,param_beta)+U0[i+1]*pow(xxp,param_beta));
	}



    printf("\n\nPDE simulation for the Bose-Einstein thin film model with\n");

    printf("============================================\n\n");
    

    
    // run PDE solver
	run(); 

}

/*******************************************************************/


int run()
{

    int i,b_count=0;

    double umin,umin0,umax, umax0;
    int imin,imax;

    double mass,ht,h_xx_max;
    double h_x_L2;
    double flux;
    int ctl_idx; // contact line position
    int ctl_idx_left; // left contact line position

    double maxerr;
    
    double xx, xxp;


    tt=0.0;
    dt=1e-7;

    umin0=20.0;
    umax0 = 0.0;
    umax = 0.0;
    

    for(i=0;i<N+1;++i)
    {	
        umax0 = max(umax0, U0[i]);		
    }
    
    for(i=0;i<N+1;++i)
    {	
        umin0 = min(umin0, U0[i]);		
    }

/*====================== Evolution starts ================*/


	while(tt<tend)
	{
		umin0=umin;
		BAD=0;
	
	
		if(timestep(tt))
		{
		
		    
			tt+=dt;
			dt*=1.01;
			dt = min(0.01,dt);

			b_count=0; /* reset badstep counter */
			
			

		}
		else
		{
			dt*=0.5;
			BAD=1;
			b_count++;
			printf("BAD=%d\t dt=%g\n",BAD,dt);
			if((b_count>8) || (dt<1e-8) )
			{
				printf("stop. something's wrong.\n");
			
				exit(1);
			}
		
			continue; /* restart loop from top */
		}
	
	    
		imax=N;
		umax=U[N];
		for(i=1;i<N+1;++i)
		{			
			if(U[i]>umax)
			{
				imax=i;
				umax=U[i];
			}
		}
		imin=0;
		umin=umax;	
		for(i=1;i<N+1;++i)
		{
			if(U[i]<umin)
			{
				imin=i; 
				umin=U[i];
			}
		}
		
		// Mass 
		mass=0.0;
        for(i=0;i<N;++i)
        {
            xx = i*dx;
            xxp = (i+1)*dx; 
                 mass+=0.5*dx*(U0[i]*pow(xx,param_beta)+U0[i+1]*pow(xxp,param_beta));
        }
		
		// right contact line
		ctl_idx = 0;
        for(i=(int)(3.0/dx);i<(int)(7.0/dx);++i)
		{
		    if (U[i] < prewet_thickness + prewet_eps)
		    {
		        ctl_idx = i;
		        break;
		    }
		}
		
	    // left contact line
		ctl_idx_left = 0;
		
for(i=(int)(7.0/dx);i>(int)(3.0/dx);--i)
		{
		    if (U[i] < prewet_thickness + prewet_eps)
		    {
		        ctl_idx_left = i;
		        break;
		    }
		}
		
		h_xx_max = (U[imax+1]-2.0*U[imax]+U[imax-1])/(dx*dx);
		
		
		if(flag_plot)
		{
		fprintf(out1,"%g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\t %g\n",
			tt,umax, dx*imax, h_xx_max, dx*imin, umin, mass, ctl_idx*dx, ctl_idx_left*dx);
			fflush(out1);
		}
			
		
		if((tt>tnext) || (umin<umin0*0.9) || (umin > umin0*1.1) || (umax<umax0*0.9) || (umax > umax0*1.1))
		{		

			tnext=tt+dtnext;

			if(flag_plot)
			{
				for(i=0;i<N+1;++i)
				{			
		
					fprintf(out,"%g\t %g\t %g\t %g\t %g\n",
					i*dx,U[i],tt,umax, dx*imax);
				}
				fprintf(out,"\n\n");
				fflush(out);
			}

			printf("%d\t %g \t %g \t %g \t %g \t  %g\n",plot_idx, tt,dt, mass,U[0],umax);
			umax0 = umax;
			umin0 = umin;
			
			plot_idx += 1;
		}
	}


	return(0);

}



int timestep(double tt) /* return 1=good, 0=bad */
{
    int i,j;
	double maxerr;

	for(i=0;i<N+1;++i)
		U[i]=U0[i];

	err_tol=1e-8;

	for(j=0;j<10;++j)
	{
		maxerr=0.0;

		build_rhs_jac(U,U0);
		
		band2lu(J,N+1,2,2,J1,indx);
		bandsolve(J,N+1,2,2,J1,indx,rhs);
		
			
		for(i=0;i<N+1;++i)
		{
			maxerr=max(maxerr,fabs(rhs[i]));
			U[i]+=rhs[i];
			
		}
		

		if(maxerr<err_tol/10.0)
			break;
	
	}
	
	
	if(maxerr<err_tol) 
	{
		for(i=0;i<N+1;++i)
		{
			U0_record[i] = U0[i];
			U0[i]=U[i];  /* store previous step */
		}
		return(1);
	}
	else
	{
		printf("maxerr=%g\n",maxerr);
		return(0);
	}

}




/*--------------------------------------------------------------------------*/

int build_rhs_jac(double *u, double *u_old)
{
	int i,j,ni;

	double u0,up,upp,um,umm,uold,uoldm;
	double x0, xp, xm;

	#include"../helper/tempvars.h"

	for(i=0;i<N+1;++i) /* zero to start with */
	{	
		rhs[i]=0.0;
		for(j=0;j<5;++j)
			J[i][j]=0.0;
	}

    for(i=0;i<N+1;++i)
    {
    /*-------------------------*/
    /*-------------------------*/
        ni=N-i;

        u0=u[i];
        up=u[i+1];
        upp=u[i+2];
        um=u[i-1];
        umm=u[i-2];
        uold=u_old[i];
        
        x0=i*dx;
		xp=x0+dx;
		xm=x0-dx;
		
		
        /*-------------------------*/
        t1 = pow(x0,1.0*param_beta);
      t6 = pow(xp,1.0*param_alpha);
      t7 = ff(up);
      t11 = dx*dx;
      t12 = 1/t11;
      t14 = dff(up);
      t16 = pow(upp-u0,2.0);
      t22 = pow(x0,1.0*param_alpha);
      t23 = ff(u0);
      t28 = dff(u0);
      t30 = pow(up-um,2.0);
      t37 = pow(xm,1.0*param_alpha);
      t38 = ff(um);
      t43 = dff(um);
      t45 = pow(u0-umm,2.0);
      t53 = 1/dt*(u0-uold)*t1+t12*((t12*(upp-2.0*up+u0)*t7-t12*t16*t14/4.0)*t6
-2.0*(t12*(up-2.0*u0+um)*t23-t12*t30*t28/4.0)*t22+(t12*(u0-2.0*um+umm)*t38-t12*
t45*t43/4.0)*t37);

        rhs[i]= t53;
        
        /*-------------------------*/
        
        
        /*--------------*/
        t1 = pow(xm,1.0*param_alpha);
      t2 = ff(um);
      t3 = dx*dx;
      t4 = 1/t3;
      t6 = dff(um);
      t13 = t4*(t4*t2+t4*(u0-umm)*t6/2.0)*t1;
        
        J[i][0]=t13;
        
        /*--------------*/
         t1 = pow(x0,1.0*param_alpha);
      t2 = ff(u0);
      t3 = dx*dx;
      t4 = 1/t3;
      t6 = dff(u0);
      t14 = pow(xm,1.0*param_alpha);
      t15 = dff(um);
      t20 = ff(um);
      t23 = ddff(um);
      t25 = pow(u0-umm,2.0);
      t32 = t4*(-2.0*(t4*t2+t4*(up-um)*t6/2.0)*t1+(t4*(u0-2.0*um+umm)*t15-2.0*
t4*t20-t4*t25*t23/4.0)*t14);

        
        J[i][1]=t32;
        
        /*--------------*/
       t1 = pow(x0,1.0*param_beta);
      t4 = pow(xp,1.0*param_alpha);
      t5 = ff(up);
      t6 = dx*dx;
      t7 = 1/t6;
      t9 = dff(up);
      t16 = pow(x0,1.0*param_alpha);
      t17 = dff(u0);
      t22 = ff(u0);
      t25 = ddff(u0);
      t27 = pow(up-um,2.0);
      t34 = pow(xm,1.0*param_alpha);
      t35 = ff(um);
      t37 = dff(um);
      t46 = 1/dt*t1+t7*((t7*t5+t7*(upp-u0)*t9/2.0)*t4-2.0*(t7*(up-2.0*u0+um)*
t17-2.0*t7*t22-t7*t27*t25/4.0)*t16+(t7*t35-t7*(u0-umm)*t37/2.0)*t34);

        J[i][2]=t46;
        
        /*--------------*/
        t1 = pow(xp,1.0*param_alpha);
      t2 = dff(up);
      t6 = dx*dx;
      t7 = 1/t6;
      t9 = ff(up);
      t12 = ddff(up);
      t14 = pow(upp-u0,2.0);
      t20 = pow(x0,1.0*param_alpha);
      t21 = ff(u0);
      t23 = dff(u0);
      t32 = t7*((t7*(upp-2.0*up+u0)*t2-2.0*t7*t9-t7*t14*t12/4.0)*t1-2.0*(t7*t21
-t7*(up-um)*t23/2.0)*t20);

        
        J[i][3]=t32;
        
        /*--------------*/
         t1 = pow(xp,1.0*param_alpha);
      t2 = ff(up);
      t3 = dx*dx;
      t4 = 1/t3;
      t6 = dff(up);
      t13 = t4*(t4*t2-t4*(upp-u0)*t6/2.0)*t1;

    
        J[i][4]=t13;

        /*-------------------------*/
        /*-------------------------*/

        // BC at x = 0
        if (i == 0)
        {
            t1 = pow(x0,1.0*param_beta);
            t6 = pow(xp,1.0*param_alpha);
            t7 = ff(up);
            t11 = dx*dx;
            t12 = 1/t11;
            t14 = dff(up);
            t16 = pow(upp-u0,2.0);
            t24 = 1/dt*(u0-uold)*t1+2.0*t12*(t12*(upp-2.0*up+u0)*t7-t12*t16*t14/4.0)*
            t6;
            
            
            rhs[i]= t24;
            
            J[i][0] = 0.0;
            J[i][1] = 0.0;
            
            t1 = pow(x0,1.0*param_beta);
            t4 = pow(xp,1.0*param_alpha);
            t5 = ff(up);
            t6 = dx*dx;
            t7 = 1/t6;
            t9 = dff(up);
            t18 = 1/dt*t1+2.0*t7*(t7*t5+t7*(upp-u0)*t9/2.0)*t4;
            
            
            J[i][2] = t18;
            
            t1 = pow(xp,1.0*param_alpha);
            t2 = dff(up);
            t6 = dx*dx;
            t7 = 1/t6;
            t9 = ff(up);
            t12 = ddff(up);
            t14 = pow(upp-u0,2.0);
            t21 = 2.0*t7*(t7*(upp-2.0*up+u0)*t2-2.0*t7*t9-t7*t14*t12/4.0)*t1;
            
            
            J[i][3] = t21;
            
            t1 = pow(xp,1.0*param_alpha);
            t2 = ff(up);
            t3 = dx*dx;
            t4 = 1/t3;
            t6 = dff(up);
            t14 = 2.0*t4*(t4*t2-t4*(upp-u0)*t6/2.0)*t1;
            
            
            J[i][4] = t14;
        
        }
        
        // BC at x = 0
        if (i == 1)
        {
            t1 = pow(x0,1.0*param_beta);
            t6 = pow(xp,1.0*param_alpha);
            t7 = ff(up);
            t11 = dx*dx;
            t12 = 1/t11;
            t14 = dff(up);
            t16 = pow(upp-u0,2.0);
            t22 = pow(x0,1.0*param_alpha);
            t23 = ff(u0);
            t28 = dff(u0);
            t30 = pow(up-um,2.0);
            t39 = 1/dt*(u0-uold)*t1+t12*((t12*(upp-2.0*up+u0)*t7-t12*t16*t14/4.0)*t6
            -2.0*(t12*(up-2.0*u0+um)*t23-t12*t30*t28/4.0)*t22);
            
            rhs[i] = t39;
            
            J[i][0] = 0;
            
            t1 = pow(x0,1.0*param_alpha);
            t2 = ff(u0);
            t3 = dx*dx;
            t4 = 1/t3;
            t6 = dff(u0);
            t15 = -2.0*t4*(t4*t2+t4*(up-um)*t6/2.0)*t1;
            
            J[i][1] = t15;
            
            t1 = pow(x0,1.0*param_beta);
            t4 = pow(xp,1.0*param_alpha);
            t5 = ff(up);
            t6 = dx*dx;
            t7 = 1/t6;
            t9 = dff(up);
            t16 = pow(x0,1.0*param_alpha);
            t17 = dff(u0);
            t22 = ff(u0);
            t25 = ddff(u0);
            t27 = pow(up-um,2.0);
            t36 = 1/dt*t1+t7*((t7*t5+t7*(upp-u0)*t9/2.0)*t4-2.0*(t7*(up-2.0*u0+um)*
            t17-2.0*t7*t22-t7*t27*t25/4.0)*t16);
            
            J[i][2] = t36;
            
            t1 = pow(xp,1.0*param_alpha);
            t2 = dff(up);
            t6 = dx*dx;
            t7 = 1/t6;
            t9 = ff(up);
            t12 = ddff(up);
            t14 = pow(upp-u0,2.0);
            t20 = pow(x0,1.0*param_alpha);
            t21 = ff(u0);
            t23 = dff(u0);
            t32 = t7*((t7*(upp-2.0*up+u0)*t2-2.0*t7*t9-t7*t14*t12/4.0)*t1-2.0*(t7*t21
            -t7*(up-um)*t23/2.0)*t20);
            
            
            J[i][3] = t32;
            
            t1 = pow(xp,1.0*param_alpha);
            t2 = ff(up);
            t3 = dx*dx;
            t4 = 1/t3;
            t6 = dff(up);
            t13 = t4*(t4*t2-t4*(upp-u0)*t6/2.0)*t1;
            
            
            J[i][4] = t13;
        
        }
        
        // BC at x = L
        if (i == N-1)
        {
            t1 = pow(x0,1.0*param_beta);
            t6 = pow(x0,1.0*param_alpha);
            t7 = ff(u0);
            t11 = dx*dx;
            t12 = 1/t11;
            t14 = dff(u0);
            t16 = pow(up-um,2.0);
            t23 = pow(xm,1.0*param_alpha);
            t24 = ff(um);
            t29 = dff(um);
            t31 = pow(u0-umm,2.0);
            t39 = 1/dt*(u0-uold)*t1+t12*(-2.0*(t12*t7*(up-2.0*u0+um)-t12*t16*t14/4.0)
            *t6+(t12*(u0-2.0*um+umm)*t24-t12*t31*t29/4.0)*t23);
            
            
            rhs[i] = t39;
            
            t1 = pow(xm,1.0*param_alpha);
            t2 = ff(um);
            t3 = dx*dx;
            t4 = 1/t3;
            t6 = dff(um);
            t13 = t4*(t4*t2+t4*(u0-umm)*t6/2.0)*t1;
            
            J[i][0] = t13;
            
            t1 = pow(x0,1.0*param_alpha);
            t2 = ff(u0);
            t3 = dx*dx;
            t4 = 1/t3;
            t6 = dff(u0);
            t14 = pow(xm,1.0*param_alpha);
            t15 = dff(um);
            t20 = ff(um);
            t23 = ddff(um);
            t25 = pow(u0-umm,2.0);
            t32 = t4*(-2.0*(t4*t2+t4*(up-um)*t6/2.0)*t1+(t4*(u0-2.0*um+umm)*t15-2.0*
            t4*t20-t4*t25*t23/4.0)*t14);
            
            
            J[i][1] = t32;
            
            t1 = pow(x0,1.0*param_beta);
            t4 = pow(x0,1.0*param_alpha);
            t5 = dff(u0);
            t9 = dx*dx;
            t10 = 1/t9;
            t12 = ff(u0);
            t15 = ddff(u0);
            t17 = pow(up-um,2.0);
            t24 = pow(xm,1.0*param_alpha);
            t25 = ff(um);
            t27 = dff(um);
            t36 = 1/dt*t1+t10*(-2.0*(t10*(up-2.0*u0+um)*t5-2.0*t12*t10-t10*t17*t15/
            4.0)*t4+(t10*t25-t10*(u0-umm)*t27/2.0)*t24);
            
            
            J[i][2] = t36;
            
            t1 = pow(x0,1.0*param_alpha);
            t2 = ff(u0);
            t3 = dx*dx;
            t4 = 1/t3;
            t6 = dff(u0);
            t15 = -2.0*t4*(t4*t2-t4*(up-um)*t6/2.0)*t1;
            
            
            J[i][3] = t15;
            
            J[i][4] = 0;
        
        }
        
        // BC at x = L
        if (i == N)
        {
            t1 = pow(x0,1.0*param_beta);
            t6 = pow(xm,1.0*param_alpha);
            t7 = ff(um);
            t11 = dx*dx;
            t12 = 1/t11;
            t14 = dff(um);
            t16 = pow(u0-umm,2.0);
            t24 = 1/dt*(u0-uold)*t1+2.0*t12*(t12*t7*(u0-2.0*um+umm)-t12*t16*t14/4.0)*
            t6;
            
            rhs[i] = t24;
            
            t1 = pow(xm,1.0*param_alpha);
            t2 = ff(um);
            t3 = dx*dx;
            t4 = 1/t3;
            t6 = dff(um);
            t14 = 2.0*t4*(t4*t2+t4*(u0-umm)*t6/2.0)*t1;
            
            
            J[i][0] = t14;
            
            t1 = pow(xm,1.0*param_alpha);
            t2 = dff(um);
            t6 = dx*dx;
            t7 = 1/t6;
            t9 = ff(um);
            t12 = ddff(um);
            t14 = pow(u0-umm,2.0);
            t21 = 2.0*t7*(t7*(u0-2.0*um+umm)*t2-2.0*t7*t9-t7*t14*t12/4.0)*t1;
            
            
            
            J[i][1] = t21;
            
            t1 = pow(x0,1.0*param_beta);
            t4 = pow(xm,1.0*param_alpha);
            t5 = ff(um);
            t6 = dx*dx;
            t7 = 1/t6;
            t9 = dff(um);
            t18 = 1/dt*t1+2.0*t7*(t7*t5-t7*(u0-umm)*t9/2.0)*t4;
            
            
            
            J[i][2] = t18;
            
            J[i][3] = 0;
            
            J[i][4] = 0;
        
        }
    }
    

	for(i=0;i<N+1;++i)
		rhs[i]*= -1.0;

	/*-------------------------*/


	return(0);

}

int memory_setup()
{
        U=myalloc(N+1);
        U0=myalloc(N+1);
        U0_record=myalloc(N+1);
        U_PDE = myalloc(N+1);

        rhs=myalloc(N+1);
        J=matrix(N+1,5);
        J1=matrix(N+1,5);
        indx=myalloci(N+1);
        
        return(0);
}

double ff(double u)
{
    return (pow(u,4.0));

}
double dff(double u)
{
    return (4.0*pow(u,3.0));
}

double ddff(double u)
{
    return (12.0*pow(u,2.0));
}

