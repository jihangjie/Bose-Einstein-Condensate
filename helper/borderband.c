/*-----------------------------------------------------------------------*/
/*-----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*------------------------------------------------------------------------*/
#include"borderband.h"
#define fabs(x) ((x)>0.0 ?(x):-(x))


/* n=number of equations [0..N] -> N+1, etc... */

void band2lu(double **a,int n,int m1,int m2,double **al,int *indx)
{
        int i,j,k,l;
        int mm;
        double temp;

        mm=m1+m2+1;
        l=m1;
        for (i=1;i<=m1;i++)
        {
                for (j=m1+2-i;j<=mm;j++)
                        a[i-1][j-l-1]=a[i-1][j-1];
                l--;
                for (j=mm-l;j<=mm;j++)
                        a[i-1][j-1]=0.0;
        }
        l=m1;
        for (k=1;k<=n;k++)
        {
                temp=a[k-1][0];
                i=k;
                if (l < n) l++;
                for (j=k+1;j<=l;j++)
                {
                        if (fabs(a[j-1][0]) > fabs(temp))
                        {
                                temp=a[j-1][0];
                                i=j;
                        }
                }
                indx[k-1]=i;
                if (temp == 0.0)
                        a[k-1][0]=TINY;
                if (i != k)
                {
                        for (j=1;j<=mm;j++)
                                swap(a[k-1][j-1],a[i-1][j-1]);
                }
                for (i=k+1;i<=l;i++)
                {
                        temp=a[i-1][0]/a[k-1][0];
                        al[k-1][i-k-1]=temp;
                        for (j=2;j<=mm;j++)
                                a[i-1][j-2]=a[i-1][j-1]-temp*a[k-1][j-1];
                        a[i-1][mm-1]=0.0;
                }
        }
}

void bandsolve(double **a,int n,int m1,int m2,double **al,int *indx,double *b)
{
        int i,k,l;
        int mm;
        double temp;

        mm=m1+m2+1;
        l=m1;
        for (k=1;k<=n;k++)
        {
                i=indx[k-1];
                if (i != k)
                        swap(b[k-1],b[i-1]);
                if (l < n)
                        l++;
                for (i=k+1;i<=l;i++)
                        b[i-1] -= al[k-1][i-k-1]*b[k-1];
        }
        l=1;
        for (i=n;i>=1;i--)
        {
                temp=b[i-1];
                for (k=2;k<=l;k++)
                        temp -= a[i-1][k-1]*b[k+i-2];
                b[i-1]=temp/a[i-1][0];
                if (l < mm)
                        l++;
        }
}

void bandmult(double **a,int n,int m1,int m2,double *x,double *b)
{
        int i,j,k,tmploop;

        for (i=1;i<=n;i++)
        {
                k=i-m1-1;
                tmploop=min(m1+m2+1,n-k);
                b[i-1]=0.0;
                for (j=max(1,1-k);j<=tmploop;j++)
                {
                        b[i-1] += a[i-1][j-1]*x[j+k-1];
                }
        }
}
/*-----------------------------------------------------------*/


/* 
matrix A = [a(n*n, band matrix with bandwidth (m1+1+m2)) b1(n*1) b2(n*1)
		  				c1(1*n)							d11		d12
		  				c2(1*n)							d21		d22 ]
B = [f(n*1)
	 g(2*1)]

solve for x = A\B. 			

*/
void border2bandsolve(double **a, double *b1,double *b2, double *c1,double *c2, 
					double **d, int n, int m1,int m2, double **al,int *indx,double *f, double *g)
{
	double *x0;
	double *x1;
	double *x2;


	double alpha1,alpha2,alpha3,alpha4,beta1,beta2;
	int i,j;
	double tmp1,tmp2;	

	x0 = myalloc(n);
	x1 = myalloc(n);
	x2 = myalloc(n);
	

	band2lu(a,n,m1,m2,al,indx);
	bandsolve(a,n,m1,m2,al,indx,f);
	for (i=1;i<=n;i++)
		x0[i-1] = f[i-1];
		
	bandsolve(a,n,m1,m2,al,indx,b1);
	for (i=1;i<=n;i++)
		x1[i-1] = b1[i-1];
	

	bandsolve(a,n,m1,m2,al,indx,b2);
	for (i=1;i<=n;i++)
		x2[i-1] = b2[i-1];
	
	alpha1 = d[0][0];
	for (i=1;i<=n;i++)
		alpha1 -= c1[i-1]*x1[i-1];
	
	alpha2 = d[0][1];
	for (i=1;i<=n;i++)
		alpha2 -= c1[i-1]*x2[i-1];
		
	alpha3 = d[1][0];
	for (i=1;i<=n;i++)
		alpha3 -= c2[i-1]*x1[i-1];
	
	alpha4 = d[1][1];
	for (i=1;i<=n;i++)
		alpha4 -= c2[i-1]*x2[i-1];
	
	beta1 = g[0];
	for (i=1;i<=n;i++)
		beta1 -= c1[i-1]*x0[i-1];
		
	beta2 = g[1];
	for (i=1;i<=n;i++)
		beta2 -= c2[i-1]*x0[i-1];
	
	tmp1 = (alpha4*beta1-alpha2*beta2)/(alpha1*alpha4-alpha2*alpha3);
	tmp2 = (alpha1*beta2-alpha3*beta1)/(alpha1*alpha4-alpha2*alpha3);
	
	
	for(i=1;i<=n;i++)
		f[i-1] = x0[i-1]-tmp1*x1[i-1]-tmp2*x2[i-1];

	
	g[0] = tmp1;
	g[1] = tmp2;
	
	free(x0);
	free(x1);
	free(x2);

}

/*---------------------------------------------------------------------*/

/* 
matrix A = [a(n*n, band matrix with bandwidth (m1+1+m2)) b(n*1)
		  				c(1*n)							 d	]
B = [f(n*1)
	 g(1*1)]

solve for x = A\B. 			

*/

void border1bandsolve(double **a, double *b,double *c, 
					double d, int n, int m1,int m2, double **al,int *indx,double *f, double *g)
{
	double *x0;
	double *x1;

	int i,j;
	double tmp;
	double alpha,beta;	

	x0 = myalloc(n);
	x1 = myalloc(n);
	

	band2lu(a,n,m1,m2,al,indx);
	bandsolve(a,n,m1,m2,al,indx,f);
	for (i=1;i<=n;i++)
		x0[i-1] = f[i-1];
		
	bandsolve(a,n,m1,m2,al,indx,b);
	for (i=1;i<=n;i++)
		x1[i-1] = b[i-1];


	alpha = d;
	for (i=1;i<=n;i++)
		alpha -= c[i-1]*x1[i-1];
	
	beta = g[0];
	for (i=1;i<=n;i++)
		beta -= c[i-1]*x0[i-1];

	
	tmp = beta/alpha;

	for(i=1;i<=n;i++)
		f[i-1] = x0[i-1]-tmp*x1[i-1];

	g[0] = tmp;

	
	free(x0);
	free(x1);

}

/* 
matrix A = [a(n*n, band matrix with bandwidth (m1+1+m2)) b1(n*1) b2(n*1) b3(n*1)
		  				c1(1*n)							d11		d12		d13
		  				c2(1*n)							d21		d22 	d23
		  				c3(1*n)							d31		d32		d33 ]
B = [f(n*1)
	 g(3*1)]

solve for x = A\B. 			

*/
void border3bandsolve(double **a, double *b1,double *b2,double *b3, double *c1,double *c2, double *c3, 
					double **d, int n, int m1,int m2, double **al,int *indx,double *f, double *g)
{
	double *x0;
	double *x1;
	double *x2;
	double *x3;


	double alpha11,alpha12,alpha13,alpha21,alpha22, alpha23, alpha31, alpha32, alpha33, beta1,beta2,beta3;
	int i,j;
	double tmp1,tmp2,tmp3;	

	x0 = myalloc(n);
	x1 = myalloc(n);
	x2 = myalloc(n);
	x3 = myalloc(n);

	band2lu(a,n,m1,m2,al,indx);
	bandsolve(a,n,m1,m2,al,indx,f);
	for (i=1;i<=n;i++)
		x0[i-1] = f[i-1];
	
		
	bandsolve(a,n,m1,m2,al,indx,b1);
	for (i=1;i<=n;i++)
		x1[i-1] = b1[i-1];
	

	bandsolve(a,n,m1,m2,al,indx,b2);
	for (i=1;i<=n;i++)
		x2[i-1] = b2[i-1];
	
	bandsolve(a,n,m1,m2,al,indx,b3);
	for (i=1;i<=n;i++)
		x3[i-1] = b3[i-1];


	alpha11 = d[0][0];
	for (i=1;i<=n;i++)
		alpha11 -= c1[i-1]*x1[i-1];
	
	alpha12 = d[0][1];
	for (i=1;i<=n;i++)
		alpha12 -= c1[i-1]*x2[i-1];
	
	alpha13 = d[0][2];
	for (i=1;i<=n;i++)
		alpha13 -= c1[i-1]*x3[i-1];	
		
	alpha21 = d[1][0];
	for (i=1;i<=n;i++)
		alpha21 -= c2[i-1]*x1[i-1];
	
	alpha22 = d[1][1];
	for (i=1;i<=n;i++)
		alpha22 -= c2[i-1]*x2[i-1];
	
	alpha23 = d[1][2];
	for (i=1;i<=n;i++)
		alpha23 -= c2[i-1]*x3[i-1];
	
	alpha31 = d[2][0];
	for (i=1;i<=n;i++)
		alpha31 -= c3[i-1]*x1[i-1];
	
	alpha32 = d[2][1];
	for (i=1;i<=n;i++)
		alpha32 -= c3[i-1]*x2[i-1];
	
	alpha33 = d[2][2];
	for (i=1;i<=n;i++)
		alpha33 -= c3[i-1]*x3[i-1];
	
	beta1 = g[0];
	for (i=1;i<=n;i++)
		beta1 -= c1[i-1]*x0[i-1];
		
	beta2 = g[1];
	for (i=1;i<=n;i++)
		beta2 -= c2[i-1]*x0[i-1];
	
	beta3 = g[2];
	for (i=1;i<=n;i++)
		beta3 -= c3[i-1]*x0[i-1];
	
	tmp1 =  (alpha12*alpha23*beta3 - alpha12*alpha33*beta2 - alpha13*alpha22*beta3 + alpha13*alpha32*beta2 + alpha22*alpha33*beta1 - alpha23*alpha32*beta1)/(alpha11*alpha22*alpha33 - alpha11*alpha23*alpha32 -alpha12*alpha21*alpha33 + alpha12*alpha23*alpha31 + alpha13*alpha21*alpha32 - alpha13*alpha22*alpha31);
	tmp2 =  - (alpha11*alpha23*beta3 - alpha11*alpha33*beta2 - alpha13 *alpha21 *beta3 + alpha13* alpha31 *beta2 + alpha21 *alpha33* beta1 - alpha23 *alpha31 *beta1)/(alpha11*alpha22*alpha33 - alpha11*alpha23*alpha32 -  alpha12*alpha21 *alpha33 + alpha12* alpha23 *alpha31 + alpha13 *alpha21 *alpha32 - alpha13 *alpha22 *alpha31);
	tmp3 = (alpha11* alpha22* beta3 - alpha11 *alpha32 *beta2 - alpha12 *alpha21 *beta3 + alpha12 *alpha31 *beta2 + alpha21 *alpha32 *beta1 - alpha22 *alpha31 *beta1) /(alpha11 *alpha22 *alpha33 - alpha11 *alpha23 *alpha32 - alpha12* alpha21* alpha33 + alpha12 *alpha23 *alpha31 + alpha13 *alpha21 *alpha32 - alpha13* alpha22 *alpha31);
	
	
	for(i=1;i<=n;i++)
		f[i-1] = x0[i-1]-tmp1*x1[i-1]-tmp2*x2[i-1]-tmp3*x3[i-1];

	
	g[0] = tmp1;
	g[1] = tmp2;
	g[2] = tmp3;
	

	free(x0);
	free(x1);
	free(x2);
	free(x3);

}

/*--------------------------------------------------------------------------*/
