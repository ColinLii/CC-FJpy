//
//  CCFJ-py
//  A Python Package for seismic ambient noise cross-correlation and Frequency-Bessel Transform Method
//
//  GNU General Public License, Version 3, 29 June 2007
//
//  Copyright © 2021 Xiaofei Chen Research Group,
//  Department of Earth and Space Sciences,
//  Southern University of Science and Technology, China.
//
#include "FJcpu.h"
#ifdef useOMP
#include <omp.h>
#endif
#include <math.h>
#define M_PI 3.14159265358

using namespace std;

static double bessj0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          0 at input x                                      */
/*------------------------------------------------------------*/
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
         +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
      ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
         +y*(59272.64853+y*(267.8532712+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6
         -y*0.934935152e-7)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
   }
   return ans;
}

static double bessj1( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of first kind and order  */
/*          1 at input x                                      */
/*------------------------------------------------------------*/
{
   double ax,z;
   double xx,y,ans,ans1,ans2;

   if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
         +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
      ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
         +y*(99447.43394+y*(376.9991397+y*1.0))));
      ans=ans1/ans2;
   } else {
      z=8.0/ax;
      y=z*z;
      xx=ax-2.356194491;
      ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
         +y*(0.2457520174e-5+y*(-0.240337019e-6))));
      ans2=0.04687499995+y*(-0.2002690873e-3
         +y*(0.8449199096e-5+y*(-0.88228987e-6
         +y*0.105787412e-6)));
      ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
      if (x < 0.0) ans = -ans;
   }
   return ans;
}

static double bessy0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of second kind and order */
/*          0 at input x.                                     */
/*------------------------------------------------------------*/
{
   double z;
   double xx,y,ans,ans1,ans2;

   if (x < 8.0) {
      y=x*x;
      ans1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6
         +y*(10879881.29+y*(-86327.92757+y*228.4622733))));
      ans2=40076544269.0+y*(745249964.8+y*(7189466.438
         +y*(47447.26470+y*(226.1030244+y*1.0))));
      ans=(ans1/ans2)+0.636619772*bessj0(x)*log(x);
   } else {
      z=8.0/x;
      y=z*z;
      xx=x-0.785398164;
      ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
         +y*(-0.2073370639e-5+y*0.2093887211e-6)));
      ans2 = -0.1562499995e-1+y*(0.1430488765e-3
         +y*(-0.6911147651e-5+y*(0.7621095161e-6
         +y*(-0.934945152e-7))));
      ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
   }
   return ans;
}



static double bessy1( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel function of second kind and order */
/*          1 at input x.                                     */
/*------------------------------------------------------------*/
{
   double z;
   double xx,y,ans,ans1,ans2;

   if (x < 8.0) {
      y=x*x;
      ans1=x*(-0.4900604943e13+y*(0.1275274390e13
         +y*(-0.5153438139e11+y*(0.7349264551e9
         +y*(-0.4237922726e7+y*0.8511937935e4)))));
      ans2=0.2499580570e14+y*(0.4244419664e12
         +y*(0.3733650367e10+y*(0.2245904002e8
         +y*(0.1020426050e6+y*(0.3549632885e3+y)))));
      ans=(ans1/ans2)+0.636619772*(bessj1(x)*log(x)-1.0/x);
   } else {
      z=8.0/x;
      y=z*z;
      xx=x-2.356194491;
      ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
         +y*(0.2457520174e-5+y*(-0.240337019e-6))));
      ans2=0.04687499995+y*(-0.2002690873e-3
         +y*(0.8449199096e-5+y*(-0.88228987e-6
         +y*0.105787412e-6)));
      ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
   }
   return ans;
}

static double STVH0(double X) {
/*      =============================================
!       Purpose: Compute Struve function H0(x)
!       Input :  x   --- Argument of H0(x) ( x ò 0 )
!       Output:  SH0 --- H0(x)
!       ============================================= */
        double A0,BY0,P0,PI,Q0,R,S,T,T2,TA0;
	int K, KM;
        double SH0;
	    
	PI=3.141592653589793;
        S=1.0;
        R=1.0;
        if (X <= 20.0) {
           A0=2.0*X/PI;
           for (K=1; K<61; K++) {
              R=-R*X/(2.0*K+1.0)*X/(2.0*K+1.0);
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e15;
           }
e15:       SH0=A0*S;
        }
        else {
           KM=int(0.5*(X+1.0));
           if (X >= 50.0) KM=25;
           for (K=1; K<=KM; K++) {
			  //R=-R*pow((2.0*K-1.0)/X,2);
			  R = -R*(2.0*K-1.0)*(2.0*K-1.0)/X/X;
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e25;
           }
e25:       T=4.0/X;
           T2=T*T;
           P0=((((-.37043e-5*T2+.173565e-4)*T2-.487613e-4)*T2+.17343e-3)*T2-0.1753062e-2)*T2+.3989422793;
           Q0=T*(((((.32312e-5*T2-0.142078e-4)*T2+0.342468e-4)*T2-0.869791e-4)*T2+0.4564324e-3)*T2-0.0124669441);
           TA0=X-0.25*PI;
           BY0=2.0/sqrt(X)*(P0*sin(TA0)+Q0*cos(TA0));
           SH0=2.0/(PI*X)*S+BY0;
        }
        return SH0;
}

double STVH1(double X) {
/*      =============================================
!       Purpose: Compute Struve function H1(x)
!       Input :  x   --- Argument of H1(x) ( x ò 0 )
!       Output:  SH1 --- H1(x)
!       ============================================= */
        double A0,BY1,P1,PI,Q1,R,S,T,T2,TA1;
	int K, KM;
        double SH1;

        PI=3.141592653589793;
        R=1.0;
        if (X <= 20.0) {
           S=0.0;
           A0=-2.0/PI;
           for (K=1; K<=60; K++) {
              R=-R*X*X/(4.0*K*K-1.0);
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e15;
           }
e15:       SH1=A0*S;
        }
        else {
           S=1.0;
           KM=int(0.5*X);
           if (X > 50.0)  KM=25;
           for (K=1; K<=KM; K++) {
              R=-R*(4.0*K*K-1.0)/(X*X);
              S=S+R;
              if (fabs(R) < fabs(S)*1.0e-12) goto e25;
           }
e25:       T=4.0/X;
           T2=T*T;
           P1=((((0.42414e-5*T2-0.20092e-4)*T2+0.580759e-4)*T2-0.223203e-3)*T2+0.29218256e-2)*T2+0.3989422819;
           Q1=T*(((((-0.36594e-5*T2+0.1622e-4)*T2-0.398708e-4)*T2+0.1064741e-3)*T2-0.63904e-3)*T2+0.0374008364);
           TA1=X-0.75*PI;
           BY1=2.0/sqrt(X)*(P1*sin(TA1)+Q1*cos(TA1));
           SH1=2.0/PI*(1.0+S/(X*X))+BY1;
        }
        return SH1;
}

int integral_J(float *U_f,float *r, float *f, float *out,float *c, int nc, int nr, int nf,int nthreads){
	int indx;
	float kernel;
	float k,fl,cl,r1,r2,g1,g2,a,b;
	float dr0,B01,B02;
	int indx_d,B0n;
#ifdef useOMP
    omp_set_num_threads(nthreads);
#endif
#ifdef useOMP
#pragma omp parallel
#pragma omp for private(k,kernel,fl,cl,r1,r2,g1,g2,a,b,dr0,B01,B02,indx_d,B0n)
#endif
	for(int i =0;i<nf;i++)
	{
		for(int j =0;j<nc;j++)
		{
			indx = i + j*nf; 
			fl = f[i];
			cl = c[j];
			k = 2*M_PI*fl/cl;
			B02 = 0;
			kernel = 0;
			for (int ir=1;ir<nr;ir++){
				indx_d = i + ir*nf;
				g1 = U_f[indx_d-nf];
				g2 = U_f[indx_d];
				r1 = r[ir-1];
				r2 = r[ir];
				dr0 = fmaxf((r2-r1),0.1);
				a = g1-r1*(g2-g1)/dr0;
				b = (g2-g1)/dr0;
				kernel += a*(r2*bessj1(k*r2)-r1*bessj1(k*r1))/k;
				kernel += b*(r2*r2*bessj1(k*r2)-r1*r1*bessj1(k*r1))/k;
				kernel += b*(r2*bessj0(k*r2)-r1*bessj0(k*r1))/k/k;
				B02 = k*r2*bessj0(k*r2)+M_PI*k*r2*(bessj1(k*r2)*STVH0(k*r2)-bessj0(k*r2)*STVH1(k*r2))/2;
                B01 = k*r1*bessj0(k*r1)+M_PI*k*r1*(bessj1(k*r1)*STVH0(k*r1)-bessj0(k*r1)*STVH1(k*r1))/2;
				kernel += -b*(B02-B01)/k/k/k;
			}
			out[indx] = kernel;
		}
	}
    return 0;
}


int integral_Y(float *U_f,float *r, float *f, float *out,float *c, int nc, int nr, int nf,int nthreads){
	int indx;
	float kernel;
	float k,fl,cl,r1,r2,g1,g2,a,b;
	float dr0,B01,B02;
	int indx_d,B0n;
#ifdef useOMP
    omp_set_num_threads(nthreads);
#endif
#ifdef useOMP
#pragma omp parallel
#pragma omp for private(k,kernel,fl,cl,r1,r2,g1,g2,a,b,dr0,B01,B02,indx_d,B0n)
#endif
	for(int i =0;i<nf;i++)
	{
		for(int j =0;j<nc;j++)
		{
			indx = i + j*nf; 
			fl = f[i];
			cl = c[j];
			k = 2*M_PI*fl/cl;
			B02 = 0;
			kernel = 0;
			for (int ir=1;ir<nr;ir++){
				indx_d = i + ir*nf;
				g1 = U_f[indx_d-nf];
				g2 = U_f[indx_d];
				r1 = r[ir-1];
				r2 = r[ir];
				dr0 = fmaxf((r2-r1),0.1);
				a = g1-r1*(g2-g1)/dr0;
				b = (g2-g1)/dr0;
				kernel += a*(r2*bessj1(k*r2)-r1*bessj1(k*r1))/k;
				kernel += b*(r2*r2*bessj1(k*r2)-r1*r1*bessj1(k*r1))/k;
				kernel += b*(r2*bessj0(k*r2)-r1*bessj0(k*r1))/k/k;
				B02 = k*r2*bessy0(k*r2)+M_PI*k*r2*(bessy1(k*r2)*STVH0(k*r2)-bessy0(k*r2)*STVH1(k*r2))/2;
                B01 = k*r1*bessy0(k*r1)+M_PI*k*r1*(bessy1(k*r1)*STVH0(k*r1)-bessy0(k*r1)*STVH1(k*r1))/2;
				kernel += -b*(B02-B01)/k/k/k;
			}
			out[indx] = kernel;
		}
	}
    return 0;
}



int trap_J(float *U_f,float *r, float *f, float *out,float *c, int nc, int nr, int nf,int nthreads){
	int indx;
	float kernel;
	float k,fl,cl,r1,r2,g1,g2,a,b;
	float dr0,B01,B02;
	int indx_d,B0n;
#ifdef useOMP
    omp_set_num_threads(nthreads);
#endif
#ifdef useOMP
#pragma omp parallel
#pragma omp for private(k,kernel,fl,cl,r1,r2,g1,g2,a,b,dr0,B01,B02,indx_d,B0n)
#endif
	for(int i =0;i<nf;i++)
	{
		for(int j =0;j<nc;j++)
		{
			indx = i + j*nf; 
			fl = f[i];
			cl = c[j];
			k = 2*M_PI*fl/cl;
			B02 = 0;
			kernel = 0;
			for (int ir=1;ir<nr;ir++){
				indx_d = i + ir*nf;
				g1 = U_f[indx_d-nf];
				g2 = U_f[indx_d];
				r1 = r[ir-1];
				r2 = r[ir];
				dr0 = fmaxf((r2-r1),0.1);
			    kernel += (g1*bessj0(k*r1)*r1+g2*(bessj0(k*r2))*r2)*dr0/2;
			}
			out[indx] = kernel;
		}
	}
    return 0;
}

int trap_Y(float *U_f,float *r, float *f, float *out,float *c, int nc, int nr, int nf,int nthreads){
	int indx;
	float kernel;
	float k,fl,cl,r1,r2,g1,g2,a,b;
	float dr0,B01,B02;
	int indx_d,B0n;
#ifdef useOMP
    omp_set_num_threads(nthreads);
#endif
#ifdef useOMP
#pragma omp parallel
#pragma omp for private(k,kernel,fl,cl,r1,r2,g1,g2,a,b,dr0,B01,B02,indx_d,B0n)
#endif
	for(int i =0;i<nf;i++)
	{
		for(int j =0;j<nc;j++)
		{
			indx = i + j*nf; 
			fl = f[i];
			cl = c[j];
			k = 2*M_PI*fl/cl;
			B02 = 0;
			kernel = 0;
			for (int ir=1;ir<nr;ir++){
				indx_d = i + ir*nf;
				g1 = U_f[indx_d-nf];
				g2 = U_f[indx_d];
				r1 = r[ir-1];
				r2 = r[ir];
				dr0 = fmaxf((r2-r1),0.1);
			    kernel += (g1*bessy0(k*r1)*r1+g2*(bessy0(k*r2))*r2)*dr0/2;
			}
			out[indx] = kernel;
		}
	}
    return 0;
}

int FJ(float *u_f,float *r, float *f,float *out, float *c,int nc,int nr,int nf,int type,int nThread)
{
    if(type==1)
		integral_J(u_f,r,out,c,f,nr,nc,nf,nThread);
	if(type==0)
		trap_J(u_f,r,out,c,f,nr,nc,nf,nThread);
    return 0;
}


int FH(float *u_f,float *r, float *f,float *out, float *c,int nc,int nr,int nf,int type,int nThread)
{
    if(type==1)
		integral_Y(u_f,r,out,c,f,nr,nc,nf,nThread);
	if(type==0)
		trap_Y(u_f,r,out,c,f,nr,nc,nf,nThread);
    return 0;
}
