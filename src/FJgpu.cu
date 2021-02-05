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
#include "cuda_helper.h"
#include "FJgpu.hh"

using namespace std;

__device__ double bessj0( double x )
/*------------------------------------------------------------*/
/* PURPOSE: Evaluate Bessel fun//ction of first kind and order  */
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

__device__ double bessy0( double x )
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

__device__ double STVH0(double X) {
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

__device__ double STVH1(double X) {
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

__global__ void trap_J(float *U_f, float *r, float *out, float *c, float *f, int nr, int nc, int nf) {
	bool validr = true;
	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	int indx = 0;
	indx = gtidx + gtidy * nf;
	float kernel;
	float k, fl, cl, r1, r2, g1, g2;
	float dr0;
	int indx_d;
	if ((gtidx >= nf) || (gtidy >= nc))
		validr = false;
	if (validr) {
		fl = f[gtidx];
		cl = c[gtidy];
		k = 2 * M_PI *fl / cl;
		kernel = 0.0;
		for (int i = 1; i < nr; i++) {
			indx_d = gtidx + i * nf;
			//g1 = U_f[indx_d-1];
			g1 = U_f[indx_d - nf];
			g2 = U_f[indx_d];
			r1 = r[i - 1];
			r2 = r[i];
			dr0 = fmaxf((r2 - r1), 0.1);
			kernel += (g1*j0(k*r1)*r1+g2*j0(k*r2)*r2)*dr0/2;
		}
		out[indx] = kernel;
	}
}

__global__ void trap_Y(float *U_f, float *r, float *out, float *c, float *f, int nr, int nc, int nf) {
	bool validr = true;
	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	int indx = 0;
	indx = gtidx + gtidy * nf;
	float kernel;
	float k, fl, cl, r1, r2, g1, g2;
	float dr0;
	int indx_d;
	if ((gtidx >= nf) || (gtidy >= nc))
		validr = false;
	if (validr) {
		fl = f[gtidx];
		cl = c[gtidy];
		k = 2 * M_PI *fl / cl;
		kernel = 0.0;
		for (int i = 1; i < nr; i++) {
			indx_d = gtidx + i * nf;
			//g1 = U_f[indx_d-1];
			g1 = U_f[indx_d - nf];
			g2 = U_f[indx_d];
			r1 = r[i - 1];
			r2 = r[i];
			dr0 = fmaxf((r2 - r1), 0.1);
			kernel += (g1*bessy0(k*r1)*r1+g2*(bessy0(k*r2))*r2)*dr0/2;
		}
		out[indx] = kernel;
	}
}


__global__ void integral_J(float *U_f, float *r, float *out, float *c, float *f, int nr, int nc, int nf) {
	bool validr = true;
	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	int indx = 0;
	indx = gtidx + gtidy * nf;
	float kernel;
	float k, fl, cl, r1, r2, g1, g2, a, b;
	float dr0, B01, B02;
	int indx_d;
	if ((gtidx >= nf) || (gtidy >= nc))
		validr = false;
	if (validr) {
		fl = f[gtidx];
		cl = c[gtidy];
		k = 2 * M_PI *fl / cl;
		kernel = 0.0;
		for (int i = 1; i < nr; i++) {
			indx_d = gtidx + i * nf;
			//g1 = U_f[indx_d-1];
			g1 = U_f[indx_d - nf];
			g2 = U_f[indx_d];
			r1 = r[i - 1];
			r2 = r[i];
			dr0 = fmaxf((r2 - r1), 0.1);
			a = g1 - r1 * (g2 - g1) / dr0;
			b = (g2 - g1) / dr0;
			kernel += a * (r2*j1(k*r2) - r1 * j1(k*r1)) / k;
			kernel += b * (r2*r2*j1(k*r2) - r1 * r1*j1(k*r1)) / k;
			kernel += b * (r2*j0(k*r2) - r1 * j0(k*r1)) / k / k;
			/* B01 = 0.0; */
			/* B0n = floorf(floorf(B01 / 2 + 10)*1.5); */
			/* for (int j = 0; j < B0n; j++) */
			/* 	B01 += 2*jn(j * 2 + 1, k*r1); */
			/* B02 = 0.0; */
			/* B0n = floorf(floorf(B02 / 2 + 10)*1.5); */
			/* for (int j = 0; j < B0n; j++) */
			/* 	B02 += 2*jn(j * 2 + 1, k*r2); */
                        B02 = k*r2*j0(k*r2)+M_PI*k*r2*(j1(k*r2)*STVH0(k*r2)-j0(k*r2)*STVH1(k*r2))/2;
                        B01 = k*r1*j0(k*r1)+M_PI*k*r1*(j1(k*r1)*STVH0(k*r1)-j0(k*r1)*STVH1(k*r1))/2;
			kernel += -b * (B02 - B01 ) /k / k / k;
		}
		out[indx] = kernel;
	}
}

__global__ void integral_Y(float *U_f, float *r, float *out, float *c, float *f, int nr, int nc, int nf) {
	bool validr = true;
	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	int indx = 0;
	indx = gtidx + gtidy * nf;
	float kernel;
	float k, fl, cl, r1, r2, g1, g2, a, b;
	float dr0, B01, B02;
	int indx_d;
	if ((gtidx >= nf) || (gtidy >= nc))
		validr = false;
	if (validr) {
		fl = f[gtidx];
		cl = c[gtidy];
		k = 2 * M_PI *fl / cl;
		kernel = 0.0;
		for (int i = 1; i < nr; i++) {
			indx_d = gtidx + i * nf;
			//g1 = U_f[indx_d-1];
			g1 = U_f[indx_d - nf];
			g2 = U_f[indx_d];
			r1 = r[i - 1];
			r2 = r[i];
			dr0 = fmaxf((r2 - r1), 0.1);
			a = g1 - r1 * (g2 - g1) / dr0;
			b = (g2 - g1) / dr0;
			kernel += a * (r2*y1(k*r2) - r1 * y1(k*r1)) / k;
			kernel += b * (r2*r2*y1(k*r2) - r1 * r1*y1(k*r1)) / k;
			kernel += b * (r2*y0(k*r2) - r1 * y0(k*r1)) / k / k;
			/* B01 = 0.0; */
			/* B0n = floorf(floorf(B01 / 2 + 10)*1.5); */
			/* for (int j = 0; j < B0n; j++) */
			/* 	B01 += 2*yn(j * 2 + 1, k*r1); */
			/* B02 = 0.0; */
			/* B0n = floorf(floorf(B02 / 2 + 10)*1.5); */
			/* for (int j = 0; j < B0n; j++) */
			/* 	B02 += 2*yn(j * 2 + 1, k*r2); */
                        B02 = k*r2*y0(k*r2)+M_PI*k*r2*(y1(k*r2)*STVH0(k*r2)-y0(k*r2)*STVH1(k*r2))/2;
                        B01 = k*r1*y0(k*r1)+M_PI*k*r1*(y1(k*r1)*STVH0(k*r1)-y0(k*r1)*STVH1(k*r1))/2;
			kernel += -b * (B02 - B01 ) /k / k / k;
		}
		out[indx] = kernel;
	}
}

__global__ void integral_c(float *U_f,float *r,float *out,float *c,float *f,int nr,int nc,int nf){
    bool validr = true;
	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	int indx = 0;
	indx = gtidx + gtidy * nf;
	float kernel;
	float k,fl,cl,r1,r2,g1,g2,a,b;
	float dr0,B01,B02;
	int indx_d,B0n;
	if ((gtidx >=nf) ||(gtidy>=nc))
		validr = false;
	if(validr){
		fl = f[gtidx];
		cl = c[gtidy];
		k = 2 * M_PI *fl/cl;
		kernel = 0.0;
		for(int i=1;i<nr;i++){
			indx_d = gtidx + i*nf;
			//g1 = U_f[indx_d-1];
			g1 = U_f[indx_d-nf];
			g2 = U_f[indx_d];
			r1 = r[i-1];
			r2 = r[i];
			dr0 = fmaxf((r2 - r1),0.1);
			a = g1 - r1*(g2-g1)/dr0;
			b = (g2-g1)/dr0;
			kernel += a*(r2*j1f(k*r2)-r1*j1f(k*r1))/k;
			kernel += b*(r2*r2*j1f(k*r2)-r1*r1*j1f(k*r1))/k;
			kernel += b*(r2*j0f(k*r2)-r1*j0f(k*r1))/k/k;
			B01 = 0.0;
			B0n = floorf(floorf(B01/2+10)*1.5);
			for(int j=0;j<B0n;j++)
				B01 += jnf(j*2+1,k*r1);
			B02 = 0.0;
			B0n = floorf(floorf(B02/2+10)*1.5);
			for(int j=0;j<B0n;j++)
				B02 += jnf(j*2+1,k*r2);
			//kernel += -b*(B02-B01+j0f(k*r1)*0.5+j0f(k*r2)*0.5)*dr/k/k;
			kernel += -b*(B02-B01)/k/k/k;
		}
		out[indx] = kernel;
	}
}

__global__ void trap_c(float *U_f,float *r,float *out,float *c,float *f,int nr,int nc,int nf){
	bool validr = true;
	const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
	const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
	int indx = 0;
	indx = gtidx + gtidy * nf;
	float kernel;
	float k,fl,cl,r1,r2,g1,g2,dr0;
	int indx_d;
	if ((gtidx >=nf) ||(gtidy>=nc))
		validr = false;
	if(validr){
		fl = f[gtidx];
		cl = c[gtidy];
		k = 2 * M_PI *fl/cl;
		kernel = 0.0;
		for(int i=1;i<nr;i++){
			indx_d = gtidx + i*nf;
			g1 = U_f[indx_d-nf];
			g2 = U_f[indx_d];
			r1 = r[i-1];
			r2 = r[i];
			dr0 = fmaxf((r2 - r1),0.1);
			kernel += 0.5*(j0f(k*r1)*g1*r1+j0f(k*r2)*g2*r2)*dr0;
		}
		out[indx] = kernel;
	}
}


int FJ(float *u_f,float *r, float *f,float *out, float *c,int nc,int nr,int nf,int type)
{	
    float *u_fc,*rc,*fc,*cc,*outc;
    dim3 block,grid;
    block.x = 32;
    block.y = 16;
	grid.x = (unsigned int)ceil((float)nf/block.x);
	grid.y = (unsigned int)ceil((float)nc/block.y);
    checkCudaErrors(cudaMalloc((void **)&u_fc,nr*nf*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&rc,nr*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&fc,nf*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&cc,nc*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&outc,nf*nc*sizeof(float)));

    checkCudaErrors(cudaMemcpy(u_fc,u_f,nr*nf*sizeof(float),cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(fc,f,nf*sizeof(float),cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(rc,r,nr*sizeof(float),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(cc,c,nc*sizeof(float),cudaMemcpyHostToDevice));
	if(type==1)
		integral_c<<<grid,block>>>(u_fc,rc,outc,cc,fc,nr,nc,nf);
	if(type==0)
		trap_c<<<grid,block>>>(u_fc,rc,outc,cc,fc,nr,nc,nf);
    checkCudaErrors(cudaMemcpy(out, outc, nf*nc*sizeof(float), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(u_fc));
    checkCudaErrors(cudaFree(rc));
    checkCudaErrors(cudaFree(fc));
    checkCudaErrors(cudaFree(cc));
    checkCudaErrors(cudaFree(outc));
    return 0;
}


int FHr(float *u_f,float *r, float *f,float *out, float *c,int nc,int nr,int nf,int type)
{
    float *u_fc,*rc,*fc,*cc,*outc;
    dim3 block,grid;
    block.x = 32;
    block.y = 16;
	grid.x = (unsigned int)ceil((float)nf/block.x);
	grid.y = (unsigned int)ceil((float)nc/block.y);
    checkCudaErrors(cudaMalloc((void **)&u_fc,nr*nf*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&rc,nr*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&fc,nf*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&cc,nc*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&outc,nf*nc*sizeof(float)));

    checkCudaErrors(cudaMemcpy(u_fc,u_f,nr*nf*sizeof(float),cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(fc,f,nf*sizeof(float),cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(rc,r,nr*sizeof(float),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(cc,c,nc*sizeof(float),cudaMemcpyHostToDevice));
	if(type==1)
		integral_J<<<grid,block>>>(u_fc,rc,outc,cc,fc,nr,nc,nf);
	if(type==0)
		trap_J<<<grid,block>>>(u_fc,rc,outc,cc,fc,nr,nc,nf);
    checkCudaErrors(cudaMemcpy(out, outc, nf*nc*sizeof(float), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(u_fc));
    checkCudaErrors(cudaFree(rc));
    checkCudaErrors(cudaFree(fc));
    checkCudaErrors(cudaFree(cc));
    checkCudaErrors(cudaFree(outc));
    return 0;
}


int FHi(float *u_f,float *r, float *f,float *out, float *c,int nc,int nr,int nf,int type)
{
    float *u_fc,*rc,*fc,*cc,*outc;
    dim3 block,grid;
    block.x = 32;
    block.y = 16;
	grid.x = (unsigned int)ceil((float)nf/block.x);
	grid.y = (unsigned int)ceil((float)nc/block.y);
    checkCudaErrors(cudaMalloc((void **)&u_fc,nr*nf*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&rc,nr*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&fc,nf*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&cc,nc*sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&outc,nf*nc*sizeof(float)));

    checkCudaErrors(cudaMemcpy(u_fc,u_f,nr*nf*sizeof(float),cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(fc,f,nf*sizeof(float),cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(rc,r,nr*sizeof(float),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(cc,c,nc*sizeof(float),cudaMemcpyHostToDevice));
	if(type==1)
		integral_Y<<<grid,block>>>(u_fc,rc,outc,cc,fc,nr,nc,nf);
	if(type==0)
		trap_Y<<<grid,block>>>(u_fc,rc,outc,cc,fc,nr,nc,nf);
    checkCudaErrors(cudaMemcpy(out, outc, nf*nc*sizeof(float), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(u_fc));
    checkCudaErrors(cudaFree(rc));
    checkCudaErrors(cudaFree(fc));
    checkCudaErrors(cudaFree(cc));
    checkCudaErrors(cudaFree(outc));
    return 0;
}
