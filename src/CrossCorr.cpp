//
//  CC-FJpy
//  A Python Package for seismic ambient noise cross-correlation and Frequency-Bessel Transform Method
//
//  GNU General Public License, Version 3, 29 June 2007
//
//  Copyright © 2021 Xiaofei Chen Research Group,
//  Department of Earth and Space Sciences,
//  Southern University of Science and Technology, China.
//
#ifdef useOMP
#include <omp.h>
#endif
#include <fftw3.h>
#include<cstdlib>
#include <math.h>
#include <algorithm>
#include "CrossCorr.h"
using namespace std;


float sign(float x)
{
    if(x>0)
        return 1.0;
    else if(x<0)
        return -1.0;
    else
        return 0.0;

}

int CrossCorrelation(CC_data d)
{   
#ifdef useOMP
    omp_set_num_threads(d.nThreads);
#endif
    int npairs = (d.nsta-1)*d.nsta/2;
    int * CCnumbers;
    CCnumbers = (int *)malloc(sizeof(int)*npairs);
    for(int i=0;i<npairs;i++)
        CCnumbers[i] = 0;
    bool * ifCC;
    ifCC = (bool *) malloc(sizeof(bool)*d.nsta);
    int nshifts;
    if(d.fftlen>=d.npts){
        d.fftlen = d.npts;
        nshifts = 1;
    }
    else{
        nshifts = 1 + floorf((d.npts-d.fftlen)/d.steplen);
    }
    float *datar, *datai;
    float ** in;
    float absv;
    fftwf_complex **out;
    fftwf_plan *p;
    in = (float **)malloc(sizeof(float *)*d.nsta);
    out = (fftwf_complex **)malloc(sizeof(fftwf_complex *)*d.nsta);
    p = (fftwf_plan *) malloc(sizeof(fftwf_plan)*d.nsta);
    for(int i =0;i<d.nsta;i++){
        //in[i] = (float *)malloc(sizeof(float)*fftlen);
        //out[i] = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex)*fftlen);
        in[i] = fftwf_alloc_real(d.fftlen);
        out[i] = fftwf_alloc_complex(d.fftlen);
        p[i] = fftwf_plan_dft_r2c_1d(d.fftlen,in[i],out[i],FFTW_FORWARD);
    }

    datar = (float *)malloc(d.nf*d.nsta*sizeof(float));
    datai = (float *)malloc(d.nf*d.nsta*sizeof(float));
    if(d.ifonebit){
#ifdef useOMP
#pragma omp parallel
#pragma omp for private(absv)
#endif
    for(int i=0;i<d.nsta;i++){
        for(int j =0;j<d.npts;j++){
            d.data[j+i*d.npts] = sign(d.data[j+i*d.npts]);
        }
    }
    }

    for(int k=0;k<nshifts;k++){
#ifdef useOMP
#pragma omp parallel
#pragma omp for private(absv)
#endif
        for(int i=0;i<d.nsta;i++){
            if(((k*d.steplen)>d.startend[i*2])&&((k*d.steplen+d.fftlen)<d.startend[i*2+1]))
            {
                ifCC[i] = 1;
            }
            else
            {
                ifCC[i] = 0;   
            }
            
            if (ifCC[i]){
                for(int j=0;j<d.fftlen;j++){
                    in[i][j] = d.data[j+k*d.steplen+d.npts*i];
                }
                fftwf_execute(p[i]);
                if(d.ifspecwhitenning){
                    for (int j=0;j<d.nf;j++){
                        absv = sqrtf(out[i][j*d.fstride][0]*out[i][j*d.fstride][0]+out[i][j*d.fstride][1]*out[i][j*d.fstride][1]);
                        absv = max(absv,1e-8f);
                        if(absv!=0){
                            datar[j+d.nf*i] = out[i][j*d.fstride][0]/absv;
                            datai[j+d.nf*i] = out[i][j*d.fstride][1]/absv;
                        }
                    }
                }
                else{
                    for(int j=0;j<d.nf;j++){
                        datar[j+d.nf*i] = out[i][j*d.fstride][0];
                        datai[j+d.nf*i] = out[i][j*d.fstride][1];
                    }
                }
            }
        }
#ifdef useOMP
#pragma omp parallel
#pragma omp for
#endif
        for(int i=0;i<npairs; i++){
            if(ifCC[d.Pairs[i*2]]&ifCC[d.Pairs[i*2+1]]){
                for(int j=0;j<d.nf;j++){
                    d.ncfsr[j+i*d.nf] += datar[j+d.Pairs[i*2]*d.nf]*datar[j+d.Pairs[i*2+1]*d.nf] + datai[j+d.Pairs[i*2]*d.nf]*datai[j+d.Pairs[i*2+1]*d.nf];
                    d.ncfsi[j+i*d.nf] += datar[j+d.Pairs[i*2]*d.nf]*datai[j+d.Pairs[i*2+1]*d.nf] - datai[j+d.Pairs[i*2]*d.nf]*datar[j+d.Pairs[i*2+1]*d.nf];
                }
                CCnumbers[i] = CCnumbers[i] + 1;
            }
        }
    }
#ifdef useOMP
#pragma omp parallel
#pragma omp for
#endif
    for (int i=0;i<npairs;i++){
        if(CCnumbers[i]>0){
            for(int j =0;j<d.nf;j++){
                d.ncfsr[j+i*d.nf] = d.ncfsr[j+i*d.nf]/CCnumbers[i];
                d.ncfsi[j+i*d.nf] = d.ncfsi[j+i*d.nf]/CCnumbers[i];
            }
        }
    }
    for (int i=0;i<d.nsta;i++){
        fftwf_destroy_plan(p[i]);
        fftwf_free(out[i]);
        fftwf_free(in[i]);
    }
    free(in);
    free(out);
    free(p);
    free(datar);
    free(datai);
    return 0;
}
