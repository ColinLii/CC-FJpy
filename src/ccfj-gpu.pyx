"""
CC-FJpy: A Python Package for seismic ambient noise cross-correlation and the frequency-Bessel transform method.

:copyright:
 Xiaofei Chen Research Group, Department of Earth and Space Sciences, SUSTech, China.
:license:
 GNU Lesser General Public License, Version 3
 (https://www.gnu.org/copyleft/lesser.html)
"""
from cython.parallel import prange, parallel, threadid
from libc.stdio cimport printf
from numpy cimport ndarray as array
from scipy import fftpack
import numpy as np
import os
cimport numpy as np
cimport cython


cdef extern from "CrossCorr.h":
    struct CC_data:    
        int nsta;
        int npts;
        int nf;
        int fstride;
        int fftlen;
        int steplen;
        int *Pairs;
        int *startend;
        float *data;
        float *ncfsr;
        float *ncfsi;
        int ifonebit;
        int ifspecwhitenning;
        int nThreads;
    int CrossCorrelation(CC_data data)

cdef extern from "FJgpu.hh":
    int FJ(float *u_f,float *r, float *f,float *out, float *c,int nc,int nr,int nf, int type)
    int FHr(float *u_f,float *r, float *f,float *out, float *c,int nc,int nr,int nf,int type)
    int FHi(float *u_f,float *r, float *f,float *out, float *c,int nc,int nr,int nf,int type)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef CC_full(
    int npts, int nsta, int nf, int fstride,
    int fftlen,int steplen,int ifonebit,int ifspecwhitenning,
    int nThreads,
    int[:] Pairs, int[:] startend,
    float[:] data,
    float[:] ncfsr, float[:] ncfsi):

    cdef CC_data d
    d.npts = npts
    d.nsta = nsta
    d.nf = nf
    d.fstride = fstride
    d.fftlen = fftlen
    d.steplen = steplen
    d.ifonebit = ifonebit
    d.ifspecwhitenning = ifspecwhitenning
    d.nThreads = nThreads
    d.Pairs = <int *>& Pairs[0]
    d.startend = <int *>&startend[0]
    d.data = <float *> & data[0]
    d.ncfsr = <float *> &ncfsr[0]
    d.ncfsi = <float *> &ncfsi[0]
    CrossCorrelation(d)
    
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef fj_full(
    float[:] U_f, 
    float[:] r, 
    float[:] f, 
    float[:] out, 
    float[:] c, 
    int nc, 
    int nr, 
    int nf,
    int type):

    U_f1 = <float *> & U_f[0]
    r1 = <float *> & r[0]
    f1 = <float *> & f[0]
    out1 = <float *> & out[0]
    c1 = <float *> & c[0]
    
    FJ(U_f1, r1, f1, out1, c1, nc, nr, nf, type)


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef fHr_full(
    float[:] U_f, 
    float[:] r, 
    float[:] f, 
    float[:] out, 
    float[:] c, 
    int nc, 
    int nr, 
    int nf,
    int type):

    U_f1 = <float *> & U_f[0]
    r1 = <float *> & r[0]
    f1 = <float *> & f[0]
    out1 = <float *> & out[0]
    c1 = <float *> & c[0]
    
    FHr(U_f1, r1, f1, out1, c1, nc, nr, nf, type)
    
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef fHi_full(
    float[:] U_f, 
    float[:] r, 
    float[:] f, 
    float[:] out, 
    float[:] c, 
    int nc, 
    int nr, 
    int nf,
    int type):

    U_f1 = <float *> & U_f[0]
    r1 = <float *> & r[0]
    f1 = <float *> & f[0]
    out1 = <float *> & out[0]
    c1 = <float *> & c[0]
    
    FHi(U_f1, r1, f1, out1, c1, nc, nr, nf, type)

def GetStationPairs(nsta):
    StationPair = []
    for ii in range(nsta):
        for jj in range(ii+1,nsta):
            StationPair.append(ii)
            StationPair.append(jj)
    StationPair = np.array(StationPair,dtype=np.int32)
    return StationPair
    
def CC(
    npts,nsta,nf,fftlen,
    Pairs,startend,data,
    overlaprate=0.0,
    nThreads=8,
    fstride=1,
    ifonebit=0,
    ifspecwhittenning=1):
    steplen = int(fftlen *(1-overlaprate))
    nPairs = int(len(Pairs)/2)
    
    if data.dtype != np.float32:
        data = np.array(data,dtype=np.float32) # Make sure the data type
    
    ncfsr = np.zeros(nf*nPairs,dtype=np.float32)
    ncfsi = np.zeros(nf*nPairs,dtype=np.float32)
    CC_full(npts,nsta,nf,fstride,fftlen,steplen,ifonebit,ifspecwhittenning,nThreads,Pairs,startend,data,ncfsr,ncfsi)
    ncfsr = ncfsr.reshape([nPairs,nf])
    ncfsi = ncfsi.reshape([nPairs,nf])
    ncfs = ncfsr + ncfsi*1j
    return ncfs
    
def win(npts,Fs,T1,T2,taper=0.8):
    n1 = max(0,np.floor(T1*Fs))
    n2 = min(npts,np.floor(T2*Fs))
    out = np.ones(npts)
    delta = np.floor((n2-n1)*(1-taper)/2)
    out[0:int(n1)]=0
    out[int(n2):int(npts)] = 0
    for i in range(int(delta)):
        out[int(n1+i)] = np.sin(np.pi/2*i/delta)
        out[int(n2-i-1)] = np.sin(np.pi/2*i/delta)
    return out


def fj(uf,r,c,f,fstride=1,itype=1):
    nr = len(r)
    nc = len(c)
    nf = len(f)
    uf = uf[:,::fstride]
    uf = uf.reshape(nr*nf)

    if r.dtype != np.float32:
        r = np.array(r,dtype=np.float32) # Make sure the data type
    if f.dtype != np.float32:
        f = np.array(f,dtype=np.float32) # Make sure the data type
    if uf.dtype !=np.float32:
        uf = np.array(uf,dtype=np.float32) # Make sure the data type
    if c.dtype != np.float32:
        c = np.array(c,dtype=np.float32)
    
    out = np.zeros(nf*nc,dtype=np.float32)
    fj_full(uf,r,f,out,c,nc,nr,nf,itype)
    out = out.reshape([nc,nf])
    return out

def fhr(uf,r,c,f,fstride=1,itype=1):
    nr = len(r)
    nc = len(c)
    nf = len(f)
    uf = uf[:,::fstride]
    uf = uf.reshape(nr*nf)
    if r.dtype != np.float32:
        r = np.array(r,dtype=np.float32) # Make sure the data type
    if f.dtype != np.float32:
        f = np.array(f,dtype=np.float32) # Make sure the data type
    if c.dtype != np.float32:
        c = np.array(c,dtype=np.float32)
    if uf.dtype !=np.float32:
        uf = np.array(uf,dtype=np.float32) # Make sure the data type
    out = np.zeros(nf*nc,dtype=np.float32)
    fHr_full(uf,r,f,out,c,nc,nr,nf,itype)
    out = out.reshape([nc,nf])
    return out
    
def fhi(uf,r,c,f,fstride=1,itype=1):
    nr = len(r)
    nc = len(c)
    nf = len(f)
    uf = uf[:,::fstride]
    uf = uf.reshape(nr*nf)
    if r.dtype != np.float32:
        r = np.array(r,dtype=np.float32) # Make sure the data type
    if f.dtype != np.float32:
        f = np.array(f,dtype=np.float32) # Make sure the data type
    if c.dtype != np.float32:
        c = np.array(c,dtype=np.float32)
    if uf.dtype !=np.float32:
        uf = np.array(uf,dtype=np.float32) # Make sure the data type

    out = np.zeros(nf*nc,dtype=np.float32)
    fHi_full(uf,r,f,out,c,nc,nr,nf,itype)
    out = out.reshape([nc,nf])
    return out

def fj_noise(uf,r,c,f,
    fstride = 1,
    itype = 1, # 0 for trapz integral, 1 for linear approximate
    func = 0, # 'B' for Bessel funciton, 'H' for Hankel function
    num=-1
    ):
    if num != -1:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(int(num))
        
    indx = np.argsort(r)
    r = r[indx]
    uf = uf[indx]
    nr = len(r)
    if func == 0:
        out = fj(uf,r,c,f,fstride,itype)
    elif func == 1:
        outr = fhr(uf,r,c,f,fstride,itype)
        ufi = np.zeros(np.shape(uf),dtype=np.float32)
        for i in range(nr):
            ufi[i,:] = fftpack.hilbert(uf[i,:])
        outi = fhi(ufi,r,c,f,fstride,itype)
        out = outr - outi
    else:
        print('set func as 0 for Bessel function 1 for Hankel function')
        return 0
    for i in range(len(f)):
        out[:,i] = out[:,i]/max(np.abs(out[:,i]))
    return out

def fj_earthquake(u,r,c,f,fstride=1,itype=1,func=0,num=-1):
    if num != -1:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(int(num))
        
    indx = np.argsort(r)
    r = r[indx]
    u1 = u[indx]
    uf = np.fft.rfft(u1)
    uf = uf[:,0:len(f)*fstride:fstride]
    ufr = np.real(uf)
    ufi = np.imag(uf)
    if func == 0:
        outr = fj(ufr,r,c,f,fstride,itype)
        outi = fj(ufi,r,c,f,fstride,itype)
        out = np.sqrt(outr**2+outi**2)
    elif func == 1:
        outr = fhr(ufr,r,c,f,fstride,itype) - fhi(ufi,r,c,f,fstride,itype)
        outi = fhr(ufi,r,c,f,fstride,itype) + fhi(ufr,r,c,f,fstride,itype)
        out = np.sqrt(outr**2+outi**2)
        #out = outi
    else:
        print('set func as 0 for Bessel function 1 for Hankel function')
        return 0
    for i in range(len(f)):
        out[:,i] = out[:,i]/max(np.abs(out[:,i]))
    return out

def MWFJ(u,r,c,f,Fs,nwin,winl,winr,
    taper =0.9,
    fstride= 1,
    itype = 1,
    func = 0,
    num = -1):
        
    nr = len(r)
    npts = len(u[0])
    out = np.zeros([nwin,len(c),len(f)])
    for i in range(0,nwin):
        u0 = np.zeros([nr,npts])
        for j  in range(nr):
            tmp = win(npts,Fs,winl[i,j],winr[i,j],taper)
            u0[j,:] = u[j,:]*tmp
        out[i,:,:] = fj_earthquake(u0,r,c,f,fstride=fstride,itype=itype,func = func,num=num)
    return out
