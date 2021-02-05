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
#ifndef _FJCPU_
#define _FJCPU_
#ifdef windows
#define EXPORTFUNC extern "C" __declspec(dllexport)
#else
#define EXPORTFUNC extern "C"
#endif

EXPORTFUNC int FJ(float *U_f,float *r, float *f, float *out,float *c, int nc, int nr, int nf,int type,int nThread);
EXPORTFUNC int FH(float *U_f,float *r, float *f, float *out,float *c, int nc, int nr, int nf,int type, int nThread);
#endif
