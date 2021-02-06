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
#define EXPORTFUNC extern "C"

EXPORTFUNC int FJ(float *u_f,float *r, float *f,float *out, float *c,int nc,int nr,int nf,int type);
EXPORTFUNC int FHr(float *u_f,float *r, float *f,float *out, float *c,int nc,int nr,int nf,int type);
EXPORTFUNC int FHi(float *u_f,float *r, float *f,float *out, float *c,int nc,int nr,int nf,int type);
