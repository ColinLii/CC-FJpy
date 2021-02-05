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

#ifndef CROSSCORR_H
#define CROSSCORR_H

typedef struct CC_data{
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
}CC_data;



int CrossCorrelation(CC_data data);


#endif
