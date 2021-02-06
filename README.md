# CC-FJpy: A Python Package for seismic ambient noise cross-correlation and the frequency-Bessel transform method

## Copyright

Xiaofei Chen Research Group

Department of Earth and Space Sciences, SUSTech, China.

## Appointed DOI for preprint: 10.1002/essoar.10506115.1

## Installation

Python 3 is required.

Anaconda environment is required for installation by `make install`.

**Make sure that the CUDA and Python in the same version for you install and run.**



### Before Installation: Download & Compile `fftw` Library

```
make fftw
```

### Installation for GPU version (default)

```
make
make install
```

**If the `nvcc` compiler or `$cudahome` cannot be found, `ccfj`will be compiled by default in the CPU version**

### Installation for CPU version

```
make cpu
make install
```


## Usage

```
import ccfj
```

### Cross-correlation (CC)


```
CCFs=ccfj.CC(
    npts,nsta,nf,fftlen,
    Pairs,startend,data,
    overlaprate=0.0,
    nThreads=8,
    fstride=1,
    ifonebit=0,
    ifspecwhittenning=1)
```

- `npts`: The number of points for data of one station read in
- `nsta`: The number of stations
- `nf`: The number of points of the frequency domain CC output
- `fftlen`: The number of points for one CC
- `Pairs`: A numpy array for Station Pairs, the `dtype` should be `np.int32`

You can use function GetStationPairs to generate Pairs.

startend: A numpy array records the start point and end point for each station, the dtype should be np.int32. For example, if the npts length is for one day, and station A have all of A day's data. Then for station A the startend is [0,npts]

- `data`: The seismic records, the `dtype` should be `np.float32`
- `overlaprate`: The rate for overlap, do not >= 1

- `nThreads`: The number of Threads for reading data and CC

- `fstride`: The output frequency stride

- `ifonebit`: if perform onebit

- `ifspecwhittening`: if perform specwhittenning

- `CCFs`: The output is the noise cross-correlation functions (CCFs)

According to our experience, `ifonebit` and `ifspecwhittening`, pick one of them is ok.


The specific example, please read example_CC.ipynb.

### frequency-Bessel transform method (F-J method)
#### For ambient noise 

```
out = ccfj.fj_noise(uf,r,c,f,fstride=1,itype=1,func=0)
```

The F-J Method for CCFs

- `uf`: the CCFs in the frequency domain

- `r`: the list of station distances (unit: m)

- `c`: the list of phase velocities you want to calculate

- `f`: the list of frequencies. The number of points of f should be consistent with the columns of `uf`

- `fstride`: stride of frequency for output

- `itype`: `0` for trapezoidal integral; `1` for linear approximate

- `func`: `0` for Bessel function; `1` for Hankel function

- `out`: the output dispersion spectrum


#### For earthquakes

```
out = ccfj.fj_earthquake(u,r,c,f,fstride=1,itype=1,func=0)
```


- `u`: the records in time domain

- `r`: the list of station distances (unit: m)

- `c`: the list of phase velocities you want to calculate

- `f`: the list of frequencies. The number of points of f should be consistent with the columns of `uf`

- `fstride`: stride of frequency for output

- `itype`: `0` for trapezoidal integral; `1` for linear approximate

- `func`: `0` for Bessel function; `1` for Hankel function

- `out`: the output dispersion spectrum

### Mutli-windows F-J method (MWFJ)

This is mainly for earthquake

```
out = ccfj.MWFJ(u,r,c,f,Fs,nwin,winl,winr,taper=0.9,fstride=1,itype=0,func=0)
```

- `u`: the records in time domain

- `r`: the list of station distances (unit: m)

- `c`: the list of phase velocities you want to calculate

- `f`: the list of frequencies. The number of points of f should be consistent with the columns of `uf`

- `Fs`: The sample frequency

- `nwin`: number of time windows

- `winl`: list of left side of time windows

- `winr`: list of right side of time windows

- `fstride`: stride of frequency for output

- `itype`: `0` for trapezoidal integral; `1` for linear approximate

- `func`: `0` for Bessel function; `1` for Hankel function

- `out`: the output dispersion spectrum


## Uninstall

```
make clean
```

