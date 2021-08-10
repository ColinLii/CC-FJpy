# CC-FJpy: A Python Package for seismic ambient noise cross-correlation (CC) and the frequency-Bessel transform (FJ)  method

# CC-FJpy: 背景噪声互相关和频率贝塞尔变换法软件包

## Copyright

Xiaofei Chen Research Group

Department of Earth and Space Sciences, SUSTech, China.

## Related article https://doi.org/10.1785/0220210042

## 相关论文 https://doi.org/10.1785/0220210042

## 版权所有

陈晓非课题组
南方科技大学，地球与空间科学系

## Installation
## 安装

Python 3 is required.

Anaconda environment is required for installation by `make install`.

目前仅支持 Python 3，推荐使用Anaconda配置Python环境。请确保您所安装和使用的Python以及编译器是同一个。

**Make sure that the CUDA and Python in the same version for you install and run.**



Before Installation: Download & Compile `fftw` Library.

在安装前请安装`fftw`库(www.fftw.org) ,可以在当前目录通过如下命令一键安装。

```
make fftw
```

After the installation of the fftw, install the CC-FJpy with the following commands

在安装fftw库之后可以在当前目录通过如下两行命令完成CC-FJpy的安装。
```
make
make install
```

**If the machine has an NVIDIA graphics card and the `nvcc` compiler is already in the path, the `GPU` version of CC-FJpy will be installed. If the above conditions are not met, the `CPU` version will be installed.** If you want to install only the `CPU` version, you can install it with the following commands. 

**若本机有英伟达显卡且`nvcc`编译器已经在路径中时，`GPU`版本的CC-FJpy将会被安装，如果不满足上述条件则会安装`CPU`版本。** 若想仅安装`CPU`版本，则可以通过如下两行命令安装。


```
make cpu
make install
```

## Usage
## 使用
Please import the package before use its functions.

请在使用前先导入。
```
import ccfj
```

### Cross-correlation (CC)
### 互相关计算


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

- `npts`: The number of points for data of one station read in（每次互相关每个台所使用数据点数，在程序中对于每个台是一致的，没有数据的点可以通过补0处理
）
- `nsta`: The number of stations（台站总数）
- `nf`: The number of points of the frequency domain CC output （输出互相关的频率点数）
- `fftlen`: The number of points for one CC （单词互相关所用点数）
- `Pairs`: A numpy array for Station Pairs, the `dtype` should be `np.int32` （要计算的台站对。一维数组,长度为nparis*2。要计算的第i个台站对互相关为`Pairs[i*2,i*2+1]`）
You can use function GetStationPairs to generate Pairs.

- `startend`: A numpy array records the start point and end point for each station, the dtype should be np.int32. For example, if the npts length is for one day, and station A have all of A day's data. Then for station A the startend is [0,npts] （记录每个台的记录开始时间和结束时间的数组。）

- `data`: The seismic records, the `dtype` should be `np.float32` （连续地震记录，长度为nsta*npts的一维数组。将每个台的记录（不满npts的点补0）依次排在该数组中。）
- `overlaprate`: The rate for overlap, do not >= 1

- `nThreads`: The number of Threads for reading data and CC （互相关计算中omp所使用的核数）

- `fstride`: The output frequency stride （输出的互相关对应的频率点数相对于原始数据傅里叶变换（fftlen）时频谱点的间隔）

- `ifonebit`: if perform onebit （是否使用onebit）

- `ifspecwhittening`: if perform specwhittenning （是否使用谱白化）

- `CCFs`: The output is the noise cross-correlation functions (CCFs) （输出为频率域互相关函数）

According to our experience, `ifonebit` and `ifspecwhittening`, pick one of them is ok. 

根据我们的经验，`ifonebit` 和`ifspecwhittening`二者选其一即可。


The specific example, please read example_CC.ipynb.

具体例子请参考 example_CC.ipynb。

### frequency-Bessel transform method (F-J method)
### 频率-贝塞尔变换法
#### For ambient noise 
#### 对背景噪声数据

```
out = ccfj.fj_noise(uf,r,c,f,fstride=1,itype=1,func=0,num=20)
```

The F-J Method for CCFs

- `uf`: the CCFs in the frequency domain （频率域互相关函数,二维数组，大小为`[互相关数，nf]`）

- `r`: the list of station distances (unit: m) （对应的台间距数组，单位m，长度为互相关数）

- `c`: the list of phase velocities you want to calculate （对要计算的相速度数组，单位m/s）

- `f`: the list of frequencies. The number of points of f should be consistent with the columns of `uf` (频率域互相关对应的频率)

- `fstride`: stride of frequency for output （默认1即可）

- `itype`: `0` for trapezoidal integral; `1` for linear approximate （积分类型:`0`梯形积分，`1`对格林函数进行线性逼近）

- `func`: `0` for Bessel function; `1` for Hankel function (使用的积分基底:`0`贝塞尔函数，`1`汉克尔函数)

- `out`: the output dispersion spectrum （输出的频散谱）

- `num`: the number of threads for cpu version (defult 20) and the device number of gpu version (default 0) （若是cpu版本则对应并行cpu核数，若是gpu版本则对应了gpu id）


#### For earthquakes
#### 对地震数据（或主动源数据）

```
out = ccfj.fj_earthquake(u,r,c,f,fstride=1,itype=1,func=0,num =20)
```


- `u`: the records in time domain （时间域地震记录）

- `r`: the list of station distances (unit: m) （震中距）

- `c`: the list of phase velocities you want to calculate （相速度）

- `f`: the list of frequencies. The number of points of f should be consistent with the columns of `uf` （频率，这里频率需要对应输入u中时间域点数对应的频率）

- `fstride`: stride of frequency for output (默认1)

- `itype`: `0` for trapezoidal integral; `1` for linear approximate （积分类型:`0`梯形积分，`1`对格林函数进行线性逼近）

- `func`: `0` for Bessel function; `1` for Hankel function (使用的积分基底:`0`贝塞尔函数，`1`汉克尔函数)

- `out`: the output dispersion spectrum （输出的频散谱）

- `num`: the number of threads for cpu version (defult 20) and the device number of gpu version (default 0) （若是cpu版本则对应并行cpu核数，若是gpu版本则对应了gpu id）

### Mutli-windows F-J method (MWFJ)
### 多窗频率贝塞尔变换法 （MWFJ）

This is mainly for earthquake

该函数主要针对地震或者主动源数据

```
out = ccfj.MWFJ(u,r,c,f,Fs,nwin,winl,winr,taper=0.9,fstride=1,itype=0,func=0, num=20)
```

- `u`: the records in time domain

- `r`: the list of station distances (unit: m)

- `c`: the list of phase velocities you want to calculate

- `f`: the list of frequencies. The number of points of f should be consistent with the columns of `uf`

- `Fs`: The sample frequency

- `nwin`: number of time windows （时间窗个数）

- `winl`: list of left side of time windows （每个记录时间窗的左端）

- `winr`: list of right side of time windows （每个记录时间窗的右端）

- `fstride`: stride of frequency for output

- `itype`: `0` for trapezoidal integral; `1` for linear approximate

- `func`: `0` for Bessel function; `1` for Hankel function

- `out`: the output dispersion spectrum

- `num`: the number of threads for cpu version (defult 20) and the device number of gpu version (default 0)

## Uninstall
## 卸载

```
make clean
make uninstall
```


## References
## 参考文献
Wang, J., Wu, G., & Chen, X. (2019). Frequency‐Bessel Transform Method for Effective Imaging of Higher‐Mode Rayleigh Dispersion Curves From Ambient Seismic Noise Data. Journal of Geophysical Research: Solid Earth, 124(4), 3708-3723. doi:10.1029/2018jb016595

Wu, G.-x., Pan, L., Wang, J.-n., & Chen, X. (2020). Shear Velocity Inversion Using Multimodal Dispersion Curves From Ambient Seismic Noise Data of USArray Transportable Array. Journal of Geophysical Research: Solid Earth, 125(1), e2019JB018213. doi:10.1029/2019jb018213

Li, Z., & Chen, X., (2020). An Effective Method to Extract Overtones of Surface Wave from Array Seismic Records of Earthquake Events. Journal of Geophysical Research: Solid Earth, 125(3), e2019JB18511. doi:10.1029/2019jb018511

Li, Z., Zhou, J., Wu, G., Wang, J., Zhang, G., Dong, S., Pan, L., Yang, Z., Gao, L., Ma, Q., Ren, H., & Chen, X. (2021). CC-FJpy: A Python Package for seismic ambient noise cross-correlation and the frequency-Bessel transform method. Seismological Research Letters. doi:10.1785/0220210042

Xi, C., Xia, J., Mi, B., Dai, T., Liu, Y., & Ning, L. (2021). Modified frequency–Bessel transform method for dispersion imaging of Rayleigh waves from ambient seismic noise. Geophysical Journal International, 225(2), 1271-1280. doi:10.1093/gji/ggab008.
