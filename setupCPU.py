# -*- coding: utf-8 -*-

"""
setupCPU.py: compile the CC-FJpy in CPU version.
CC-FJpy: A Python Package for seismic ambient noise cross-correlation and the frequency-Bessel transform method.

:copyright:
 Xiaofei Chen Research Group, Department of Earth and Space Sciences, SUSTech, China.
:license:
 GNU Lesser General Public License, Version 3
 (https://www.gnu.org/copyleft/lesser.html)
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
import os

# If you have the Intel MKL installed,
# you can use the Intel Link Line advisor
# and replace the fftw3 library with MKL
# once you have compiled the FFTW3 interface
# called fftw3xc

libs = ['m', 'fftw3f']
args = ['-std=c99', '-O3','-DuseOMP', '-fopenmp']
sources = ['src/ccfj.pyx', 'src/CrossCorr.cpp','src/FJcpu.cpp']
include = ['include',numpy.get_include()]
linkerargs = ['-Wl,-rpath,lib','-fopenmp']
libdirs = ['lib']


extensions = [
    Extension("ccfj",
              sources=sources,
              include_dirs=include,
              libraries=libs,
              library_dirs=libdirs,
              extra_compile_args=args,
              language='c++',
              
              extra_link_args=linkerargs)
]

setup(name='ccfj',
      version='0.1',
      author='Zhengbo Li et al.',
      packages=['ccfj'],
      ext_modules=cythonize(extensions),
      )
