# -*- coding: utf-8 -*-

"""
setup.py: compile the CCFJ-Py in GPU or CPU version 
CCFJ-Py: A Python Package for seismic ambient noise cross-correlation and the frequency-Bessel transform method.

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
from os.path import join as pjoin

#os.environ["CC"] = 'gcc-9'
#os.environ["CXX"] = 'g++-9'

sources = ['src/FJgpu.cu','src/ccfj-gpu.pyx','src/CrossCorr.cpp']


def find_in_path(name, path):
    """Find a file in a search path"""

    # Adapted fom http://code.activestate.com/recipes/52224
    for dir in path.split(os.pathsep):
        binpath = pjoin(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None


def locate_cuda():
    """Locate the CUDA environment on the system

    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.

    Starts by looking for the CUDAHOME env variable. If not found,
    everything is based on finding 'nvcc' in the PATH.
    """

    # First check if the CUDAHOME env variable is in use
    if 'CUDAHOME' in os.environ:
        home = os.environ['CUDAHOME']
        nvcc = pjoin(home, 'bin', 'nvcc')
    else:
        # Otherwise, search the PATH for NVCC
        nvcc = find_in_path('nvcc', os.environ['PATH'])
        if nvcc is None:
            print('The nvcc binary could not be '
                'located in your $PATH. Either add it to your path, '
                'or set $CUDAHOME')
            return None
        home = os.path.dirname(os.path.dirname(nvcc))

    cudaconfig = {'home': home, 'nvcc': nvcc,
                  'include': pjoin(home, 'include'),
                  'lib64': pjoin(home, 'lib64')}
    for k, v in iter(cudaconfig.items()):
        if not os.path.exists(v):
            raise EnvironmentError('The CUDA %s path could not be '
                                    'located in %s' % (k, v))

    return cudaconfig


def customize_compiler_for_nvcc(self):
    """Inject deep into distutils to customize how the dispatch
    to gcc/nvcc works.

    If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going
    the OO route, I have this. Note, it's kindof like a wierd functional
    subclassing going on.
    """

    # Tell the compiler it can processes .cu
    self.src_extensions.append('.cu')

    # Save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile

    # Now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cu':
            # use the cuda for .cu files
            self.set_executable('compiler_so', CUDA['nvcc'])
            # use only a subset of the extra_postargs, which are 1-1
            # translated from the extra_compile_args in the Extension class
            postargs = extra_postargs['nvcc']
        else:	
            self.set_executable('compiler_so', CUDA['nvcc'])
            postargs = extra_postargs['nvcc']

        super(obj, src, ext, cc_args, postargs, pp_opts)
        # Reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # Inject our redefined _compile method into the class
    self._compile = _compile



# Run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)



CUDA = locate_cuda()

if CUDA == None:
    print('-------------------------------------')
    print('nvcc from CUDA not found.............')
    print('Compiling the CPU Version of ccfj....')
    print('-------------------------------------')
    libs = ['m', 'fftw3f']
    args = ['-std=c99', '-O3','-DuseOMP']
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
    
else:
    print('-------------------------------------')
    print('Compiling the GPU Version of ccfj....')
    print('-------------------------------------')
    # Obtain the numpy include directory. This logic works across numpy versions.
    try:
        numpy_include = numpy.get_include()
    except AttributeError:
        numpy_include = numpy.get_numpy_include()


    ext = Extension('ccfj',
        sources = sources,
        library_dirs = [CUDA['lib64'],'lib'],
        libraries = ['cudart','fftw3f'],
        language = 'c++',
        runtime_library_dirs = [CUDA['lib64']],
        extra_link_args=['-Wl,-rpath,lib'],
        # This syntax is specific to this build system
        # we're only going to use certain compiler args with nvcc
        # and not with gcc the implementation of this trick is in
        # customize_compiler()
        extra_compile_args= {
          'gcc': ['-std=c99', '-O3','-DuseOMP','-fopenmp'],
          'nvcc': [
            '--ptxas-options=-v', '-c',
            '--compiler-options', "'-fPIC'"
            ]
          },
          include_dirs = [numpy_include, CUDA['include'], 'src','include']
        )



    setup(name='ccfj',
      version='0.1',
      author='Zhengbo Li et al.',
      packages=['ccfj'],

      ext_modules = [ext],
      # Inject our custom trigger
      cmdclass = {'build_ext': custom_build_ext},

      # Since the package has c code, the egg cannot be zipped
      zip_safe = False
      )

