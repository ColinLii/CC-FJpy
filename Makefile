FFTW_VERSION=3.3.9
FFTW_FOLDER=fftw-$(FFTW_VERSION)
PREFIX=$(CONDA_PREFIX)

all:
	python setup.py build_ext --inplace

cpu:
	python setupCPU.py build_ext --inplace
fftw:
	wget http://fftw.org/fftw-$(FFTW_VERSION).tar.gz
	tar -xvf fftw-$(FFTW_VERSION).tar.gz
	rm fftw-$(FFTW_VERSION).tar.gz
	cd $(FFTW_FOLDER) \
	&& ./configure --prefix=$(PWD) \
		--enable-single --enable-shared --with-pic \
		--enable-avx2 --enable-avx --enable-sse --enable-sse2 \
	&& make -j4 install
install:
	install -m 645 ccfj**.so $(PREFIX)/lib/python3*
	install -m 645 lib/libfftw3f* $(PREFIX)/lib

uninstall:
	rm $(PREFIX)/lib/python3*/ccfj**.so
	rm $(PREFIX)/lib/libfftw3f*

clean:
	rm -rf src/CrossCorrelation.c
	rm -rf src/ccfj-gpu.cpp
	rm -rf src/ccfj.cpp
	rm -rf include/fftw3*
	rm -rf build
	rm -rf lib
	rm -rf bin
	rm -rf share
	rm -rf *so
	rm -rf fftw-$(FFTW_VERSION).tar.gz
	rm -rf fftw-$(FFTW_VERSION)