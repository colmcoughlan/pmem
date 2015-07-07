# Set compiler to g++.
CC=g++-mp-4.9
LINK=g++-mp-4.9
# Set options for the compiler
CCFLAGS= -c -O3 -fopenmp -I/Users/admin/git/quickfits
LINKOPTS= -L/Users/admin/fftw/lib -L/opt/local/lib/ -L/Users/admin/git/quickfits  -lm -lfftw3 -fopenmp -llapack -lblas -lfftw3_threads -lcfitsio -lquickfits -O3

all: pmem link

pmem: pmem.cpp
	$(CC) $(CCFLAGS) pmem.cpp

link:
	$(LINK) -o pmem pmem.o $(LINKOPTS)

clean:
	rm -rf *.o pmem
