# Set compiler to g++.
CC=g++-mp-4.9
LINK=g++-mp-4.9
# Set options for the compiler
CCFLAGS= -c -O3 -fopenmp -I/Users/admin/git/quickfits
LINKOPTS= -L/Users/admin/fftw/lib -L/opt/local/lib/ -L/Users/admin/git/quickfits  -lm -lfftw3 -fopenmp -llapack -lblas -lfftw3_threads -lcfitsio -lquickfits -O3

all: read_driver FT_convolution pmem update_Lagrangian CE_functions map_stats link

read_driver: pmem_src/read_driver.cpp
	$(CC) $(CCFLAGS) pmem_src/read_driver.cpp


FT_convolution: pmem_src/FT_convolution.cpp
	$(CC) $(CCFLAGS) pmem_src/FT_convolution.cpp
	
update_Lagrangian: pmem_src/update_Lagrangian.cpp
	$(CC) $(CCFLAGS) pmem_src/update_Lagrangian.cpp
	
CE_functions: pmem_src/CE_functions.cpp
	$(CC) $(CCFLAGS) pmem_src/CE_functions.cpp
	
map_stats: pmem_src/map_stats.cpp
	$(CC) $(CCFLAGS) pmem_src/map_stats.cpp

pmem: pmem_src/pmem.cpp
	$(CC) $(CCFLAGS) pmem_src/pmem.cpp

link:
	$(LINK) -o pmem pmem.o read_driver.o FT_convolution.o update_Lagrangian.o CE_functions.o map_stats.o $(LINKOPTS)

clean:
	rm -rf *.o nmem
