# Set compiler to g++.
CC=g++-5
LINK=g++-5
# Set options for the compiler
CCFLAGS= -c -O3 -fopenmp -I/usr/local/Cellar/cfitsio/3.370/include/ -I/Users/admin/git/quickfits
LINKOPTS= -L/usr/local/Cellar/cfitsio/3.370/lib -lcfitsio -L/Users/admin/git/quickfits -lm -lfftw3 -fopenmp -llapack -lblas -lfftw3_threads -lquickfits -O3

all: read_driver FT_convolution pmem update_Lagrangian CE_functions map_stats link

read_driver: src/read_driver.cpp
	$(CC) $(CCFLAGS) src/read_driver.cpp


FT_convolution: src/FT_convolution.cpp
	$(CC) $(CCFLAGS) src/FT_convolution.cpp
	
update_Lagrangian: src/update_Lagrangian.cpp
	$(CC) $(CCFLAGS) src/update_Lagrangian.cpp
	
CE_functions: src/CE_functions.cpp
	$(CC) $(CCFLAGS) src/CE_functions.cpp
	
map_stats: src/map_stats.cpp
	$(CC) $(CCFLAGS) src/map_stats.cpp

pmem: src/pmem.cpp
	$(CC) $(CCFLAGS) src/pmem.cpp

link:
	$(LINK) -o pmem pmem.o read_driver.o FT_convolution.o update_Lagrangian.o CE_functions.o map_stats.o $(LINKOPTS)
	rm -rf *.o

clean:
	rm -rf *.o pmem
