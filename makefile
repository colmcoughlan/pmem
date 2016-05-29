# Set compiler to g++.
CC=g++-5
LINK=g++-5
# Set options for the compiler
CCFLAGS= -c -O3 -fopenmp -I/usr/local/Cellar/cfitsio/3.370/include/ -I/Users/admin/git/quickfits
LINKOPTS= -L/usr/local/Cellar/cfitsio/3.370/lib -lcfitsio -L/Users/admin/git/quickfits -lm -lfftw3 -fopenmp -llapack -lblas -lfftw3_threads -lquickfits -O3

all: read_driver FT_convolution update_Lagrangian CE_functions map_stats variable_metric_methods gradient_methods pmem link

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
	
variable_metric_methods: src/variable_metric_methods.cpp
	$(CC) $(CCFLAGS) src/variable_metric_methods.cpp
	
gradient_methods: src/gradient_methods.cpp
	$(CC) $(CCFLAGS) src/gradient_methods.cpp

pmem: src/pmem.cpp
	$(CC) $(CCFLAGS) src/pmem.cpp

link:
	$(LINK) -o pmem pmem.o read_driver.o FT_convolution.o update_Lagrangian.o CE_functions.o map_stats.o gradient_methods.o variable_metric_methods.o $(LINKOPTS)
	rm -rf *.o

clean:
	rm -rf *.o pmem
