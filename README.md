# Usage

To use just call the executable: ./pmem, or install through CASA (https://github.com/colmcoughlan/pmem/tree/master/CASA_interface)

This software deconvolves FITS files using the Maximum Entropy Method with support for polarisation.

# Installation instructions

Dependencies: a C compiler, fftw3-dev, blas, cfitsio and quickfits (https://github.com/colmcoughlan/quickfits).

On Ubuntu cfitsio can be found in the cfitsio-dev package. On a mac the best way to install it might be through homebrew (brew install cfitsio).
More mac instructions can be found at https://heasarc.gsfc.nasa.gov/fitsio/fitsio_macosx.html. You can read about homebrew at http://brew.sh/.

If you're using linux you probably already have a c compiler (gcc). On a mac you can install this through homebrew with "brew install gcc5", to install gcc version 5.

Once you have the dependencies installed, make sure they are on your library path. On linux this involves adding the directories to $LD_LIBRARY_PATH.

To build, edit the makefiles such that:

  1. link=\<your C compiler>
  2. CCFLAGS= -c -O3 -fopenmp -I/\<path to .h directories> -I/Users/admin/git/quickfits
  3. LINKOPTS=-L/\<path to library directories> -lcfitsio -L/Users/admin/git/quickfits -lm -lfftw3 -fopenmp -llapack -lblas -lfftw3_threads -lquickfits -O3

Note this assumes that cfitsio etc. have been installed onto your path already. If this is not the case, you should specify additional -I and 
-L paths to their location.

Finally, run make -f makefile to generate the executable.

PMEM v1.04
Deconvolves dirty maps using the Maximum Entropy Method. Supports polarisation deconvolution.
Requires quickfits v1.1

1.04

Added in two quasi-newton solvers : DFP and BFGS. Also added a conjugate gradient solver. The newton solver (option 0) remains the fastest, though the quasi newton solvers are also quite good. Reordered the header files and some of the code. Made sure multithreading was enabled everywhere. Added comments.

1.03
Added in option to use a steepest descent-based solver. This is usually slower than the Newton-Raphson solver, but can be faster for the first few iterations. Also includes a bugfix for when the user specified a non-standard default map.

1.02
Fixed bugs related to quickfits v1.1
Turned off residual scaling by beam in code, may expose option to user later
Now uses proper masking on output
