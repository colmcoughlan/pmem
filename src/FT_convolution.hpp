/*
    pmem_heads.hpp is a header file used by pmem and some programs it links with
    Copyright (C) 2012  Colm Coughlan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <string>
#include <fftw3.h>
#include <omp.h>
#include <cmath>

void arrange_ft(double* arr, int imsize);
int ft(fftw_complex* in, fftw_complex* out, int imsize, int direction);
int convolve(double* data, fftw_complex* response, int imsize, int pad_factor, double* output , fftw_plan& forward_transform, fftw_plan& backward_transform, double* double_buff, fftw_complex* complex_buff);
int gen_gauss(double* matrix, int imsize, double cellsize, double bmaj, double bmin, double bpa, int* peak);
int ft_beam(double* beam, fftw_complex* ft_beam, int imsize, int pad_factor, fftw_plan& plan, double* double_buff, fftw_complex* complex_buff);



