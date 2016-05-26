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
#include <fstream>
#include <string>
#include <string.h>
#include <cmath>
#include <stdlib.h>
#include <sstream>
#include <fftw3.h>
#include <fitsio.h>
#include <stdio.h>
#include <omp.h>

#ifndef grad_def
	#include "gradient_structure.hpp"
#endif
#include "read_driver.hpp"
#include "update_Lagrangian.hpp"
#include "CE_functions.hpp"
#include "FT_convolution.hpp"
#include "map_stats.hpp"
#include "steepest_descent.hpp"

extern "C" {
	#include "quickfits.h"
};

// warning : make sure to pass pointers to current and new models and residuals by reference, since they swap sometimes!

int newton_raphson(double**& current_model, double**& new_model, double**& current_residuals, double**& new_residuals, 
double* default_map2, double* mask, double** convolved_model, fftw_complex* dirty_beam_ft, fftw_complex* complex_buff, 
double* double_buff, int pad_factor, fftw_plan& forward_transform, fftw_plan& backward_transform, double** dirty_map, 
gradient_structure& grad, double& alpha, double& beta, double& gamma, double& alpha_old, double& beta_old, double& gamma_old, 
bool force_chi2_method, double zsf, int imsize, int ignore_edge_pixels, double& imin, double& imax, double min_flux, int npol, 
double& q, double pol_upweight_factor, double* chi2_rms, double* rms_theoretical, double* current_rms, int* noise_box, 
double acceleration_factor, double& old_step_length1, double& old_rescale_step, bool conserve_flux, bool debug, int ctr);
