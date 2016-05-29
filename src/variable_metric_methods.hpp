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


#include <cmath>
#include <iostream>
#include "CE_functions.hpp"
#include "FT_convolution.hpp"
#include "map_stats.hpp"
#include <fftw3.h>

// define BLAS dot product

extern "C" {
	double ddot_( int * , double * , int * , double * , int * );
};

// warning : make sure to pass pointers to current and new models and residuals by reference, since they swap sometimes!

int newton_raphson(double**& current_model, double**& new_model, double**& current_residuals, double**& new_residuals, 
double* default_map2, double* mask, double** convolved_model, fftw_complex* dirty_beam_ft, fftw_complex* complex_buff, 
double* double_buff, int pad_factor, fftw_plan& forward_transform, fftw_plan& backward_transform, double** dirty_map, 
gradient_structure& grad, double& alpha, double& beta, double& gamma, double& alpha_old, double& beta_old, double& gamma_old, 
bool force_chi2_method, double zsf, int imsize, int ignore_edge_pixels, double& imin, double& imax, double min_flux, int npol, 
double& q, double pol_upweight_factor, double* chi2_rms, double* rms_theoretical, double* current_rms, int* noise_box, 
double acceleration_factor, double& old_step_length1, double& old_rescale_step, bool conserve_flux, bool debug, int ctr);

int quasi_newton(double**& current_model, double**& new_model, double**& current_residuals, double**& new_residuals, 
double* default_map2, double* mask, double** convolved_model, fftw_complex* dirty_beam_ft, fftw_complex* complex_buff, 
double* double_buff, int pad_factor, fftw_plan& forward_transform, fftw_plan& backward_transform, double** dirty_map, 
gradient_structure& grad, double& alpha, double& beta, double& gamma, double& alpha_old, double& beta_old, double& gamma_old, 
bool force_chi2_method, double zsf, int imsize, int ignore_edge_pixels, double& imin, double& imax, double min_flux, int npol, 
double& q, double pol_upweight_factor, double* chi2_rms, double* rms_theoretical, double* current_rms, int* noise_box, 
double acceleration_factor, double& old_step_length1, double& old_rescale_step, bool conserve_flux, bool debug, int ctr, 
double** gradold, double** gradnew, double** hessian, bool do_bfgs);

double dot_product(double** vec1, double** vec2, int imsize2, int npol);
double metric_product(double** vec1, double** vec2, double** matrix, int imsize2, int npol);
int sub_vec(double** vec1, double** vec2, int imsize2, int npol);
int update_hessian(double** model , double** residual, double* mask, double* default_map2 , gradient_structure& grad, 
double& total_flux , double alpha , double beta , double gamma , double& imin, double& imax, double** model_step, 
double** grad_step, double** hessian, int imsize, int ignore_edge_pixels , int npol, bool do_bfgs);
