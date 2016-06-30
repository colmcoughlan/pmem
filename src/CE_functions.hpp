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
#ifndef grad_def
	#include "gradient_structure.hpp"
#endif
	using namespace std;

double get_residual_map(double* dirty_map , double* convolved_model , double* residual_map, int imsize, int ignore_pixels);
int get_info(double** model , double** residual, double* mask, double* default_map2 , gradient_structure& grad, 
double& total_flux , double alpha , double beta , double gamma , double& imin, double& imax , int imsize, int ignore_pixels , 
int npol, double q, double** hessian_buffer, bool keep_hessian);
int cal_step(double** model , double** residual, double* mask, double* default_map2 , double alpha , double beta , double gamma , int imsize, int ignore_pixels , int npol , double q , double& J0 , double** step_map);
int take_step(double** model , double** step , double step_length , double step_limit , int imsize, int ignore_pixels , int npol, double total_flux, double max_flux_change, double zsp, double min_flux);
int take_step_nd(double** model , double** step , double** new_model , double step_length , double step_limit , int imsize, int ignore_pixels , int npol, double total_flux, double max_flux_change, double zsp, double min_flux);
double check_step(double** old_model , double** new_model , double** new_residual, double* mask, double* default_map2 , double alpha , double beta , double gamma , int imsize, int ignore_pixels , int npol , double q);
int interpolate_models(double** current_model , double** new_model , double frac_new , int imsize, int ignore_pixels , int npol, double min_flux);
int interpolate_models_residuals(double** current_model , double** new_model, double** current_residuals, double**new_residuals , double frac_new , int imsize, int ignore_pixels , int npol, double min_flux, double* chi2_rms);
int interpolate_residuals(double** current_residuals , double** new_residuals, double* chi2_rms , double frac_new , int imsize2 , int npol);
int copy_model(double**& model1, double**& model2);
double find_j0(double** model , double** residual, double* mask, double* default_map2 , double alpha , double beta , double gamma , int imsize, int ignore_pixels , int npol, double q , double** new_mod, double** old_mod);
int find_gradj(double** model , double** residual, double* mask, double* default_map2 , double alpha , double beta , double gamma , int imsize, int ignore_pixels , int npol, double** gradj);