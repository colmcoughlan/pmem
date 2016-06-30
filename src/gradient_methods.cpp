#include "gradient_methods.hpp"
	using namespace std;
	
// return the entropy

double h_func(double** model, double* mask, double* default_map, int npol, int imsize2)
{
	int i,k;
	double m, h;
	
	h = 0.0;
	for(i=0;i< imsize2;i++)
	{
		if(mask[i] > 0)
		{
			m = 0.0;
			for(k=1;k<npol;k++)
			{
				m +=  pow( model[k][i], 2.0);
			}
			m = sqrt(m) / model[0][i];

			h += model[0][i] * ( log(2.0*model[0][i]/(default_map[i]*exp(1.0))) + 0.5*(1.0+m)*log(0.5*(1.0+m)) + 0.5*(1.0-m)*log(0.5*(1.0-m)) );
		}
	}
	
	return(-h);
}

// calculate the step without using the Hessian

double cal_step_sd(double** model , double** residual, double* mask, double* default_map2 , double alpha , double beta , double gamma , int imsize, int ignore_pixels , int npol, double q , double** step_map)
{
	double p , plog , m , dh;
	int i, j, k, ctr;
	int imsize2 = imsize * imsize;
	int right_pixel_limit = imsize - ignore_pixels;
	double j0;
	
	j0 = 0.0;
	#pragma omp parallel for  collapse(2) reduction(+:j0) private(ctr, k,p,m,plog,dh)
	for(i=ignore_pixels;i<right_pixel_limit;i++)
	{
		for(j=ignore_pixels;j<right_pixel_limit;j++)
		{
			ctr = i * imsize + j;
			// Find the total polarised flux, P
			if(mask[ctr] > 0)
			{
				p = 0.0;
				for(k=1;k<npol;k++)
				{
					p +=  pow( model[k][ctr], 2.0);
				}
				p = sqrt(p);

			/*
					Need to implement dH and gradJ according to the MEM formulae
					To do this efficiently, make use of plog and m to prevent calculating the same thing more than once
					If the total polarised flux is sufficiently small, use the Taylor expansion of the logs to save computation
			*/

				if(p > 0.01 * model[0][ctr] )
				{
					plog = 0.5 * log( ( model[0][ctr] - p ) / ( model[0][ctr] + p) ) / p;
				}
				else
				{
					m = p / model[0][ctr];
					plog = - (1.0 + (m * m / 3.0) ) / model[0][ctr];
				}

				for(k=0;k<npol;k++)
				{
					if(k==0)
					{
						dh = -0.5 * log( (model[0][ctr] * model[0][ctr] - p * p ) / (default_map2[ctr]));
						step_map[k][ctr] = dh - 2.0 * alpha * residual[0][ctr] - gamma;
					}
					else
					{
						dh = model[k][ctr] * plog;
						step_map[k][ctr] = dh - 2.0 * beta * residual[k][ctr];
					}
					j0 += step_map[k][ctr] * step_map[k][ctr];
				}
			}
		}
	}
	return(sqrt(j0));
}


// return the total flux

double get_total_flux(double* model, double* mask, int imsize, int ignore_pixels)
{
	int i,j,ctr;
	int right_pixel_limit = imsize - ignore_pixels;
	double flux;
	
	flux = 0.0;
	#pragma omp parallel for  collapse(2) reduction(+:flux) private(ctr)
	for(i=ignore_pixels;i<right_pixel_limit;i++)
	{
		for(j=ignore_pixels;j<right_pixel_limit;j++)
		{
			ctr = i * imsize + j;
			flux += model[ctr];
		}
	}
	
	return(flux);
}

// Find the conjugate gradient (see numerical recipes)

int find_conj_grad(double** hnew, double** hold, double** gnew, double** gold, double step_length, int imsize2, int npol)
{
	double gamma_nom;
	double gamma_denom;
	double gamma_old[4];
	int i, k;
	
	int imsize = sqrt(imsize2);

	for(k=0;k<npol;k++)
	{
		gamma_nom = 0.0;
		gamma_denom = 0.0;
		#pragma omp parallel for reduction(+:gamma_nom, gamma_denom)
		for(i=0;i<imsize2;i++)
		{
			gamma_nom += (gnew[k][i]-gold[k][i])*gnew[k][i];//gnew[k][i] * gnew[k][i];
			gamma_denom += gold[k][i] * gold[k][i];
		}
		gamma_old[k] = gamma_nom / gamma_denom;
	}

	#pragma omp parallel for collapse(2)
	for(k=0;k<npol;k++)
	{
		for(i=0;i<imsize2;i++)
		{
			hnew[k][i] = gnew[k][i] + gamma_old[k] * hold[k][i];
		}
	}
	
	return(0);
}

// steepest descent method (no Hessian, very slow, but stable)

int steepest_descent(double** current_model, double** new_model, double** current_residuals, double** new_residuals, double* default_map2, double* mask, double** convolved_model, fftw_complex* dirty_beam_ft, fftw_complex* complex_buff, double* double_buff, int pad_factor, fftw_plan& forward_transform, fftw_plan& backward_transform, double** dirty_map, double& alpha, double& beta, double& gamma, double& alpha_old, double& beta_old, double& gamma_old, double zsf, int imsize, int ignore_edge_pixels, double& imin, double& imax, double min_flux, int npol, double& q, double pol_upweight_factor, double* chi2_rms, double* rms_theoretical, double* current_rms, int* noise_box, double& old_step_length1, double& old_rescale_step, bool conserve_flux, bool debug, int ctr)
{
	double delta_E, delta_F, delta_G;
	double J0, J1;
	double step_length1, step_limit, rescale_step;
	double total_flux;
	int i, k, err;
	int imsize2 = imsize * imsize;
	
	double eps = 0.01;
	step_limit = 1.0;

	cout<<"SD: Iteration "<<ctr<<endl;
	J0 = cal_step_sd( current_model , current_residuals, mask , default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol , q , new_model);
	cout<<"\tJ0 = "<<J0<<endl;
	step_length1 = eps;
	J0 = J0 * step_length1;
	cout<<"\tCurrent centre = "<<current_model[0][(imsize2+imsize)/2]<<endl;
	cout<<"\tRaw suggested change = "<<new_model[0][(imsize2+imsize)/2]<<endl;
	err = take_step( current_model , new_model , step_length1 , 0.0 , imsize, ignore_edge_pixels , npol, total_flux, 0.0, zsf, min_flux);	// take the step calculated, but scaled by step_length1
	for(k=0;k<npol;k++)
	{
		convolve( new_model[k] , dirty_beam_ft , imsize , pad_factor  , convolved_model[k] , forward_transform , backward_transform , double_buff , complex_buff);
		chi2_rms[k] = get_residual_map(dirty_map[k] , convolved_model[k] , new_residuals[k] , imsize, ignore_edge_pixels);
	}
	cout<<"\tStokes I partial chi2 = "<<chi2_rms[0]<<endl;
	J1 = check_step( current_model , new_model , new_residuals, mask , default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol , q );
	cout<<"\tJ1 = "<<J1<<endl;


	cout<<"\tAlpha = "<<alpha<<endl;
	cout<<"\tBeta = "<<beta<<endl;
	cout<<"\tGamma = "<<gamma<<endl;
	
	// see what rescaling might help dJ = 0
	
	if( J0 - J1 != 0.0 )//	JX = gradJ at point X. Need gradJ = 0; look at J0 and J1 and try predict a rescaling factor that will result in gradJ=0
	{
		rescale_step = J0 / (J0 - J1);	// Derive this formula by finding the slope of J0 to J1 w.r.t. step, and using it to predict the step change which will set gradJ = JN = 0
		if(rescale_step < 0.0)
		{
			rescale_step = 1.0;
		}
		else
		{
			rescale_step = min(rescale_step, 2.0);
		}
	}
	else
	{
		rescale_step = 1.0;
	}
	rescale_step = 0.5 * (rescale_step + old_rescale_step);
	rescale_step = min( rescale_step , step_limit / step_length1 );
	old_rescale_step = rescale_step;
	cout<<"\tRescale step = "<<rescale_step<<endl;
	
	// rescale if worth it
	if( fabs( rescale_step - 1.0) > 0.05 )	// if step 1 was okay, just use it, but if step 2 offers a decent advantage take the average of step1 and step2
	{
		err = interpolate_models( current_model , new_model , rescale_step , imsize, ignore_edge_pixels , npol, min_flux);	// interpolate between old and new models, scaling with step2
		if(err!=0)
		{
			cout<<endl<<"Error detected in interpolation of new and old models, err = "<<err<<endl<<endl;
			return(1);
		}

		err = interpolate_residuals( current_residuals , new_residuals, chi2_rms , rescale_step , imsize2 , npol);	// interpolate residuals too
		if(err!=0)
		{
			cout<<endl<<"Error detected in interpolation of new and old residuals, err = "<<err<<endl<<endl;
			return(1);
		}
	}
	else
	{
		err = copy_model( current_model , new_model );	// just copy new model into old model
		if(err!=0)
		{
			cout<<endl<<"Error detected replacing new and old models, err = "<<err<<endl<<endl;
			return(1);
		}

		err = copy_model( current_residuals , new_residuals );	// just copy new residuals into old residuals
		if(err!=0)
		{
			cout<<endl<<"Error detected replacing new and old residuals, err = "<<err<<endl<<endl;
			return(1);
		}
	}
	// update alpha, beta, gamma. Make sure delta_E, delta_F, delta_G are in the same units
	
	total_flux =  get_total_flux(current_model[0], mask, imsize, ignore_edge_pixels);
	for(k=0;k<npol;k++)
	{
		current_rms[k] = rms_region( current_residuals[k] , noise_box[0] , noise_box[1] , noise_box[2] , noise_box[3] , imsize);	// this might actually go up at the start, as the residual map is not a noise map initially
	}
	cout<<"\tTotal flux = "<<total_flux<<endl;

	delta_E = ( current_rms[0] - rms_theoretical[0]);	// Alpha should "encourage" chi2_I to have rms_theorertical (dalpha = 0 at Jmax)
	delta_F = 0.0;
	for(k=1;k<npol;k++)
	{
		delta_F += pow(current_rms[k] - rms_theoretical[k], 2.0);
	}
	delta_F = sqrt( delta_F);
	delta_F *= pol_upweight_factor;	// upweight polarisation error
	
	delta_G = (total_flux - zsf);	// delta G is just the flux offset. gamma slows down changes when the flux is off
	
//		hval = h_func(current_model, mask, default_map, npol, imsize2);
//		jval = hval - (alpha * chi2_rms[0]/q) - (beta * pol_upweight_factor * temp / q) - gamma * delta_G;

	alpha = max(alpha + delta_E,0.0);
	beta = max(beta + delta_F, 0.0);
	gamma = gamma + delta_G;	// gamma should be allowed go negative (gamma = 0 at ideal solution)

	return(0);
}

// conjugate gradient method (no Hessian, very slow, but stable)
// should be faster than steepest descent, but in this case only by a tiny bit
// this is likely due to all the clipping that needs to be applied to steps (keeping Stokes I >0, 0<m<1 etc.)

int conj_gradient(double** current_model, double** new_model, double** current_residuals, double** new_residuals, 
double* default_map2, double* mask, double** convolved_model, fftw_complex* dirty_beam_ft, fftw_complex* complex_buff, 
double* double_buff, int pad_factor, fftw_plan& forward_transform, fftw_plan& backward_transform, double** dirty_map, 
double& alpha, double& beta, double& gamma, double& alpha_old, double& beta_old, double& gamma_old, double zsf, int imsize, 
int ignore_edge_pixels, double& imin, double& imax, double min_flux, int npol, double& q, double pol_upweight_factor, 
double* chi2_rms, double* rms_theoretical, double* current_rms, int* noise_box, double& old_step_length1, 
double& old_rescale_step, bool conserve_flux, bool debug, int ctr, double** hold, double** hnew, double** gold, double** gnew)
{
	double delta_E, delta_F, delta_G;
	double J0, J1;
	double step_length1, step_limit, rescale_step;
	double total_flux;
	int i, k, err;
	int imsize2 = imsize * imsize;
	
	double eps = 0.01;
	step_limit = 1.0;

	cout<<"CG: Iteration "<<ctr<<endl;
	J0 = cal_step_sd( current_model , current_residuals, mask , default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol , q , gnew);
	
	cout<<"\tJ0 = "<<J0<<endl;
	step_length1 = eps;
	J0 = J0 * step_length1;
	if(ctr!=1)
	{
		find_conj_grad(hnew, hold, gnew, gold, step_length1, imsize2, npol);	// get the conjugate gradient (use this instead of the usual one!)
	}
	else
	{
		for(k=0;k<npol;k++)	// this needs to be a hard copy to keep gnew initialised
		{
			for(i=0;i<imsize2;i++)
			{
				hnew[k][i] = gnew[k][i];
			}
		}
	}
	
	
	cout<<"\tCurrent centre = "<<current_model[0][(imsize2+imsize)/2]<<endl;
	cout<<"\tRaw suggested change = "<<hnew[0][(imsize2+imsize)/2]<<endl;
	err = take_step_nd( current_model, hnew , new_model , step_length1 , 0.0 , imsize, ignore_edge_pixels , npol, total_flux, 0.0, zsf, min_flux);	// take the step calculated, but scaled by step_length1
	for(k=0;k<npol;k++)
	{
		convolve( new_model[k] , dirty_beam_ft , imsize , pad_factor  , convolved_model[k] , forward_transform , backward_transform , double_buff , complex_buff);
		chi2_rms[k] = get_residual_map(dirty_map[k] , convolved_model[k] , new_residuals[k] , imsize, ignore_edge_pixels);
	}
	cout<<"\tStokes I unweighted chi2 = "<<chi2_rms[0]<<endl;
	J1 = check_step( current_model , new_model , new_residuals, mask , default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol , q );
	cout<<"\tJ1 = "<<J1<<endl;


	cout<<"\tAlpha = "<<alpha<<endl;
	cout<<"\tBeta = "<<beta<<endl;
	cout<<"\tGamma = "<<gamma<<endl;
	
	// see what rescaling might help dJ = 0
	
	if( J0 - J1 != 0.0 )//	JX = gradJ at point X. Need gradJ = 0; look at J0 and J1 and try predict a rescaling factor that will result in gradJ=0
	{
		rescale_step = J0 / (J0 - J1);	// Derive this formula by finding the slope of J0 to J1 w.r.t. step, and using it to predict the step change which will set gradJ = JN = 0
		if(rescale_step < 0.0)
		{
			rescale_step = 1.0;
		}
		else
		{
			rescale_step = min(rescale_step, 2.0);
		}
	}
	else
	{
		rescale_step = 1.0;
	}
	rescale_step = 0.5 * (rescale_step + old_rescale_step);
	rescale_step = min( rescale_step , step_limit / step_length1 );
	old_rescale_step = rescale_step;
	cout<<"\tRescale step = "<<rescale_step<<endl;
	
	// rescale if worth it
	if( fabs( rescale_step - 1.0) > 0.05 )	// if step 1 was okay, just use it, but if step 2 offers a decent advantage take the average of step1 and step2
	{
		err = interpolate_models( current_model , new_model , rescale_step , imsize, ignore_edge_pixels , npol, min_flux);	// interpolate between old and new models, scaling with step2
		if(err!=0)
		{
			cout<<endl<<"Error detected in interpolation of new and old models, err = "<<err<<endl<<endl;
			return(1);
		}

		err = interpolate_residuals( current_residuals , new_residuals, chi2_rms , rescale_step , imsize2 , npol);	// interpolate residuals too
		if(err!=0)
		{
			cout<<endl<<"Error detected in interpolation of new and old residuals, err = "<<err<<endl<<endl;
			return(1);
		}
	}
	else
	{
		err = copy_model( current_model , new_model );	// just copy new model into old model
		if(err!=0)
		{
			cout<<endl<<"Error detected replacing new and old models, err = "<<err<<endl<<endl;
			return(1);
		}

		err = copy_model( current_residuals , new_residuals );	// just copy new residuals into old residuals
		if(err!=0)
		{
			cout<<endl<<"Error detected replacing new and old residuals, err = "<<err<<endl<<endl;
			return(1);
		}
	}
	// update alpha, beta, gamma. Make sure delta_E, delta_F, delta_G are in the same units
	
	total_flux =  get_total_flux(current_model[0], mask, imsize, ignore_edge_pixels);
	for(k=0;k<npol;k++)
	{
		current_rms[k] = rms_region( current_residuals[k] , noise_box[0] , noise_box[1] , noise_box[2] , noise_box[3] , imsize);	// this might actually go up at the start, as the residual map is not a noise map initially
	}
	cout<<"\tTotal flux = "<<total_flux<<endl;

	delta_E = ( current_rms[0] - rms_theoretical[0]);	// Alpha should "encourage" chi2_I to have rms_theorertical (dalpha = 0 at Jmax)
	delta_F = 0.0;
	for(k=1;k<npol;k++)
	{
		delta_F += pow(current_rms[k] - rms_theoretical[k], 2.0);
	}
	delta_F = sqrt( delta_F);
	delta_F *= pol_upweight_factor;	// upweight polarisation error
	
	delta_G = (total_flux - zsf);	// delta G is just the flux offset. gamma slows down changes when the flux is off
	
//		hval = h_func(current_model, mask, default_map, npol, imsize2);
//		jval = hval - (alpha * chi2_rms[0]/q) - (beta * pol_upweight_factor * temp / q) - gamma * delta_G;

	alpha = max(alpha + delta_E,0.0);
	beta = max(beta + delta_F, 0.0);
	gamma = gamma + delta_G;	// gamma should be allowed go negative (gamma = 0 at ideal solution)
	
	// record hold, gold (take step length into account in gamma function)
	
	#pragma omp parallel for collapse(2)
	for(k=0;k<npol;k++)
	{
		for(i=0;i<imsize2;i++)
		{
			hold[k][i] = hnew[k][i];
			gold[k][i] = gnew[k][i];
		}
	}

	return(0);
}