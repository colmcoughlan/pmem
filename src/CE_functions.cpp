#include "CE_functions.hpp"
	using namespace std;

/*
	"get_residual_map"

	gets the residual map of the convolved model - the dirty map (units = Jy/beam)

	inputs:

	dirty_map = dirty map, nothing done to it (units = Jy/beam)
	convolved_model = the current model convolved with the dirty beam (units = Jy/beam)
	size = the number of pixels

	outputs:

	residual_map = difference between convolved model and dirty map
	on return, chi2 for the current Stokes parameter
*/

double get_residual_map(double* dirty_map , double* convolved_model , double* residual_map , int imsize, int ignore_pixels)
{
	int i , j , ctr;
	int right_pixel_limit = imsize - ignore_pixels;
	double chi2 = 0.0;

	for(i=ignore_pixels;i<right_pixel_limit;i++)
	{
		for(j=ignore_pixels;j<right_pixel_limit;j++)
		{
			ctr = i * imsize + j;
			residual_map[ctr] = convolved_model[ctr] - dirty_map[ctr];
			chi2 += (residual_map[ctr] * residual_map[ctr]);
		}
	}

	return(chi2);
}

/*
	"get_info"

	Calculates inner products. Finds the current total flux and total rms residuals

	Inputs:

	model = Current model imap
	residual = Current residual maps (I,Q,U,V)
	default_map2 = Default imap squared (normally flat)
	alpha = Lagrange parameter for chi2 for Stokes I
	beta = Lagrange parameter for chi2 for polarisation Stokes parameters
	gamma = Lagrange paramter for flux conservation
	imsize2 = The number of pixels in the images
	npol = The number of polarisations
	q = Factor converting Jy/pix to Jy/beam

	outputs:

	grad = Structure holding all the inner products that need to be calculated
	imin = Minimum Stokes I flux
	imax = Maximum Stokes I flux
	on return = error (0 if none)
*/

int get_info(double** model , double** residual, double* mask, double* default_map2 , gradient_structure& grad, double& total_flux , double alpha , double beta , double gamma , double& imin, double& imax , int imsize, int ignore_pixels , int npol, double q)
{
	double metric,dh,dh2,de,df,dg,temp,p,pexp,plog,m;
	double EE , EF , EG , EH , FF , FG , FH , GG , GH , HH , II;
	int i ,j , k , ctr , err;
	int imsize2 = imsize * imsize;
	int right_pixel_limit = imsize - ignore_pixels;

	EE = 0.0;
	EF = 0.0;
	EG = 0.0;
	EH = 0.0;
	FF = 0.0;
	FG = 0.0;	// note EF and FG will stay as zero
	FH = 0.0;
	GG = 0.0;
	GH = 0.0;
	HH = 0.0;
	II = 0.0;

	temp=0.0;

	err=0;

//	#pragma omp parallel for collapse(2) reduction( +: EE , EG , EH , FF , FH , GG , GH , HH ,temp) private(ctr,k,p,m,plog,pexp,dh,dh2,de,df,dg,metric)
	for(i=ignore_pixels;i<right_pixel_limit;i++)
	{
		for(j=ignore_pixels;j<right_pixel_limit;j++)
		{
			ctr = i * imsize + j;
			if(mask[ctr]>0)
			{
				p = 0.0;
				for(k=1;k<npol;k++)
				{
					p += (model[k][ctr] * model[k][ctr]);	// find polarised intensity for each pixel
				}
				p = sqrt(p);


				if(p > 0.01 * model[0][ctr])
				{
					plog = 0.5 * log( ( model[0][ctr] - p ) / ( model[0][ctr] + p) ) / p;
					pexp = ( ( model[0][ctr] / (p * p - model[0][ctr] * model[0][ctr] ) ) - plog ) / (p * p);	// calculate some expressions needed for evaluating the derivative of the entropy
				}
				else
				{
					m = p / model[0][ctr];
					plog = - (1.0 + (m * m / 3.0) ) / model[0][ctr];
					pexp = - ( (2.0 / 3.0) + 0.8 * m * m) / (model[0][ctr] * model[0][ctr] * model[0][ctr]);	// Taylor series expansion of above for small m (see Mathematica)
				}


				for(k=0;k<npol;k++)
				{
					if(k==0)	// Stokes I only
					{
						dh = -0.5 * log( ( model[0][ctr] * model[0][ctr] - p * p ) / ( default_map2[ctr] ) );	// derivative of the entropy contribution of pixel
						dh2 = - model[0][ctr] / ( model[0][ctr] * model[0][ctr] - p * p);					// second derivative of the entropy contribution of pixel
						de = 2.0 * residual[0][ctr];								// derivative of the Stokes I chi2 contribution of the pixel
						df = 0.0;										// derivative of the Stokes Q,U,V chi2 contribution of the pixel
						dg = 1.0;										// derivative of the flux conservation term contribution of the pixel
						metric = 1.0/( 2.0 * q * alpha - dh2);							// metric of the minimisation at this point
/*
						if(ctr == 524800)
						{
							cout<<"dh , dh2, de, metric = "<<dh<<","<<dh2<<" , "<<de<<" , "<<metric<<endl;
						}
*/
					}
					else	// Other Stokes parameters (same idea)
					{
						dh = model[k][ctr] * plog;
						dh2 = plog + model[k][ctr] * model[k][ctr] * pexp;
						de = 0.0;
						df = 2.0 * residual[k][ctr];
						dg = 0.0;
						metric = 1.0/( 2.0 * q * beta - dh2);
/*
						if(ctr == 524800)
						{
							cout<<"k = "<<k<<": dh , dh2, df, metric = "<<dh<<","<<dh2<<" , "<<de<<" , "<<metric<<endl;
						}
*/
					}

					EE += de*metric*de;
					EG += de*metric*dg;
					EH += de*metric*dh;	// contributions to the inner products
					FF += df*metric*df;
					FH += df*metric*dh;
					GG += dg*metric*dg;
					GH += dg*metric*dh;
					HH += dh*metric*dh;
					temp += metric;	// inner product 1.1
				}
			}
		}
	}



	total_flux=0.0;
	imin=9999.0;
	imax=-9999.0;

	for(i = 0;i<imsize2;i++)	// find max, min and total flux
	{
			total_flux += model[0][i];
			imin=min(imin,model[0][i]);
			imax=max(imax,model[0][i]);
	}


	// save inner products into the structure

	grad.EJ = EH - alpha * EE - beta * EF - gamma * EG;
	grad.FJ = FH - alpha * EF - beta * FF - gamma * FG;
	grad.GJ = GH - alpha * EG - beta * FG - gamma * GG;
	grad.JJ = HH + alpha * alpha * EE + beta * beta * FF + gamma * gamma * GG - 2.0 * alpha * EH - 2.0 * beta * FH - 2.0 * gamma * GH + 2.0 * alpha * beta * EF + 2.0 * alpha * gamma * EG + 2.0 * beta * gamma * FG;
	II = HH + alpha * alpha * EE + beta * beta * FF + gamma * gamma * GG;

	if(II < 0)
	{
		II=temp;
	}

	grad.EE = EE;
	grad.EF = EF;
	grad.EG = EG;
	grad.EH = EH;
	grad.FF = FF;
	grad.FG = FG;
	grad.FH = FH;
	grad.GG = GG;
	grad.GH = GH;
	grad.HH = HH;
	grad.II = II;
	

	return(err);
}

/*
	cal_step

	cal_step finds the next step to take for the image using the Newton-Ralphson method
	It also uses the plog, m and pexp terms to try to be computationally effecient

	inputs:
		model = the model Stokes I, Q and U maps
		residuals = the residual Stokes I, Q and U maps
		default_map2 = the square of the default map
		imsize2 = the number of pixels in the new map
		alpha = the Lagrangian parameter associated with minimising Stokes I residuals
		beta = the Lagrangian parameter associated with minimising Stokes Q and U residuals
		gamma = the Lagrangian parameter associated with conserving the Stokes I flux
		q = the number of pixels in a beam
		npol = the number of polarisations being deconvolved

	outputs:
		J0 =  a measure of the change in J due to the step taken
		step_map = the steps for the Stokes I,Q,U model maps
		on return, 0
*/

int cal_step(double** model , double** residual, double* mask, double* default_map2 , double alpha , double beta , double gamma , int imsize, int ignore_pixels , int npol, double q , double& J0 , double** step_map)
{
	double p , plog , pexp , m , dh , dh2 , gradJ , metric , step;
	int i, j, k, ctr;
	int imsize2 = imsize * imsize;
	int right_pixel_limit = imsize - ignore_pixels;

	double J0_temp = 0.0;
	double step_sum;
	
	for(k=0;k<npol;k++)
	{
		for(i=0;i<imsize2;i++)
		{
			step_map[k][i] = 0.0;
		}
	}

	#pragma omp parallel for  collapse(2) reduction(+:J0_temp) private(ctr, k,p,m,plog,pexp,dh,dh2,gradJ,metric,step)	// eff : rewrite to reduce operations (default map^2, 2* q )
	for(i=ignore_pixels;i<right_pixel_limit;i++)
	{
		for(j=ignore_pixels;j<right_pixel_limit;j++)
		{
			ctr = i * imsize + j;
			if(mask[ctr] > 0)
			{
				p = 0.0;
				if(npol>1)
				{
					for(k=1;k<npol;k++)
					{
						p += (model[k][ctr] * model[k][ctr]);
					}
					p = sqrt(p);

					if(p > 0.01 * model[0][ctr])
					{
						plog = 0.5 * log( ( model[0][ctr] - p ) / ( model[0][ctr] + p) ) / p;
						pexp = ( ( model[0][ctr] / (p * p - model[0][ctr] * model[0][ctr] ) ) - plog ) / (p * p);
					}
					else
					{
						m = p / model[0][ctr];
						plog = - (1.0 + (m * m / 3.0) ) / model[0][ctr];
						pexp = - ( (2.0 / 3.0) + 0.8 * m * m) / (model[0][ctr] * model[0][ctr] * model[0][ctr]);	// Taylor series expansion of above for small m (see Mathematica)
					}
				}

				for(k=0;k<npol;k++)
				{
					if(k==0)
					{
						dh = - 0.5 * log( (model[0][ctr] * model[0][ctr] - p * p) / (default_map2[ctr]));
						dh2 = - model[0][ctr] / (model[0][ctr] * model[0][ctr] - p * p);
						gradJ = dh - 2.0 * alpha * residual[0][ctr] - gamma;
						metric = 1.0 / ( 2.0 * q * alpha - dh2);
					}
					else
					{
						dh = model[k][ctr] * plog;
						dh2 = plog + model[k][ctr] * model[k][ctr] * pexp;
						gradJ = dh - 2.0 * beta * residual[k][ctr];
						metric = 1.0/( 2.0 * q * beta - dh2);
					}

					step = metric * gradJ;
			
					// clip step if necessary
		//			step = min(0.00001 , max(-0.00001,step));
			
					J0_temp += gradJ * step;
					step_map[k][ctr] = step;
				}
			}
		}
	}

	J0 = J0_temp;

	return(0);
}

double find_j0(double** model , double** residual, double* mask, double* default_map2 , double alpha , double beta , double gamma , int imsize, int ignore_pixels , int npol, double q , double** new_mod, double** old_mod)
{
	double p , plog , pexp , m , dh , dh2 , gradJ , metric , step;
	int i, j, k, ctr;
	int imsize2 = imsize * imsize;
	int right_pixel_limit = imsize - ignore_pixels;

	double J0 = 0.0;

	#pragma omp parallel for  collapse(2) reduction(+:J0) private(ctr, k,p,m,plog,pexp,dh,dh2,gradJ,metric,step)	// eff : rewrite to reduce operations (default map^2, 2* q )
	for(i=ignore_pixels;i<right_pixel_limit;i++)
	{
		for(j=ignore_pixels;j<right_pixel_limit;j++)
		{
			ctr = i * imsize + j;
			if(mask[ctr] > 0)
			{
				p = 0.0;
				if(npol>1)
				{
					for(k=1;k<npol;k++)
					{
						p += (model[k][ctr] * model[k][ctr]);
					}
					p = sqrt(p);

					if(p > 0.01 * model[0][ctr])
					{
						plog = 0.5 * log( ( model[0][ctr] - p ) / ( model[0][ctr] + p) ) / p;
						pexp = ( ( model[0][ctr] / (p * p - model[0][ctr] * model[0][ctr] ) ) - plog ) / (p * p);
					}
					else
					{
						m = p / model[0][ctr];
						plog = - (1.0 + (m * m / 3.0) ) / model[0][ctr];
						pexp = - ( (2.0 / 3.0) + 0.8 * m * m) / (model[0][ctr] * model[0][ctr] * model[0][ctr]);	// Taylor series expansion of above for small m (see Mathematica)
					}
				}

				for(k=0;k<npol;k++)
				{
					if(k==0)
					{
						dh = - 0.5 * log( (model[0][ctr] * model[0][ctr] - p * p) / (default_map2[ctr]));
						dh2 = - model[0][ctr] / (model[0][ctr] * model[0][ctr] - p * p);
						gradJ = dh - 2.0 * alpha * residual[0][ctr] - gamma;
						metric = 1.0 / ( 2.0 * q * alpha - dh2);
					}
					else
					{
						dh = model[k][ctr] * plog;
						dh2 = plog + model[k][ctr] * model[k][ctr] * pexp;
						gradJ = dh - 2.0 * beta * residual[k][ctr];
						metric = 1.0/( 2.0 * q * beta - dh2);
					}
			
					J0 += gradJ * (new_mod[k][ctr] - old_mod[k][ctr]);
				}
			}
		}
	}

	return(J0);
}

/*
	take_step

	take_step actually takes the step that was generated by cal_step
	it also enforces a minimum flux and makes sure that the fractional polarisation does not leave the region 0<m<1

	inputs:
		model = the model Stokes I, Q and U maps
		step = the steps for the Stokes I, Q and U model maps
		step_length = the factor by which the steps should be scaled
		step_limit = the maximum allowed change
		imsize2 = the number of pixels in a map
		npol = the number of polarisations being deconvolved

	outputs:
		step = the new model I, Q and U maps (after the step)
		on return, 0
*/

int take_step(double** model , double** step , double step_length , double step_limit , int imsize, int ignore_pixels , int npol, double total_flux, double max_flux_change, double zsp, double min_flux)
{

	// First declare variables

	double pnew , inew , factor;
	int i, j, k , ctr ;
	int imsize2 = imsize * imsize;
	double step_sum, max_pixel_change, min_pixel_change;
	int right_pixel_limit = imsize - ignore_pixels;

	// Start openmp for loop if possible to efficiently loop over every pixel


	// check to see how much the total flux will change by
/*
	max_pixel_change = 0.0;
	min_pixel_change = 0.0;
	for(ctr=0; ctr < imsize2; ctr++)
	{
		step_sum += max( model[0][ctr] + step_length * step[0][ctr] , min_flux );
		max_pixel_change = max(max_pixel_change , step[0][ctr]);
		min_pixel_change = min(min_pixel_change , step[0][ctr]);
	}
	
	
		// scale if necessary
	factor = fabs(step_sum - total_flux) / max_flux_change;
	cout<<"fabs(step_sum - total_flux) = "<<fabs(step_sum - total_flux)<<endl;
	cout<<"max_flux_change = "<<max_flux_change<<endl;
	cout<<"max_pixel_change = "<<max_pixel_change<<endl;
	cout<<"min_pixel_change = "<<min_pixel_change<<endl;
	cout<<"Factor = "<<factor<<endl;
	
	if( factor > 1.0)
	{
		step_length /= factor;
	}
*/
	#pragma omp parallel for  collapse(2) private( ctr, pnew, k, inew, factor)
	for(i=ignore_pixels;i<right_pixel_limit;i++)
	{
		for(j=ignore_pixels;j<right_pixel_limit;j++)
		{
			ctr = i * imsize + j;
	/*
			First take the Stokes I step
			Record the old value
			Set istep = step (which is the step length times the recommended step). If the step is negative, make sure it's not too negative
			Update Stokes I value with new flux, or min_flux - whichever is greater
	*/
			step[0][ctr] = max( model[0][ctr] + step_length * step[0][ctr] , min_flux );
			inew = step[0][ctr];

		//		Polarisation section

			if(npol > 1 )
			{
				// Find P, the total polarised flux, for the old and new cases

				pnew = 0.0;

				for(k=1;k<npol;k++)
				{
					pnew += pow( model[k][ctr] + step_length * step[k][ctr] , 2);
				}

				pnew = sqrt(pnew);

				// Set a factor to ensure that m, the fractional polarisation, is an element of [0,1]

				if( pnew < inew )
				{
					factor = 1.0;
				}
				else
				{
					factor = 0.8 * inew / pnew;
				}

				// Update the Q and U maps with their new values

				for(k=1;k<npol;k++)
				{
					step[k][ctr] = factor * ( model[k][ctr] + step_length * step[k][ctr] );
				}
			}
		}
	}

	return(0);
}

/*
	check_step

	check_step checks to see if the previous step length was near optimal

	inputs:
		old_model = the model Stokes I, Q and U maps from before the last step
		new_model = the model Stokes I, Q and U maps from after the last step
		new_residuals = the residual Stokes I, Q and U maps from after the last step
		default_map2 = the square of the default map
		imsize2 = the number of pixels in the new map
		alpha = the Lagrangian parameter associated with minimising Stokes I residuals
		beta = the Lagrangian parameter associated with minimising Stokes Q and U residuals
		gamma = the Lagrangian parameter associated with conserving the Stokes I flux
		q = the number of pixels in a beam
		npol = the number of polarisations being deconvolved

	outputs:
		J1 =  a measure of the change in J due to the step taken. To be compared with J0 to see if the step should have been bigger or smaller.
		on return, 0 if no error, 1 if error detected
*/

double check_step(double** old_model , double** new_model , double** new_residual, double* mask, double* default_map2 , double alpha , double beta , double gamma , int imsize, int ignore_pixels , int npol, double q)
{
	// Declare some variables

	double p , plog , m , dh , gradJ , step;
	int i, j, k , ctr;
	int imsize2 = imsize * imsize;
	int right_pixel_limit = imsize - ignore_pixels;

	double J1 = 0.0;

	// Start openmp for loop if possible to efficiently loop over every pixel

//	#pragma omp parallel for collapse(2) reduction(+:J1_temp) private(ctr,k,m,p,plog,dh,gradJ,step)
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
					p += (new_model[k][ctr] * new_model[k][ctr]);
				}
				p = sqrt(p);

			/*
					Need to implement dH and gradJ according to the MEM formulae
					To do this efficiently, make use of plog and m to prevent calculating the same thing more than once
					If the total polarised flux is sufficiently small, use the Taylor expansion of the logs to save computation
			*/

				if(p > 0.01 * new_model[0][ctr])
				{
					plog = 0.5 * log( ( new_model[0][ctr] - p ) / ( new_model[0][ctr] + p) ) / p;
				}
				else
				{
					m = p / new_model[0][ctr];
					plog = - (1.0 + (m * m / 3.0) ) / new_model[0][ctr];
				}

				for(k=0;k<npol;k++)
				{
					if(k==0)
					{
						dh = -0.5 * log( (new_model[0][ctr] * new_model[0][ctr] - p * p ) / (default_map2[ctr]));
						gradJ = dh - 2.0 * alpha * new_residual[0][ctr] - gamma;
					}
					else
					{
						dh = new_model[k][ctr] * plog;
						gradJ = dh - 2.0 * beta * new_residual[k][ctr];
					}

					step = new_model[k][ctr] - old_model[k][ctr];	// Step = difference between the two maps
					J1 += gradJ * step;			// The change in J due to this is gradJ
				}
			}
		}
	}

	return(J1);
}

/*
	interpolate_models

	interpolate_models interpolates current models and new models if necessary.
	It also enforces a minimum flux and makes sure the polarisation doesn't do anything strange

	inputs:

		current_model = the current model Stokes I, Q, U maps
		new_model = the new model Stokes I, Q, U maps
		frac_new = the weighting factor of the new Stokes I, Q, U maps
		imsize2 = the number of pixels in a map
		npol = the number of polarisations being deconvolved

	outputs:
		current_model = the interpolated model maps
		on return, 0
*/

int interpolate_models(double** current_model , double** new_model , double frac_new , int imsize, int ignore_pixels , int npol, double min_flux)
{
	double frac_old , inew , pnew , factor;
	int i , j, k, ctr;
	int imsize2 = imsize * imsize;
	int right_pixel_limit = imsize - ignore_pixels;

	frac_old = 1 - frac_new;

	#pragma omp parallel for collapse(2) private( ctr , inew , pnew , factor , k)
	for(i=ignore_pixels;i<right_pixel_limit;i++)
	{
		for(j=ignore_pixels;j<right_pixel_limit;j++)
		{
			ctr = i * imsize + j;
			// First interpolate the Stokes I part

			current_model[0][ctr] = max( frac_old * current_model[0][ctr] + frac_new * new_model[0][ctr] , min_flux);
			inew = current_model[0][ctr];

			// Now interpolate Stokes Q and U, making sure that the resulting m is an element of [0,1]

			if(npol > 1)
			{
				pnew = 0.0;
				for(k=1;k<npol;k++)
				{
					current_model[k][ctr] = frac_old * current_model[k][ctr] + frac_new * new_model[k][ctr];
					pnew += current_model[k][ctr] * current_model[k][ctr];
				}
				pnew = sqrt(pnew);

				// Set a factor to ensure that m, the fractional polarisation, is an element of [0,1]

				if( pnew < inew )
				{
					factor = 1.0;
				}
				else
				{
					factor = 0.8 * inew / pnew;
				}

				// apply interpolation

				if(factor < 1.0)
				{
					for(k=1;k<npol;k++)	
					{
							current_model[k][ctr] = factor * current_model[k][ctr];
					}
				}
			}
		}
	}

	return(0);
}

/*
	interpolate_residuals

	interpolate_residuals interpolates current resdiuals and new resdiuals if necessary.

	inputs:

		current_resdiuals = the current model Stokes I, Q, U maps
		new_resdiuals = the new model Stokes I, Q, U maps
		frac_new = the weighting factor of the new Stokes I, Q, U maps
		imsize2 = the number of pixels in a map
		npol = the number of polarisations being deconvolved

	outputs:
		current_resdiuals = the interpolated resdiual maps
		on return, 0
*/

int interpolate_residuals(double** current_residuals , double** new_residuals, double* chi2_rms , double frac_new , int imsize2 , int npol)
{
	double frac_old;

	frac_old = 1.0 - frac_new;

	for(int j=0;j<npol;j++)
	{
		chi2_rms[j] = 0.0;
		for(int i=0;i<imsize2;i++)
		{
			current_residuals[j][i] = frac_old * current_residuals[j][i] + frac_new * new_residuals[j][i];
			chi2_rms[j]+=(current_residuals[j][i]*current_residuals[j][i]);
		}
	}

	return(0);
}

/*
	copy_model

	copy_model updates either the model or residual maps by just copying the new maps into the old
	All it actually does is swap the pointers
	model2 points at model1's old memory afterwards

	inputs:
		model1 = a pointer to the map to be updated
		model2 = a pointer to the map to be copied

	outputs:
		model1 = model2
		model2 = model1
*/

int copy_model(double**& model1, double**& model2)
{
	double** temp;

	temp = model1;
	model1 = model2;
	model2 = temp;	// this actually swaps model 1 and model 2

	return(0);
}
