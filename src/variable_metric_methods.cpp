#include "pmem.hpp"
	using namespace std;
	
	
// The Newton raphson method uses the diagonal of the exact Hessian to calculate the metric

int newton_raphson(double**& current_model, double**& new_model, double**& current_residuals, double**& new_residuals, double* default_map2, double* mask, double** convolved_model, fftw_complex* dirty_beam_ft, fftw_complex* complex_buff, double* double_buff, int pad_factor, fftw_plan& forward_transform, fftw_plan& backward_transform, double** dirty_map, gradient_structure& grad, double& alpha, double& beta, double& gamma, double& alpha_old, double& beta_old, double& gamma_old, bool force_chi2_method, double zsf, int imsize, int ignore_edge_pixels, double& imin, double& imax, double min_flux, int npol, double& q, double pol_upweight_factor, double* chi2_rms, double* rms_theoretical, double* current_rms, int* noise_box, double acceleration_factor, double& old_step_length1, double& old_rescale_step, bool conserve_flux, bool debug, int ctr)
{
//	printf("Inner Address of x is %p\n", (void *)current_model);  

	double delta_E, delta_F, delta_G;
	double J0, J1;
	double step_length1, step_limit, rescale_step;
	double total_flux, temp;
	int i, k, err;
	
	
	// find the current metric (no need to save it, and find all the useful norms for the CE algorithm)

	err = get_info(  current_model , current_residuals, mask , default_map2 , grad , total_flux , alpha , beta , gamma , imin, imax , imsize, ignore_edge_pixels , npol , q, new_model, false);	// get grads etc.
	if(err!=0)
	{
		cout<<endl<<"Error detected in get_info, err = "<<err<<endl<<endl;
		return(1);
	}

	delta_E = ( current_rms[0] - rms_theoretical[0]);	// Alpha should "encourage" chi2_I to have rms_theorertical (dalpha = 0 at Jmax)
	delta_F = 0.0;
	for(k=1;k<npol;k++)
	{
		delta_F += pow(current_rms[k] - rms_theoretical[k], 2.0);
	}
	delta_F = sqrt( delta_F);
	delta_F *= pol_upweight_factor;	// upweight polarisation error
	
	delta_G = (total_flux - zsf);	// delta G is just the flux offset. gamma slows down changes when the flux is off
	
	alpha_old = alpha;
	beta_old = beta;
	gamma_old = gamma;
	
	// update ABG using the CE formulae


	err = new_ABG( grad , delta_E , delta_F , delta_G , alpha , beta , gamma , conserve_flux , npol , force_chi2_method);	// update values of alpha, beta, gamma
	if(err!=0)
	{
		cout<<endl<<"Error detected in new_ABG, err = "<<err<<endl<<endl;
		return(1);
	}
	
	
	// add some magic that helps stability when close to a solution
	// This restricts changes in Lagrangian parameters to inside a small interval


	alpha = max( alpha , max(beta,gamma));
	if(rms_theoretical[0] !=0.0 and ctr > 3 )
	{
		if( current_rms[0]/rms_theoretical[0] < 5.0)
		{
			if(debug)
			{
				cout<<"Restricting changes in Lagrangian parameters"<<endl;
			}
			force_chi2_method = true;
			alpha = max (alpha , alpha_old);
			beta = max( beta, beta_old);
			gamma = max( gamma, gamma_old);
		}
		else
		{
			force_chi2_method = false;
		}
	}

	if(debug || (ctr%100 == 0) )
	{
		cout<<"Iteration number "<<ctr<<endl;
		cout<<"Total flux, max and min = "<<total_flux<<" , "<<imax<<" , "<<imin<<endl;
	}
	
	// find the best step from the current gradient + metric

	err = cal_step( current_model , current_residuals, mask , default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol , q , J0 , new_model);	// find a good step
	if(err!=0)
	{
		cout<<endl<<"Error detected in cal_step, err = "<<err<<endl<<endl;
		return(1);
	}

	// make sure it's not too big

	step_limit = 1.0;
	if( grad.JJ > 0 )
	{
		step_limit = min( 2.0 , acceleration_factor * 0.15 * grad.II / grad.JJ );	// this line is very very very important... (especially the 0.15 factor)
	}
	step_length1 = min( 0.5 * (1.0 + old_step_length1) , step_limit);	// step length etc. is measured as a fraction (between 0 and 1)
	old_step_length1 = step_length1;
	J0 *= step_length1;
	
	// take it (with a scaling factor)

	err = take_step( current_model , new_model , step_length1 , step_limit , imsize, ignore_edge_pixels , npol, total_flux, max(fabs(total_flux - zsf),0.1), zsf, min_flux);	// take the step calculated, but scaled by step_length1
	if(err!=0)
	{
		cout<<endl<<"Error detected in take_step, err = "<<err<<endl<<endl;
		return(1);
	}
	
	// Get a better estimate of J0

	temp = find_j0( current_model , current_residuals, mask , default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol , q , new_model, current_model);
//		cout<<"Old J0 = "<<J0<<".Fancy J0 = "<<temp<<endl;
	J0 = temp;

	// convolve new model maps with dirty beam and get new residual maps

	for(i=0;i<npol;i++)
	{
		convolve( new_model[i] , dirty_beam_ft , imsize , pad_factor  , convolved_model[i] , forward_transform , backward_transform , double_buff , complex_buff);
		chi2_rms[i] = get_residual_map( dirty_map[i] , convolved_model[i] , new_residuals[i] , imsize, ignore_edge_pixels );
	}
	
	// See if the step actually helped (ideally gradJ at the new position (gradJnew called J1) should be zero)

	J1 = check_step( current_model , new_model , new_residuals, mask , default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol , q);	// check if the step was near optimal

	if( J0 - J1 != 0.0 )//	JX = gradJ at point X. Need gradJ = 0; look at J0 and J1 and try predict a rescaling factor that will result in gradJ=0
	{
		rescale_step = J0 / (J0 - J1);	// Derive this formula by finding the slope of J0 to J1 w.r.t. step, and using it to predict the step change which will set gradJ = JN = 0
	}
	else
	{
		rescale_step = 1.0;
	}
	rescale_step = 0.5 * (rescale_step + old_rescale_step);
	rescale_step = min( rescale_step , step_limit / step_length1 );
	old_rescale_step = rescale_step;

	if(debug)
	{
		cout<<"Initial step, Rescaling factor, step limit = "<<step_length1<<" , "<<rescale_step<<" , "<<step_limit<<endl;	// output some info
	}

	// line minimisation based on the value of rescale_step

	if( fabs( rescale_step - 1.0) > 0.05 )	// if step 1 was okay, just use it, but if step 2 offers a decent advantage take the average of step1 and step2
	{
		err = interpolate_models( current_model , new_model , rescale_step , imsize, ignore_edge_pixels , npol, min_flux);	// interpolate between old and new models, scaling with step2
		if(err!=0)
		{
			cout<<endl<<"Error detected in interpolation of new and old models, err = "<<err<<endl<<endl;
			return(1);
		}

		err = interpolate_residuals( current_residuals , new_residuals, chi2_rms , rescale_step , imsize*imsize , npol);	// interpolate residuals too
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
	
	// change q if things are slow
	
	if(fabs(step_length1-1.0) < 0.05)
	{
		q *= sqrt((1.0/max(0.5,min(2.0,step_length1*rescale_step))+3.0)/4.0);
//			cout<<"Q updated to "<<q<<endl;
	}

	if(debug)
	{
		cout<<"Alpha, beta, gamma =  "<<alpha<<" , "<<beta<<" , "<<gamma<<endl;
		cout<<"GradJ.J, Grad1.1, J0, J1 = "<<grad.JJ<<" , "<<grad.II<<" , "<<J0<<" , "<<J1<<endl;
		cout<<"Delta E, F, G = "<<delta_E<<" , "<<delta_F<<" , "<<delta_G<<endl;
		cout<<"GradE.E, GradF.F, GradG.G = "<<grad.EE<<" , "<<grad.FF<<" , "<<grad.GG<<endl;
		cout<<"GradE.F, GradF.G, GradE.G = "<<grad.EF<<" , "<<grad.FG<<" , "<<grad.EG<<endl;	// output even more info
		cout<<"GradE.H, GradF.H, GradG.H = "<<grad.EH<<" , "<<grad.FH<<" , "<<grad.GH<<endl;
		cout<<"GradE.J, GradF.J, GradG.J = "<<grad.EJ<<" , "<<grad.FJ<<" , "<<grad.GJ<<endl;
	}





	// update the rms that has been reached

	for(i=0;i<npol;i++)
	{
		current_rms[i] = rms_region( current_residuals[i] , noise_box[0] , noise_box[1] , noise_box[2] , noise_box[3] , imsize);	// this might actually go up at the start, as the residual map is not a noise map initially
	}
	
	return(0);
}

// note : using ddot is definitely much faster than single threaded looping! (and that even with the double** structure)
// Collapsed multithreaded looping is about the same speed - maybe marginally slower
// Use blas and let someone else worry about implementing the most efficient multithreaded solution!

double dot_product(double** vec1, double** vec2, int imsize2, int npol)
{
	double sum = 0.0;
	int inc = 1;

	for(int k=0;k<npol;k++)
	{
		sum += ddot_(&imsize2, vec1[k], &inc, vec2[k], &inc);
	}

	return(sum);
}

// there is probably a BLAS function appropriate for this somewhere, but based on the results above, not much more efficient
// Return vec1*vec2*metric -> assuming a diagonal metric

double metric_product(double** vec1, double** vec2, double** matrix, int imsize2, int npol)
{
	double sum = 0.0;

	#pragma omp parallel for collapse(2) reduction( +: sum)
	for(int k=0;k<npol;k++)
	{
		for(int i=0;i<imsize2;i++)
		{
			sum += vec1[k][i] * vec2[k][i] * matrix[k][i];
		}
	}
	
	return(sum);
}

// subtract two vectors ans = vec1 - vec2, storing result in vec2

int sub_vec(double** vec1, double** vec2, int imsize2, int npol)
{
	#pragma omp parallel for collapse(2)
	for(int k=0;k<npol;k++)
	{
		for(int i=0;i<imsize2;i++)
		{
			vec2[k][i] = vec1[k][i] - vec2[k][i];
		}
	}
	
	return(0);
}

// update the hessian according to the DFP or BFGS formulae
// store the new Hessian, but also update relevant norms for the CE algorithm

int update_hessian(double** model , double** residual, double* mask, double* default_map2 , gradient_structure& grad, 
double& total_flux , double alpha , double beta , double gamma , double& imin, double& imax, double** model_step, 
double** grad_step, double** hessian, int imsize, int ignore_edge_pixels , int npol, bool do_bfgs)
{
	int right_pixel_limit = imsize - ignore_edge_pixels;
	int i,j,k, ctr;
	int imsize2 = imsize* imsize;
	double t1, t2, t3, u, dx_dot_df, df_dot_hessian_dot_df;
	double p, plog, pexp, m, dh, de, df, dg, metric, temp;
	double EE , EF , EG , EH , FF , FG , FH , GG , GH , HH , II;
	
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
	
	dx_dot_df = dot_product(model_step, grad_step, imsize2, npol);
	df_dot_hessian_dot_df = metric_product(grad_step, grad_step, hessian, imsize2, npol);
	
//	cout<<"dx_dot_df = "<<dx_dot_df<<endl;
//	cout<<"df_dot_hessian_dot_df = "<<df_dot_hessian_dot_df<<endl;
	
	if(dx_dot_df * df_dot_hessian_dot_df == 0.0)	// make sure these are non-zero
	{
		return(1);
	}

	t3 = 0.0;	// only changes if solving via BFGS
	#pragma omp parallel for collapse(3) private(k, i, j, ctr, t1, t2, t3, u)
	for(k=0;k<npol;k++)	// note only the diagonal terms are being kept!
	{
		for(i=ignore_edge_pixels;i<right_pixel_limit;i++)
		{
			for(j=ignore_edge_pixels;j<right_pixel_limit;j++)	// See numerical recipes for description of process
			{
				ctr = i * imsize + j;
				t1 = (model_step[k][ctr]*model_step[k][ctr])/dx_dot_df;
				t2 = - pow(hessian[k][ctr] * grad_step[k][ctr], 2.0) / df_dot_hessian_dot_df;
				if(do_bfgs)	// 3rd term needed if doing BFGS
				{
					u = (model_step[k][ctr]/dx_dot_df) - (hessian[k][ctr] * grad_step[k][ctr])/df_dot_hessian_dot_df;
					t3 = df_dot_hessian_dot_df * u * u;
				}
				hessian[k][ctr] = hessian[k][ctr] + t1 + t2 +t3;
			}
		}
	}
	
	
				
	// now use hessian to find new norms etc.
	
	#pragma omp parallel for collapse(2) reduction( +: EE , EG , EH , FF , FH , GG , GH , HH ,temp) private(ctr,k,p,m,plog,pexp,dh,de,df,dg,metric)
	for(i=ignore_edge_pixels;i<right_pixel_limit;i++)
	{
		for(j=ignore_edge_pixels;j<right_pixel_limit;j++)
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
					metric = -1.0/hessian[k][ctr];	// remember that the metric is -1.0/hessian
					if(k==0)	// Stokes I only
					{
						dh = -0.5 * log( ( model[0][ctr] * model[0][ctr] - p * p ) / ( default_map2[ctr] ) );	// derivative of the entropy contribution of pixel
						de = 2.0 * residual[0][ctr];								// derivative of the Stokes I chi2 contribution of the pixel
						df = 0.0;										// derivative of the Stokes Q,U,V chi2 contribution of the pixel
						dg = 1.0;										// derivative of the flux conservation term contribution of the pixel
					}
					else	// Other Stokes parameters (same idea)
					{
						dh = model[k][ctr] * plog;
						de = 0.0;
						df = 2.0 * residual[k][ctr];
						dg = 0.0;
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
	
	return(0);
}


// Quasi newton methods are very like the newton method, but make up the Hessian instead of calculating it directly
// Note the Hessian is still diagonal
// It is initialised to the exact Newtonian diagonal Hessian

int quasi_newton(double**& current_model, double**& new_model, double**& current_residuals, double**& new_residuals, 
double* default_map2, double* mask, double** convolved_model, fftw_complex* dirty_beam_ft, fftw_complex* complex_buff, 
double* double_buff, int pad_factor, fftw_plan& forward_transform, fftw_plan& backward_transform, double** dirty_map, 
gradient_structure& grad, double& alpha, double& beta, double& gamma, double& alpha_old, double& beta_old, double& gamma_old, 
bool force_chi2_method, double zsf, int imsize, int ignore_edge_pixels, double& imin, double& imax, double min_flux, int npol, 
double& q, double pol_upweight_factor, double* chi2_rms, double* rms_theoretical, double* current_rms, int* noise_box, 
double acceleration_factor, double& old_step_length1, double& old_rescale_step, bool conserve_flux, bool debug, int ctr, 
double** gradold, double** gradnew, double** hessian, bool do_bfgs)
{
	double delta_E, delta_F, delta_G;
	double J0, J1;
	double step_length1, step_limit, rescale_step;
	double total_flux, temp;
	int i, k, err;
	int imsize2 = imsize * imsize;
	
	
	// find current gradient
	
	err = find_gradj(current_model , current_residuals, mask, default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol, gradnew);

	
	if(ctr == 1)	// on the first iteration we need to find an estimate of the Hessian (Newton-Raphson)
	{
		err = get_info(  current_model , current_residuals, mask , default_map2 , grad , total_flux , alpha , beta , gamma , imin, 
				imax , imsize, ignore_edge_pixels , npol , q, hessian, true);	// get grads etc.
		if(err!=0)
		{
			cout<<endl<<"Error detected in get_info, err = "<<err<<endl<<endl;
			return(1);
		}
	}
	else
	{
		#pragma omp parallel for collapse(2)	// Otherwise use DFP or BFGS to update current Hessian
		for(k=0;k<npol;k++)
		{
			for(i=0;i<imsize2;i++)	// need to find recent steps in model and gradJ
			{
				new_model[k][i] = current_model[k][i] - new_model[k][i]; //new model set to model step: x_{i+1} - x_{i}
			}
		}
		err = sub_vec(gradnew, gradold, imsize2, npol);	// gradold is now set to gradient step
		err = update_hessian( current_model , current_residuals, mask , default_map2 , grad , total_flux , alpha , beta , gamma , imin, imax, 
					new_model, gradold, hessian, imsize, ignore_edge_pixels , npol, do_bfgs);
		if(err!=0)
		{
			cout<<"Error updating Hessian"<<endl;
			return(1);
		}
	}

	delta_E = ( current_rms[0] - rms_theoretical[0]);	// Alpha should "encourage" chi2_I to have rms_theorertical (dalpha = 0 at Jmax)
	delta_F = 0.0;
	for(k=1;k<npol;k++)
	{
		delta_F += pow(current_rms[k] - rms_theoretical[k], 2.0);
	}
	delta_F = sqrt( delta_F);
	delta_F *= pol_upweight_factor;	// upweight polarisation error
	
	delta_G = (total_flux - zsf);	// delta G is just the flux offset. gamma slows down changes when the flux is off
	
	alpha_old = alpha;
	beta_old = beta;
	gamma_old = gamma;
	
	
	// update ABG (CE formulae)


	err = new_ABG( grad , delta_E , delta_F , delta_G , alpha , beta , gamma , conserve_flux , npol , force_chi2_method);	// update values of alpha, beta, gamma
	if(err!=0)
	{
		cout<<endl<<"Error detected in new_ABG, err = "<<err<<endl<<endl;
		return(1);
	}
	
	// magic to help convergence
	// This restricts changes in Lagrangian parameters to inside a small interval

	alpha = max( alpha , max(beta,gamma));
	if(rms_theoretical[0] !=0.0 and ctr > 3 )
	{
		if( current_rms[0]/rms_theoretical[0] < 5.0)
		{
			if(debug)
			{
				cout<<"Restricting changes in Lagrangian parameters"<<endl;
			}
			force_chi2_method = true;
			alpha = max (alpha , alpha_old);
			beta = max( beta, beta_old);
			gamma = max( gamma, gamma_old);
		}
		else
		{
			force_chi2_method = false;
		}
	}

	if(debug || (ctr%100 == 0) )
	{
		cout<<"Iteration number "<<ctr<<endl;
		cout<<"Total flux, max and min = "<<total_flux<<" , "<<imax<<" , "<<imin<<endl;
	}
	
	
	//
	

	err = cal_step( current_model , current_residuals, mask , default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol , q , J0 , new_model);	// find a good step
	if(err!=0)
	{
		cout<<endl<<"Error detected in cal_step, err = "<<err<<endl<<endl;
		return(1);
	}

	step_limit = 1.0;
	if( grad.JJ > 0 )
	{
		step_limit = min( 2.0 , acceleration_factor * 0.15 * grad.II / grad.JJ );	// this line is very very very important... (especially the 0.15 factor)
	}
	step_length1 = min( 0.5 * (1.0 + old_step_length1) , step_limit);	// step length etc. is measured as a fraction (between 0 and 1)
	old_step_length1 = step_length1;
	J0 *= step_length1;

	err = take_step( current_model , new_model , step_length1 , step_limit , imsize, ignore_edge_pixels , npol, total_flux, max(fabs(total_flux - zsf),0.1), zsf, min_flux);	// take the step calculated, but scaled by step_length1
	if(err!=0)
	{
		cout<<endl<<"Error detected in take_step, err = "<<err<<endl<<endl;
		return(1);
	}

	temp = find_j0( current_model , current_residuals, mask , default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol , q , new_model, current_model);
//		cout<<"Old J0 = "<<J0<<".Fancy J0 = "<<temp<<endl;
	J0 = temp;

	// convolve new model maps with dirty beam and get new residual maps

	for(i=0;i<npol;i++)
	{
		convolve( new_model[i] , dirty_beam_ft , imsize , pad_factor  , convolved_model[i] , forward_transform , backward_transform , double_buff , complex_buff);
		chi2_rms[i] = get_residual_map( dirty_map[i] , convolved_model[i] , new_residuals[i] , imsize, ignore_edge_pixels );
	}

	J1 = check_step( current_model , new_model , new_residuals, mask , default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol , q);	// check if the step was near optimal

	if( J0 - J1 != 0.0 )//	JX = gradJ at point X. Need gradJ = 0; look at J0 and J1 and try predict a rescaling factor that will result in gradJ=0
	{
		rescale_step = J0 / (J0 - J1);	// Derive this formula by finding the slope of J0 to J1 w.r.t. step, and using it to predict the step change which will set gradJ = JN = 0
	}
	else
	{
		rescale_step = 1.0;
	}
	rescale_step = 0.5 * (rescale_step + old_rescale_step);
	rescale_step = min( rescale_step , step_limit / step_length1 );
	old_rescale_step = rescale_step;

	if(debug)
	{
		cout<<"Initial step, Rescaling factor, step limit = "<<step_length1<<" , "<<rescale_step<<" , "<<step_limit<<endl;	// output some info
	}


	if( fabs( rescale_step - 1.0) > 0.05 )	// if step 1 was okay, just use it, but if step 2 offers a decent advantage take the average of step1 and step2
	{
		err = interpolate_models( current_model , new_model , rescale_step , imsize, ignore_edge_pixels , npol, min_flux);	// interpolate between old and new models, scaling with step2
		if(err!=0)
		{
			cout<<endl<<"Error detected in interpolation of new and old models, err = "<<err<<endl<<endl;
			return(1);
		}

		err = interpolate_residuals( current_residuals , new_residuals, chi2_rms , rescale_step , imsize*imsize , npol);	// interpolate residuals too
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
	
	if(fabs(step_length1-1.0) < 0.05)
	{
		q *= sqrt((1.0/max(0.5,min(2.0,step_length1*rescale_step))+3.0)/4.0);
//			cout<<"Q updated to "<<q<<endl;
	}

	if(debug)
	{
		cout<<"Alpha, beta, gamma =  "<<alpha<<" , "<<beta<<" , "<<gamma<<endl;
		cout<<"GradJ.J, Grad1.1, J0, J1 = "<<grad.JJ<<" , "<<grad.II<<" , "<<J0<<" , "<<J1<<endl;
		cout<<"Delta E, F, G = "<<delta_E<<" , "<<delta_F<<" , "<<delta_G<<endl;
		cout<<"GradE.E, GradF.F, GradG.G = "<<grad.EE<<" , "<<grad.FF<<" , "<<grad.GG<<endl;
		cout<<"GradE.F, GradF.G, GradE.G = "<<grad.EF<<" , "<<grad.FG<<" , "<<grad.EG<<endl;	// output even more info
		cout<<"GradE.H, GradF.H, GradG.H = "<<grad.EH<<" , "<<grad.FH<<" , "<<grad.GH<<endl;
		cout<<"GradE.J, GradF.J, GradG.J = "<<grad.EJ<<" , "<<grad.FJ<<" , "<<grad.GJ<<endl;
	}


	// update the rms that has been reached

	for(i=0;i<npol;i++)
	{
		current_rms[i] = rms_region( current_residuals[i] , noise_box[0] , noise_box[1] , noise_box[2] , noise_box[3] , imsize);	// this might actually go up at the start, as the residual map is not a noise map initially
	}
	
	// update the step taken, and the change in gradient
	
	err = copy_model( gradnew , gradold );	// gradold is now the latest grad
	
	return(0);
}