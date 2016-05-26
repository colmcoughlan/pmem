#include "pmem.hpp"
	using namespace std;	

int newton_raphson(double**& current_model, double**& new_model, double**& current_residuals, double**& new_residuals, double* default_map2, double* mask, double** convolved_model, fftw_complex* dirty_beam_ft, fftw_complex* complex_buff, double* double_buff, int pad_factor, fftw_plan& forward_transform, fftw_plan& backward_transform, double** dirty_map, gradient_structure& grad, double& alpha, double& beta, double& gamma, double& alpha_old, double& beta_old, double& gamma_old, bool force_chi2_method, double zsf, int imsize, int ignore_edge_pixels, double& imin, double& imax, double min_flux, int npol, double& q, double pol_upweight_factor, double* chi2_rms, double* rms_theoretical, double* current_rms, int* noise_box, double acceleration_factor, double& old_step_length1, double& old_rescale_step, bool conserve_flux, bool debug, int ctr)
{
//	printf("Inner Address of x is %p\n", (void *)current_model);  

	double delta_E, delta_F, delta_G;
	double J0, J1;
	double step_length1, step_limit, rescale_step;
	double total_flux, temp;
	int i, k, err;

	err = get_info(  current_model , current_residuals, mask , default_map2 , grad , total_flux , alpha , beta , gamma , imin, imax , imsize, ignore_edge_pixels , npol , q);	// get grads etc.
	if(err!=0)
	{
		cout<<endl<<"Error detected in get_info, err = "<<err<<endl<<endl;
		return(1);
	}

	k = (imsize - 2.0 * ignore_edge_pixels) * (imsize - 2.0* ignore_edge_pixels);	// number of pixels being used

	delta_E = (chi2_rms[0] - (rms_theoretical[0] * rms_theoretical[0] * k) ) / q;	// find differences in E, F and G

	delta_F = 0.0;
	for(i=1;i<npol;i++)
	{
		delta_F += ( ( chi2_rms[i] - (rms_theoretical[i] * rms_theoretical[i] * k) ) / q);
	}
	delta_F *= pol_upweight_factor;	// upweight polarisation error

	delta_G =  total_flux - zsf;
	
	alpha_old = alpha;
	beta_old = beta;
	gamma_old = gamma;


	err = new_ABG( grad , delta_E , delta_F , delta_G , alpha , beta , gamma , conserve_flux , npol , force_chi2_method);	// update values of alpha, beta, gamma
	if(err!=0)
	{
		cout<<endl<<"Error detected in new_ABG, err = "<<err<<endl<<endl;
		return(1);
	}
	
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
	
	return(0);
}