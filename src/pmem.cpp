// Version 1.01.

/*
 Copyright (c) 2014, Colm Coughlan
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


	Polarised Maxmimum Entropy Method
	Colm Coughlan - colmcoughlanirl   %%%at%%% gmail.com
*/


#include "pmem.hpp"

	using namespace std;


int write_csv(string filename, double* array, int imsize);
string int2str (int num);
int clip_edges(double* map, double replacement_value, int edge_limit, int imsize);
int zero_array(double* array, int imsize);

int main()
{
	cout<<"PMEM version 1.01"<<endl;

	fitsinfo_map fitsi[4];

	char history[] = "MEM deconvolution performed by PMEM. See PMEM logs for details.";

	string line;
	string filename_dirty_map[4];
	string filename_dirty_beam;
	string output_name;
	string filename_default_map;
	string filename_mask;
	gradient_structure grad;	// room to transport gradient information
	fstream fout;

	bool converged;
	bool converged_temp;
	bool conserve_flux;
	bool estimate_flux;
	bool debug;

	double current_rms[4];	// current Stokes I, Q, U , V residual rms
	double min_rms[] = {1e99,1e99,1e99,1e99};	// minimum Stokes I residual rms achieved (should be updated almost every iteration)
	int min_rms_i_ctr;	// counter value at which the minimum Stokes I residual rms has been achieved
	int convergence_limit;	// the number of allowed iterations where the Stokes I residuals are not decreasing before the program quits

	int i,j , k,ctr;
	int err;
	int imsize;
	int imsize2;
	int npol;		// no. of polarisations being deconvolved
	int niter;		// max no. of iterations
	int pad_factor = 2;	// 1 for no padding, 2 for zero padding
	double restoring_beam[3];	// [bmaj, bmin, bpa]
	int noise_box[4];	// [blcx, blcy, trcx, trcy]
	int ignore_edge_pixels;	// number of pixels at the edge of the map to ignore (useful for aliasing problems)

	double* null_double;
	int peak[2];

	double acceleration_factor;
	double q_factor;
	double pol_upweight_factor;
	double temp,temp2;
	double pixels_per_beam;
	double total_flux;
	double alpha , beta , gamma;
	double imin , imax;
	double step_length1 , rescale_step;
	double old_step_length1 , old_rescale_step;
	double step_limit;
	double delta_E, delta_F, delta_G;
	double J0 , J1;
	double rms_tolerance , flux_tolerance , convergence_tolerance;
	double alpha_old, beta_old, gamma_old;
	bool force_chi2_method = false;
	bool have_mask;
	int imsize2_unmasked;



	double** dirty_map;
	double** current_model;
	double** new_model;
	double** current_residuals;
	double** new_residuals;
	double** convolved_model;

	double* dirty_beam;
	fftw_complex* dirty_beam_ft;
	fftw_complex* complex_buff;
	double* double_buff;
	double* default_map2;
	double* mask;
	double zsf;
	double* chi2_rms;
	double rms_theoretical[4];
	int outnum=0;

	double q = 0.0;
	double sum, oldsum;
	const double min_flux = 1e-15;	// set the minimum flux allowable for the Stokes I models (must be > 0 for model)
	
	rms_tolerance=1.1;	// allow a tolerance of 10 % in guess for noise
	flux_tolerance=0.1;
	convergence_tolerance=0.1;
	
	cout<<"Reading driver"<<endl;
	err = read_driver("mempy_driver.dat", npol, filename_dirty_map, filename_dirty_beam, filename_default_map, filename_mask, zsf, conserve_flux
, rms_theoretical, niter, restoring_beam, noise_box, acceleration_factor, q_factor, pol_upweight_factor, output_name, ignore_edge_pixels, debug);
	if( err !=0 )
	{
		cout<<"Error reading from "<<"mempy_driver.dat"<<endl;
		cout<<"Exiting program..."<<endl;
		return(1);
	}
	else
	{
		cout<<"Read complete"<<endl;
	}


//	err = quickfits_read_map_header( filename_dirty_map[0].c_str() , &imsize , &cell , &ra , &dec , centre_shift , rotations , &freq , &freq_delta , &stokes[0] , object , observer , telescope , &equinox , date_obs , &bmaj , &bmin , &bpa , &ncc ,-1 );
	err = quickfits_read_map_header( filename_dirty_map[0].c_str() , &fitsi[0] );
	if( err != 0 )
	{
		cout<<"Error reading header from "<<filename_dirty_map[0]<<endl;
		cout<<"Exiting program..."<<endl;
		return(1);
	}

	imsize = fitsi[0].imsize_ra;
	imsize2= imsize * imsize;

	//check that inputs are valid
	
	noise_box[0] = ignore_edge_pixels;
	noise_box[1] = noise_box[0];
	noise_box[2] = imsize - ignore_edge_pixels;
	noise_box[3] = noise_box[2];
	
	for(i=0;i<4;i++)
	{
		if((noise_box[i] < 0) or (noise_box[i] >= imsize))
		{
			cout<<"Error - invalid noise box."<<endl;
			cout<<"Coordinate "<<noise_box[i]<<" is invalid"<<endl;
			cout<<"Exiting program..."<<endl;
			return(1);
		}
	}

	j = pad_factor * imsize * ( pad_factor * imsize / 2 + 1 );	// allowing for padding in buffers and reducing memory needed for r2c and c2r transform by taking Herm. conjugacy into account

	// Declare memory


	dirty_map=new double*[npol];
	current_model=new double*[npol];
	new_model=new double*[npol];
	current_residuals=new double*[npol];
	new_residuals=new double*[npol];
	convolved_model=new double*[npol];
	chi2_rms=new double[npol];


	dirty_beam=new double[imsize2];
	dirty_beam_ft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * j );	// saving memory with herm. conjugacy
	complex_buff = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * j );	// allowing room for padding in buffers
	double_buff = ( double* ) fftw_malloc( sizeof( double ) * pad_factor * pad_factor * imsize2 );
	default_map2=new double[imsize2];
	mask = new double[imsize2];


	for(i=0;i<npol;i++)
	{
		dirty_map[i]=new double[imsize2];
		current_model[i]=new double[imsize2];
		new_model[i]=new double[imsize2];
		current_residuals[i]=new double[imsize2];
		new_residuals[i]=new double[imsize2];
		convolved_model[i]=new double[imsize2];

		chi2_rms[i] = 0.0;
	}

	if(err!=0)
	{
		goto free_mem_exit;
	}

	// prepare FT plans
	
	// start up some of the fourier tranform settings and variables from FFTW


	fftw_init_threads();	// initialise fftw parallisation
	fftw_plan_with_nthreads(omp_get_max_threads());
	fftw_plan forward_transform, backward_transform;

	forward_transform = fftw_plan_dft_r2c_2d(imsize * pad_factor , imsize * pad_factor , double_buff , complex_buff , FFTW_MEASURE );	// optimise FFT

	backward_transform = fftw_plan_dft_c2r_2d(imsize * pad_factor , imsize * pad_factor , complex_buff , double_buff , FFTW_MEASURE );	// r2c is always a forward transform etc



	if( restoring_beam[0] == 0 || restoring_beam[1] ==0)	// if no restoring beam is given, use the beam from the dirty maps
	{
		restoring_beam[0] = fitsi[0].bmaj;
		restoring_beam[1] = fitsi[0].bmin;
		restoring_beam[2] = fitsi[0].bpa;	// storing in degrees
	}
	else
	{
		restoring_beam[0] /= 3600.0;
		restoring_beam[1] /= 3600.0;
	}


	// print out some information

	cout<<"Running PMEM."<<endl<<endl;
	cout<<"Target files:"<<endl;
	cout<<"\t Stokes I : "<<filename_dirty_map[0]<<" with an estimated flux of "<<zsf<<" Jy."<<endl;
	for( i = 1; i < npol ; i++ )
	{
		cout<<"\t Stokes "<<i+1<<" : "<<filename_dirty_map[i]<<"."<<endl;
	}
	cout<<"\t Dirty beam : "<<filename_dirty_beam<<"."<<endl;
	if( filename_default_map.length() > 1)
	{
		cout<<"\t Default map : "<<filename_default_map<<"."<<endl;
	}
	else
	{
		cout<<"\t No default map used."<<endl;
	}
	
	if( filename_mask.length() > 1)
	{
		cout<<"\t Using mask : "<<filename_mask<<"."<<endl;
		have_mask = true;
	}
	else
	{
		cout<<"\t No masking will be used."<<endl;
		have_mask = false;
	}
	
	if(conserve_flux)
	{
		cout<<"\t Flux will be conserved to the estimate of "<<zsf<<endl;
	}
	else
	{
		cout<<"\t Flux will not be conserved"<<endl;
	}
	cout<<endl<<"Restoring beam information:"<<endl;
	cout<<"\t BMAJ = "<<restoring_beam[0]*3600.0<<" as, BMIN = "<<restoring_beam[1]*3600.0<<" as, BPA = "<<restoring_beam[2]<<" deg."<<endl<<endl;
	cout<<"Map details :"<<endl;
	cout<<"\t Imsize = "<<imsize<<" pixels."<<endl;
	cout<<"\t Cellsize = "<<fitsi[0].cell_ra*(3600.0*1000.0)<<" mas."<<endl<<endl;
	cout<<"Running parameters :"<<endl;
	cout<<"\t Acceleration factor = "<<acceleration_factor<<endl;
	cout<<"\t Q factor = "<<q_factor<<endl;
	cout<<"\t Polarisation upweight factor = "<<pol_upweight_factor<<endl;
	cout<<endl<<endl;


	// Now read in all maps

	// turn off any clean component handling
	for(i=0;i<npol;i++)
	{
		fitsi[0].ncc = 0;
	}
	
	err = quickfits_read_map(filename_dirty_map[0].c_str() , fitsi[0], dirty_map[0] , null_double , null_double , null_double);
	if(err!=0)
	{
		cout<<endl<<"Error detected reading map from "<<filename_dirty_map[0]<<", err = "<<err<<endl<<endl;
		cout<<"Program closing"<<endl;
		goto free_mem_exit;
	}


	for(i=1;i<npol;i++)	// read in all the polarisation maps (if any)
	{
		err = quickfits_read_map_header( filename_dirty_map[i].c_str() , &fitsi[i]);	// note the only new piece of information here is the stokes value (read in temp to avoid overwriting bmaj etc.)
		if(err!=0)
		{
			cout<<endl<<"Error detected reading header from "<<filename_dirty_map[i]<<", err = "<<err<<endl<<endl;
			cout<<"Program closing"<<endl;
			goto free_mem_exit;
		}


		err = quickfits_read_map( filename_dirty_map[i].c_str(), fitsi[i] , dirty_map[i] , null_double , null_double , null_double);
		if(err!=0)
		{
			cout<<endl<<"Error detected reading map from "<<filename_dirty_map[i]<<", err = "<<err<<endl<<endl;
			cout<<"Program closing"<<endl;
			goto free_mem_exit;
		}
	}


	err = quickfits_read_map( filename_dirty_beam.c_str(), fitsi[0] , dirty_beam , null_double , null_double , null_double );
	if(err!=0)
	{
		cout<<endl<<"Error detected reading beam from "<<filename_dirty_beam<<", err = "<<err<<endl<<endl;
		cout<<"Program closing"<<endl;
		goto free_mem_exit;
	}

	// estimate number of pixels per beam

	pixels_per_beam = fitsi[0].bmaj * fitsi[0].bmin * M_PI / ( 4.0 * log(2) * fitsi[0].cell_ra * fitsi[0].cell_dec );
	if(pixels_per_beam <= 0 )
	{
		cout<<"Error in estimating beamize. Strictly positive number of pixels per beam required."<<endl;
		goto free_mem_exit;
	}
	cout<<"Pixels per beam = "<<pixels_per_beam<<endl;


	// initialize maps

	if( imsize - 2 * ignore_edge_pixels <= 0)
	{
		cout<<"Ignore edge pixels = "<<ignore_edge_pixels<<" is too large for image of size "<<imsize<<endl;
		goto free_mem_exit;
	}
	
	
	// Read in mask if present. If not, set entire image as unmasked. Binary mask assumed (0 = masked, > 0 unmasked)

	if( have_mask )
	{
		fitsi[0].ncc = 0;
		err = quickfits_read_map( filename_mask.c_str() , fitsi[0] , mask , null_double , null_double , null_double );	// load default map if given
		if(err!=0)
		{
			cout<<"Error loading "<<filename_mask<<endl;
			goto free_mem_exit;
		}
		
		imsize2_unmasked = 0;
		for(i=0;i<imsize2;i++)
		{
			if(mask[i] > 0.0)
			{
				imsize2_unmasked++;
			}
		}
		cout<<"Percentage of image unmasked : "<<(float(imsize2_unmasked)/float(imsize2))*100.0<<"%."<<endl;
		cout<<"Note mask should include any zero edges."<<endl;
	}
	else
	{
		for(i=0;i<imsize2;i++)
		{
			mask[i] = 1.0;
		}
		imsize2_unmasked = imsize2;
	}


	// initialize default map

	if( filename_default_map.length() > 1)
	{
		cout<<"Using "<<filename_default_map<<" as default map."<<endl;

		fitsi[0].ncc = 0;
		err = quickfits_read_map( filename_default_map.c_str() , fitsi[0] , default_map2 , null_double , null_double , null_double );	// load default map if given
		if(err!=0)
		{
			cout<<"Error loading "<<filename_default_map<<endl;
			goto free_mem_exit;
		}
	}
	else
	{
		err = zero_array(default_map2, imsize);
		
		if(have_mask)
		{
			temp = zsf / imsize2_unmasked;
			temp *=temp;
			alpha = min_flux*min_flux;
			for(i=0;i<imsize2;i++)
			{
				if(mask[i] > 0)
				{
					default_map2[i] = temp;
				}
				else
				{
					default_map2[i] = alpha;
				}
			}
		}
		else
		{
			ctr = (imsize - ignore_edge_pixels);
			temp = temp = zsf / float(pow(imsize - 2 * ignore_edge_pixels,2));
			cout<<"Making flat default map with constant value = "<<temp<<" Jy."<<endl;
			temp *= temp;
			for( i = ignore_edge_pixels ; i < ctr ; i++ )
			{
				for( j = ignore_edge_pixels ; j < ctr ; j++ )
				{
					default_map2[i*imsize + j] = temp;
				}
			}
		}
	}
	
	// COLM : At the moment this only works for the flat default bias map.
/*
	ctr = imsize - ignore_edge_pixels;
//	temp = zsf / float(pow(imsize - 2 * ignore_edge_pixels,2));
	temp = zsf / imsize2;
	temp *= temp;
	for(i = 0; i < imsize ; i++)
	{
		for(j = 0; j < imsize ; j++)
		{
			if( (i < ignore_edge_pixels or j < ignore_edge_pixels) or (i>ctr or j > ctr) )
			{
				default_map2[i*imsize + j] = 0;	// set edges of default map squared to min_flux squared
			}
			else
			{
				default_map2[i*imsize + j] = temp;
			}
		}
	}
*/
	


	// initialize Stokes I model to flat minimum flux. Everything else to zero.
	
	for(i=0; i<npol;i++)
	{
		err = zero_array(new_residuals[i], imsize);
		err = zero_array(current_residuals[i], imsize);
		err = zero_array(current_model[i], imsize);
		err = zero_array(new_model[i], imsize);
	}
	
	for(i=0;i<imsize2;i++)
	{
		current_model[0][i] = sqrt(default_map2[i]);
	}
/*
	ctr = (imsize - ignore_edge_pixels);
	for( i = ignore_edge_pixels ; i < ctr ; i++ )
	{
		for( j = ignore_edge_pixels ; j < ctr ; j++ )
		{
			k = i*imsize + j;
			if(mask[k] > 0)
			{
				current_model[0][k] = min_flux;
			}
			else
			{
				current_model[0][k] = min_flux;
			}
		}
	}
*/
	cout<<"Set to ignore "<<ignore_edge_pixels<<" edge pixels."<<endl;
	cout<<"Acceleration factor set to "<<acceleration_factor<<endl;

	// find ft of the dirty beam

	if( pad_factor == 1)
	{
		arrange_ft( dirty_beam , imsize );	// rearrange beam to prepare for use in convolution
	}

	ft_beam(dirty_beam , dirty_beam_ft , imsize , pad_factor , forward_transform , double_buff , complex_buff);	// get ft


	// convolve initial model maps with dirty beam and get residual maps

	for(i=0;i<npol;i++)
	{
		convolve( current_model[i] , dirty_beam_ft , imsize , pad_factor  , convolved_model[i] , forward_transform , backward_transform , double_buff , complex_buff);
		chi2_rms[i] = get_residual_map( dirty_map[i] , convolved_model[i] , current_residuals[i] ,  imsize, ignore_edge_pixels );
	}

	// Caculate Q, a factor used in the approximation of the Hessian matrix as diagonal. Important value. (see Cornwell & Evans 1984, Sault 1990)

	q = 0.0;
	#pragma omp parallel for reduction( +: q)
	for( i = 0 ; i < imsize2 ; i++ )
	{
		q += dirty_beam[i]*dirty_beam[i];
	}
	q = q_factor * sqrt( q );
	cout<<"Q = "<<q<<endl;


	for(k=0; k<npol; k++)
	{
		min_rms[k] = 999999.9999;
	}

	converged=false;
	ctr=0;
	old_step_length1=0.0;
	old_rescale_step=0.0;
	min_rms_i_ctr = 0;
	convergence_limit=200;
	alpha_old = 0.0;
	beta_old = 0.0;
	gamma_old = 0.0;
	
	alpha = 0.0;
	if(npol > 1)
	{
		beta = 0;
	}
	else
	{
		beta = 0.0;
	}
	if(conserve_flux)
	{
		gamma = 0.0;
	}
	else
	{
		gamma = 0.0;
	}





/*
#########################################################################################################################
#
#
#	// start iterations
#
#
#########################################################################################################################
*/



	
	

	while(!converged && ctr < niter)
	{
		ctr++;
		
		err = get_info(  current_model , current_residuals, mask , default_map2 , grad , total_flux , alpha , beta , gamma , imin, imax , imsize, ignore_edge_pixels , npol , q);	// get grads etc.
		if(err!=0)
		{
			cout<<endl<<"Error detected in get_info, err = "<<err<<endl<<endl;
			goto free_mem_exit;
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


//		force_chi2_method = true;

		err = new_ABG( grad , delta_E , delta_F , delta_G , alpha , beta , gamma , conserve_flux , npol , force_chi2_method);	// update values of alpha, beta, gamma
		if(err!=0)
		{
			cout<<endl<<"Error detected in new_ABG, err = "<<err<<endl<<endl;
			break;
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
			break;
		}
/*
		line.assign("StokesIstep.fits");
		line.append(int2str(outnum));
		err = quickfits_write_map( line.c_str() , new_model[0] , imsize , cell , ra , dec , centre_shift , rotations , freq , freq_delta , stokes[i] , object , observer , telescope , equinox , date_obs , history , restoring_beam[0] , restoring_beam[1] , restoring_beam[2] , ctr , false);
		if(outnum == 0)
		{
			outnum =1;
		}
		else
		{
			outnum=0;
		}
*/
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
			break;
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

		err = check_step( current_model , new_model , new_residuals, mask , default_map2 , alpha , beta , gamma , imsize, ignore_edge_pixels , npol , q , J1);	// check if the step was near optimal
		if(err!=0)
		{
			cout<<endl<<"Error detected in check_step, err = "<<err<<endl<<endl;
			break;
		}

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
				break;
			}

			err = interpolate_residuals( current_residuals , new_residuals, chi2_rms , rescale_step , imsize2 , npol);	// interpolate residuals too
			if(err!=0)
			{
				cout<<endl<<"Error detected in interpolation of new and old residuals, err = "<<err<<endl<<endl;
				break;
			}
		}
		else
		{
			err = copy_model( current_model , new_model );	// just copy new model into old model
			if(err!=0)
			{
				cout<<endl<<"Error detected replacing new and old models, err = "<<err<<endl<<endl;
				break;
			}

			err = copy_model( current_residuals , new_residuals );	// just copy new residuals into old residuals
			if(err!=0)
			{
				cout<<endl<<"Error detected replacing new and old residuals, err = "<<err<<endl<<endl;
				break;
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


		// Check if the current iteration is the best (best average polarisation, if available)

		temp = current_rms[0];
		temp2 = min_rms[0];
		for(k = 1; k < npol ; k++)
		{
			temp += pol_upweight_factor*current_rms[k];
			temp2 += pol_upweight_factor*min_rms[k];
		}


		if( temp < temp2 )
		{
			for(k=0; k<npol; k++)
			{
				min_rms[k] = current_rms[k];
			}
			min_rms_i_ctr = ctr;
		}


		// Test for convergence

		converged = true;
		for(i=0;i<npol;i++)
		{
			converged_temp = ( current_rms[i]  < rms_theoretical[i] );

			if(debug or (ctr%100 == 0) )
			{
				if(converged_temp)
				{
					cout<<"Stokes "<<i+1<<" has converged. "<<current_rms[i]<<" < "<<rms_theoretical[i]<<endl;
				}
				else
				{
					cout<<"Convergence test S"<<i+1<<": "<<current_rms[i]<<" needs to be < "<< rms_theoretical[i] <<endl;
				}
			}

			converged = converged && converged_temp;
		}


		// Add flux condition if required

		if(conserve_flux && (!estimate_flux) )
		{
			converged = converged && ( total_flux - zsf < flux_tolerance * zsf);
			if(debug&&converged)
			{
				cout<<"Flux has converged."<<endl;
			}
		}
//		converged = converged && (grad.JJ / grad.II < convergence_tolerance);

		if(alpha!=alpha)
		{
			err=1;
			cout<<endl<<"Error : Infinity detected in alpha."<<endl<<endl;	// a check to make sure no infinities are around
			goto free_mem_exit;
		}
		if(total_flux < 0)
		{
			err=1;
			cout<<endl<<"Error : Negative total Stokes I flux detected in model."<<endl<<endl;
			goto free_mem_exit;
		}

		if( ctr - min_rms_i_ctr > convergence_limit )
		{
			cout<<endl<<"Convergence appears to have stopped."<<endl;
			cout<<"If you have not achieved the desired convergence try changing some parameters."<<endl;
			break;
		}
		if(debug or (ctr%100 == 0) )
		{
			cout<<endl<<endl;
		}
	}
	
	
	
	
	
	
	
	
/*	
#######################################################################################################################################


	// Exited from main loop
	
########################################################################################################################################
*/







	if(converged)
	{
		cout<<endl<<"Successful convergence after "<<ctr<<" iterations."<<endl;
		for(i=0;i<npol;i++)
		{
			cout<<"Stokes "<<i+1<<" has converged. "<<current_rms[i]<<" < "<<rms_theoretical[i]<<endl;
		}
	}
	else
	{
		cout<<endl<<"Failed to converge after "<<ctr<<" iterations."<<endl;
		for(i=0;i<npol;i++)
		{
			cout<<"Convergence test S"<<i+1<<": "<<current_rms[i]<<" needs to be < "<< rms_theoretical[i] <<endl;
		}
	}

	cout<<endl<<endl<<"Final iteration number "<<ctr<<endl;
	cout<<"Total flux, max and min = "<<total_flux<<" , "<<imax<<" , "<<imin<<endl;
	cout<<"Alpha, beta, gamma =  "<<alpha<<" , "<<beta<<" , "<<gamma<<endl;

	if(debug)
	{
		cout<<"First step, second step, step limit = "<<step_length1<<" , "<<rescale_step<<" , "<<step_limit<<endl;	// output some info
		cout<<"GradJ.J, Grad1.1, J0, J1 = "<<grad.JJ<<" , "<<grad.II<<" , "<<J0<<" , "<<J1<<endl;
		cout<<"Delta E, F, G = "<<delta_E<<" , "<<delta_F<<" , "<<delta_G<<endl;
		cout<<"GradE.E, GradF.F, GradG.G = "<<grad.EE<<" , "<<grad.FF<<" , "<<grad.GG<<endl;
		cout<<"GradE.F, GradF.G, GradE.G = "<<grad.EF<<" , "<<grad.FG<<" , "<<grad.EG<<endl;	// output even more info
		cout<<"GradE.H, GradF.H, GradG.H = "<<grad.EH<<" , "<<grad.FH<<" , "<<grad.GH<<endl;
		cout<<"GradE.J, GradF.J, GradG.J = "<<grad.EJ<<" , "<<grad.FJ<<" , "<<grad.GJ<<endl;
	}


	// find peak of dirty beam for Gaussian clean beam
	
	temp = 0.0;
	for(i=0;i<imsize;i++)
	{
		for(j=0;j<imsize;j++)
		{
			if(dirty_beam[i*imsize+j]>temp)
			{
				temp = dirty_beam[i*imsize+j];
				peak[0] = i;
				peak[1] = j;
			}
		}
	}

	gen_gauss(new_model[0], imsize , fitsi[0].cell_ra, restoring_beam[0] , restoring_beam[1] , restoring_beam[2], peak);	// make restoring beam

	if( pad_factor == 1)
	{
		arrange_ft( new_model[0] , imsize );	// arrange beam so that the convolution will work well
	}

	ft_beam(new_model[0] , dirty_beam_ft , imsize , pad_factor  , forward_transform , double_buff , complex_buff);	// get ft of restoring beam


	for(i=0;i<npol;i++)	// calculate convolved final models
	{
		convolve( current_model[i] , dirty_beam_ft , imsize , pad_factor  , convolved_model[i] , forward_transform , backward_transform , double_buff , complex_buff);
	}


	temp = restoring_beam[0] * restoring_beam[1] * M_PI / ( 4.0 * log(2) * fitsi[0].cell_ra * fitsi[0].cell_dec );	// set time =  number of pixels in restoring beam
	// there is a good possiblity that temp = pixels per beam, but if the restoring beam is different to the initial beam, then the units of Jy/Beam in the residuals need to be converted to Jy/restoring beam.

	temp = temp / pixels_per_beam;

	cout<<endl<<endl<<"Ratio of restoring beam to \"natural\" beam = "<<temp<<endl;

	#pragma omp parallel for collapse(2)
	for(i=0;i<npol;i++)	// add in residuals
	{
		for(j=0;j<imsize2;j++)
		{
			current_residuals[i][j] = - current_residuals[i][j] * temp;	// residuals defined as dirty model - dirty map, redefine to map - model and rescale
			new_model[i][j] = convolved_model[i][j] + current_residuals[i][j];			// save final map in new_model
		}
	}


	cout<<"Writing out maps..."<<endl;


	for(i=0;i<npol;i++)	// write out convolved models
	{
		line.assign(output_name);
		line.append("_Model_S");
		line.append(int2str(i+1));
		line.append(".fits");
		err = quickfits_write_map( line.c_str() , current_model[i] , fitsi[i] , history);
		if(err!=0)
		{
			cout<<endl<<"Error detected in attempting to write to "<<line<<", err = "<<err<<endl<<endl;
		}


		line.assign(output_name);
		line.append("_Convolved_S");
		line.append(int2str(i+1));
		line.append(".fits");
		err = quickfits_write_map( line.c_str() , convolved_model[i] , fitsi[i] , history);
		if(err!=0)
		{
			cout<<endl<<"Error detected in attempting to write to "<<line<<", err = "<<err<<endl<<endl;
		}


		line.assign(output_name);
		line.append("_Residual_S");
		line.append(int2str(i+1));
		line.append(".fits");
		err = quickfits_write_map( line.c_str() , current_residuals[i] , fitsi[i] , history);
		if(err!=0)
		{
			cout<<endl<<"Error detected in attempting to write to "<<line<<", err = "<<err<<endl<<endl;
		}


		line.assign(output_name);
		line.append("_Final_S");
		line.append(int2str(i+1));
		line.append(".fits");
		err = quickfits_write_map( line.c_str() , new_model[i] , fitsi[i] , history);
		if(err!=0)
		{
			cout<<endl<<"Error detected in attempting to write to "<<line<<", err = "<<err<<endl<<endl;
		}

	}







	// free up memory

	free_mem_exit:
	cout<<"Freeing up memory..."<<endl;

	fftw_destroy_plan(forward_transform);
	fftw_destroy_plan(backward_transform);

	fftw_cleanup_threads();

	for(i=0;i<npol;i++)
	{
		delete[] dirty_map[i];
		delete[] current_model[i];
		delete[] new_model[i];
		delete[] current_residuals[i];
		delete[] new_residuals[i];
		delete[] convolved_model[i];
	}
	delete[] dirty_map;
	delete[] current_model;
	delete[] new_model;
	delete[] current_residuals;
	delete[] new_residuals;
	delete[] convolved_model;
	delete[] chi2_rms;

	delete[] default_map2;
	delete[] dirty_beam;
	delete[] mask;
	fftw_free(dirty_beam_ft);
	fftw_free(double_buff);
	fftw_free(complex_buff);


	cout<<"End of program"<<endl;

	if(err == 0 && !converged)
	{
		return(1);
	}
	else
	{
		if(err==0)
		{
			return(0);
		}
		else
		{
			return(2);
		}
	}
}

/*
	write_csv

	write_csv writes out a matrix into a comma separated volume file

	inputs:
		filename = the name of the output file
		array = the matrix to be written out
		imsize = the dimension of the matrix
	outputs:
		on return = err (0 if no error)
*/


int write_csv(string filename, double* array, int imsize)
{
	int i,j,k;
	int err;
	ofstream fout;
	char* cfilename;

	cfilename=new char[filename.length()];		// get name into C format
	strcpy(cfilename,filename.c_str());

	fout.open(cfilename,ios::out);
	err=fout.is_open();

	k=0;
	for(i=0;i<imsize;i++)
	{
		fout<<array[k];
		k++;
		for(j=1;j<imsize;j++)
		{
			fout<<","<<array[k];
			k++;
		}
		fout<<endl;
	}
	fout.close();

	delete[] cfilename;
	return(err);
}

/*
	int2str converts an integer to a string
*/

string int2str (int num)
{
	stringstream ss;

	ss<<num;

	return(ss.str());
}


/*
	"zero_array"

	sets array to zero. useful for initial declaration of residuals

	inputs:

	array = map
	imsize = length of one side of the array

	outputs:

	array = map, now zero everywhere
*/

int zero_array(double* array, int imsize)
{
	#pragma omp parallel for
	for(int i = 0; i < imsize ; i++ )
	{
		array[0] = 0.0;
	}

	return(0);
}

/*
	clip_edges

	clip_edges clips the edges of a map, replacing the values with a single given value

	inputs:
		map = the map in question
		imsize = the dimension of the map (one side)
		replacement_value = the value with which to replace all pixels in the edge region
		edge_limit = the number of pixels from the edge of the map to consider "edge" pixels

	outputs:
		map = the map, now clipped
		return value = 0
*/

int clip_edges(double* map, double replacement_value, int edge_limit, int imsize)
{
	int i,j,k;

	k = imsize - edge_limit;

	#pragma omp parallel for private(j)
	for(i=0;i<imsize;i++)
	{
		for(j=0;j<edge_limit;j++)
		{
			map[ i * imsize + j] = replacement_value;
		}

		for(j = k; j<imsize;j++)
		{
			map[ i * imsize + j] = replacement_value;
		}
	}


	#pragma omp parallel for private(i)
	for(j=edge_limit;j<k;j++)
	{
		for(i=0;i<edge_limit;i++)
		{
			map[ i * imsize + j] = replacement_value;
		}

		for(i = k; i<imsize;i++)
		{
			map[ i * imsize + j] = replacement_value;
		}
	}

	return(0);
}

