#include "update_Lagrangian.hpp"

/*
	new_ABG

	new_ABG gets new values for alpha, beta and gamma
	This function uses LAPACK functions dgetrf (factorise into LU) and dgetrs (solve A x = b when A is LU factorised)

	inputs:
		grad = a structure containing all the norms calculated in get_info
		delta_E = the difference between current residuals and expected residuals for Stokes I
		delta_F = the difference between current residuals and expected residuals for Stokes Q and U
		delta_G = the difference between current flux and expected flux
		conserve_flux = a boolean that is true if flux is being conserved, false otherwise
		npol = the number of polarisations being deconvolved
	outputs:
		alpha = the Lagrangian parameter associated with minimising Stokes I residuals
		beta = the Lagrangian parameter associated with minimising Stokes Q and U residuals
		gamma = the Lagrangian parameter associated with conserving the Stokes I flux
		on return, 0
*/


int new_ABG(gradient_structure grad , double delta_E , double delta_F , double delta_G , double& alpha , double& beta , double& gamma , bool conserve_flux , int npol, bool force_chi2_method)
{
	double l , arg;
	double alpha1 , beta1 , gamma1;
	double alpha2 , beta2 , gamma2;
	double dalpha , dbeta , dgamma;

	double quad_tol = 0.01;	// tolerance for checking if the changes in alpha etc. are "small enough". This is the MEM's equivalent of Gain. Keep very well under one.
	double chi2_mode_tol = 0.05;	// MEM gain tolerance. If lower than this MEM is progressing slowly.

	int err;

	int n;	// parameters for LAPACK functions dgetrf (factorise into LU) and dgetrs (solve A x = b when A is LU factorised)
	char trans='T';	// indicate that the transpose should be used in LAPACK functions (row major to column major)
	int nrhs = 1;
	
	beta2 = 0.0;
	gamma2 = 0.0;


	l = fabs( grad.JJ / grad.II );	// a measure of how quickly the cost function is changing scaled with how quickly the I map is changing

	if(alpha <= 0)
	{
		l=0.0;
	}

	if(npol ==1 )
	{
		if(conserve_flux)
		{
			n=2;
			double A[4] = { grad.EE , grad.EG , grad.EG , grad.GG};
			double b[2] = { grad.EH , grad.GH};
			int ipiv[2];


			dgetrf_(&n , &n , A , &n , ipiv, &err );
			if(err != 0)
			{
				cout<<"Error : new_ABG : dgetrf_ : Error performing LU factorisation."<<endl;
				return(1);
			}


			dgetrs_( &trans , &n , &nrhs , A , &n , ipiv , b , &n , &err );
			if(err != 0)
			{
				cout<<"Error : new_ABG : dgetrs_ : Error solving linear equation."<<endl;
				return(1);
			}

			alpha1 = b[0];
			beta1 = 0.0;
			gamma1 = b[1];


			b[0] = delta_E + grad.EJ;	//	reset b
			b[1] = delta_G + grad.GJ;

			dgetrs_( &trans , &n , &nrhs , A , &n , ipiv , b , &n , &err );
			if(err != 0)
			{
				cout<<"Error : new_ABG : dgetrs_ : Error solving linear equation."<<endl;
				return(1);
			}


			dalpha = b[0];
			dbeta = 0.0;
			dgamma = b[1];
		}
		else
		{
			// no beta or gamma
			// need only to solve gradE.gradJ = 0, i.e. gradE.(gradH - alpha * grad E)=0, or alpha = grad HE / grad EE

			alpha1 = grad.EH / grad.EE;
			beta1 = 0.0;
			gamma1 = 0.0;

			// need to solve grad EE * dalpha = delta E + grad EJ

			dalpha = ( delta_E + grad.EJ ) / grad.EE;
			dbeta = 0.0;
			dgamma = 0.0;
		}
	}
	else
	{
		if(conserve_flux)
		{
			// If conserving flux solve 3x3 matrix to find changes in alpha, beta, gamma


			n=3;
			double A[9] = { grad.EE , grad.EF , grad.EG , grad.EF , grad.FF , grad.FG , grad.EG , grad.FG , grad.GG };
			double b[3] = { grad.EH , grad.FH , grad. GH};
			int ipiv[3];


			dgetrf_(&n , &n , A , &n , ipiv, &err );
			if(err != 0)
			{
				cout<<"Error : new_ABG : dgetrf_ : Error performing LU factorisation."<<endl;
				return(1);
			}

			dgetrs_( &trans , &n , &nrhs , &A[0] , &n , ipiv , &b[0] , &n , &err );
			if(err != 0)
			{
				cout<<"Error : new_ABG : dgetrs_ : Error solving linear equation."<<endl;
				return(1);
			}

			alpha1 = b[0];
			beta1 = b[1];
			gamma1 = b[2];


			b[0] = (delta_E + grad.EJ);	//	reset b
			b[1] = (delta_F + grad.FJ);
			b[2] = (delta_G + grad.GJ);

			dgetrs_( &trans , &n , &nrhs , A , &n , ipiv , b , &n , &err );

			if(err != 0)
			{
				cout<<"Error : new_ABG : dgetrs_ : Error solving linear equation."<<endl;
				return(1);
			}

			dalpha = b[0];
			dbeta = b[1];
			dgamma = b[2];
		}
		else
		{
			n=2;
			double A[4] = { grad.EE , grad.EF , grad.EF , grad.FF };
			double b[2] = { grad.EH , grad.FH};
			int ipiv[2];



			dgetrf_(&n , &n , A , &n , ipiv, &err );
			if(err != 0)
			{
				cout<<"Error : new_ABG : dgetrf_ : Error performing LU factorisation."<<endl;
				return(1);
			}
		

			dgetrs_( &trans , &n , &nrhs , A , &n , ipiv , b , &n , &err );
			if(err != 0)
			{
				cout<<"Error : new_ABG : dgetrs_ : Error solving linear equation."<<endl;
				return(1);
			}


			alpha1 = b[0];
			beta1 = b[1];
			gamma1 = 0.0;


			b[0] = delta_E + grad.EJ;	//	reset b
			b[1] = delta_F + grad.FJ;

			dgetrs_( &trans , &n , &nrhs , A , &n , ipiv , b , &n , &err );
			if(err != 0)
			{
				cout<<"Error : new_ABG : dgetrs_ : Error solving linear equation."<<endl;
				return(1);
			}

			dalpha = b[0];
			dbeta = b[1];
			dgamma = 0.0;
		}
	}



	// check if change in alpha is "small enough"

	arg = grad.EJ * grad.EJ - (grad.JJ - quad_tol * grad.II) * grad.EE;
	if( arg > 0 )
	{
		arg = sqrt(arg);
		dalpha = max( (grad.EJ - arg) / grad.EE , min( (grad.EJ + arg) / grad.EE , dalpha ) );
	}
	else
	{
		if(alpha != 0)
		{
			dalpha = 0.0;
		}
		else
		{
			dalpha = 1.0;
			cout<<"Alpha = 0 detected. Accepting unconstrained new value of alpha = "<<dalpha<<endl;
		}
	}
	alpha2 = alpha + dalpha;

	// check if change in beta is "small enough"

	if(npol>1)
	{
		arg = grad.FJ * grad.FJ - (grad.JJ - quad_tol * grad.II) * grad.FF;
		if( arg > 0 )
		{
			arg = sqrt(arg);
			dbeta = max( (grad.FJ - arg) / grad.FF , min( (grad.FJ + arg) / grad.FF , dbeta ) );
		}
		else
		{
			dbeta = 0.0;
		}
		beta2 = beta + dbeta;
	}

	// check if change in gamma is "small enough"

	if(conserve_flux)
	{
		arg = grad.GJ * grad.GJ - (grad.JJ - quad_tol * grad.II) * grad.GG;
		if( arg > 0 )
		{
			arg = sqrt(arg);
			dgamma = max( (grad.GJ - arg) / grad.GG , min( (grad.GJ + arg) / grad.GG , dgamma ) );
		}
		else
		{
			dgamma = 0.0;
		}
		gamma2 = gamma + dgamma;
	}



	// decide on alpha
/*
	cout<<"l = "<<l<<", chi2_mode_tol = "<<chi2_mode_tol<<endl;
	cout<<"alpha1 = "<<alpha1<<", alpha2 = "<<alpha2<<endl;
	cout<<"beta1 = "<<beta1<<", beta2 = "<<beta2<<endl;
	cout<<"gamma1 = "<<gamma1<<", gamma2 = "<<gamma2<<endl;
*/

	if( (l < chi2_mode_tol or alpha1 <=0.0) or force_chi2_method )
	{
		alpha=max(alpha2,0.0);
	}
	else
	{
		alpha=max(alpha1,0.0);	// watch out if falling back to alpha2 when alpha 1 is zero. Alpha 2 can be very large and throw things out.
	}


	// decide on beta


	if(npol>1)
	{
		if( (l < chi2_mode_tol or beta1 <=0.0) or force_chi2_method)
		{
			beta=max(beta2,0.0);
		}
		else
		{
			beta=max(beta1,0.0);
		}
	}

	// decide on gamma

	if(conserve_flux)
	{
		if( (l < chi2_mode_tol or gamma1 <=0.0) or force_chi2_method )
		{
			gamma=max(gamma2, 0.0);
		}
		else
		{
			gamma=max(gamma1, 0.0);
		}
	}

	return(0);
}