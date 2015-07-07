#include "FT_convolution.hpp"
	using namespace std;

/*
	arrange_ft

	arrange_ft changes an image to "wrap-around" order by swapping quadrants 1 --> 3, and 2 --> 4
	This can be useful in performing an FFT
	At the moment this function is needed only for FFTs with no zero padding (not the default case)

	inputs:

	arr = the image (stored in doubles) that needs to be rearranged
	imsize = the dimension of the image (length of one side)

	outputs:

	arr = the image now in "wrap-around" order
*/


void arrange_ft(double* arr, int imsize)
{
	int i,j,k,l;
	int half=(imsize/2);	// 0.5*dimension

	double temp;


	#pragma omp parallel for collapse(2) private(k,l,temp)
	for(i=0;i<half;i++)	// swap values from q1 to q3, and q2 to q4
	{
		for(j=0;j<half;j++)
		{
			k=i*imsize+j;	// location in q1
			l=(i+half)*imsize+half+j;	// location in q3

			temp = arr[k];	// re
			arr[k] = arr[l];
			arr[l] = temp;


			k=i*imsize+j+half;	// location in q2
			l=(i+half)*imsize+j;	// location in q4

			temp = arr[k];	// re
			arr[k] = arr[l];
			arr[l] = temp;


		}
	}
}


/*
	"convolve"

	Convolves data with the repsonse function response and copies output to output.

	Inputs:

	data = data to be convolved
	response = response function (same size as data)
	imsize = size of data, response and output
	pad_factor = amount of zero padding (1 = no zero padding, 2 = imsize zero padding, etc.)
	forward_transform = fftw plan for forward transform
	backward_transform = fftw plan for backward transform
	double_buff = work space for fftw plans
	complex_buff = work space for fftw plans

	Outputs:

	output = the convolved data
	on return, 0
	
*/

int convolve(double* data, fftw_complex* response, int imsize, int pad_factor, double* output , fftw_plan& forward_transform, fftw_plan& backward_transform, double* double_buff, fftw_complex* complex_buff)
{
	int i, j, k;

	int imsize_pad = pad_factor * imsize;
	double temp1 , temp2;



	#pragma omp parallel for collapse(2) private(k)	// initialise (padded) working array for forward transform
	for(i=0;i<imsize_pad;i++)
	{
		for(j=0;j<imsize_pad;j++)
		{
			k =  (i * imsize_pad) + j ;
			if( ( i < imsize ) && ( j <imsize ) )
			{
				double_buff[k] = data[ (i * imsize) + j ];
			}
			else
			{
					double_buff[k] = 0.0;
			}
		}
	}



	fftw_execute( forward_transform );	// forward transform



	j = pad_factor * imsize * ( pad_factor * imsize / 2 + 1 );	// this is the point where Hermitian conjugacy means the FFTW routine has stopped spitting out more data
									// no need to process the additional data, as the c2r transform assumes conjugacy too

	k = imsize_pad * imsize_pad;

	#pragma omp parallel for private(temp1,temp2)
	for( i=0 ; i < j ; i++ )	// convolution theorem
	{
		temp1 = ( complex_buff[i][0] * response[i][0] - complex_buff[i][1] * response[i][1] ) / double( k );
		temp2 = ( complex_buff[i][1] * response[i][0] + complex_buff[i][0] * response[i][1] ) / double( k );	// the k scales data_pad and response_pad correctly

		complex_buff[i][0] = temp1;
		complex_buff[i][1] = temp2;
	}


	fftw_execute( backward_transform );	// inverse transform


	if(pad_factor==1)	// k is the displacement necessary to read from the centre of the padded array
	{
		k = 0;
	}
	else
	{
		k = (pad_factor - 1) * imsize / 2;
	}


	#pragma omp parallel for collapse(2)	// copy data to result, reading from the "centre" of the padded output
	for(i=0;i<imsize;i++)
	{
		for(j=0;j<imsize;j++)
		{
			output[ (i * imsize) + j ] = double_buff[ ( (i + k) * imsize_pad) + j + k];
		}
	}



	return(0);
}

/*
	ft_beam

	ft_beam FFTs a beam (dirty or restoring) and saves the result in a complex array
	As the FFT is real to complex, the array does not need the same dimensions as the input due to the Hermitian conjugacy of the resulting FT frequencies
	This FFT is done outside the convolution function because the beam only needs to be FFTed once. This saves processing time.

	inputs:

	beam = the beam image to be FFTed (dirty or restoring)
	imsize = the length/width of the image
	pad_factor = amount of zero padding (1 = no zero padding, 2 = imsize zero padding, etc.)
	plan = fftw plan for forward transform
	double_buff = work space for fftw plans
	complex_buff = work space for fftw plans

	outputs:

	ft_beam = the FFT of the beam
	on return, 0

*/

int ft_beam(double* beam, fftw_complex* ft_beam, int imsize, int pad_factor, fftw_plan& plan, double* double_buff, fftw_complex* complex_buff)
{
	int imsize_pad = pad_factor * imsize;
	int i, j, k;

	#pragma omp parallel for collapse(2) private(k)	// initialise (padded) working array for transform
	for(i=0;i<imsize_pad;i++)
	{
		for(j=0;j<imsize_pad;j++)
		{
			k =  (i * imsize_pad) + j ;
			if( ( i < imsize ) && ( j < imsize ) )
			{
				double_buff[k] = beam[ (i * imsize) + j ];
			}
			else
			{
				double_buff[k] = 0.0;
			}
		}
	}

//	pad_image( beam , double_buff , imsize, pad_factor);	// initialise (padded) working array for forward transform, using padding correctly

	fftw_execute(plan);	// forward transform

	j = pad_factor * imsize * ( pad_factor * imsize / 2 + 1 );	// only read in as much as necessary, taking the Herm. cong. into account

	for(i=0;i<j;i++)	// copy result to output array
	{
		ft_beam[i][0] = complex_buff[i][0];
		ft_beam[i][1] = complex_buff[i][1];
	}

	return(0);
}

/*
	gen_gauss

	gen_gauss makes a Gaussian distribution using the given beam parameters

	inputs:
		imsize = the dimensions of the matrix
		cellsize = the size of a single cell in degrees
		bmaj = the major axis of the ellipse at the FWHM of the beam in degrees
		bmin = the major axis of the ellipse at the FWHM of the beam in degrees
		bpa = the position angle of the ellipse at the FWHM of the beam in degrees

	output
		matrix = a matrix containing the generated Gaussian
		on return, 0
*/

int gen_gauss(double* matrix, int imsize, double cellsize, double bmaj, double bmin, double bpa, int* peak)
{
	int imsize2=imsize*imsize;
	int i,j,k;
	double x,y,a,b,c;
	int xorigin = peak[0];	// centre of distribution (pixels). Assumed centre is of the form (255,256) for a 512 map.
	int yorigin = peak[1];

	bpa = bpa + 90.0; // convert from astronomy measures (in astro bmaj is on the yaxis)
	
	x = bmaj;
	bmaj = bmin;
	bmin = x;

	bpa*=(M_PI/180.0);	// convert to radiens
	bmaj*=((M_PI/180.0)/(2.354820045));	// convert to radiens from degrees and from FWHM to sigma (2*sqrt(2*log(2)))) = 2.354820045
	bmin*=((M_PI/180.0)/(2.354820045));

	cellsize*=M_PI/180.0;	// convert from degrees to radiens

	a=0.5*(pow(cos(bpa)/bmaj,2)+pow(sin(bpa)/bmin,2));
	b=0.25*sin(2.0*bpa)*(-1.0/pow(bmaj,2)+1.0/pow(bmin,2));
	c=0.5*(pow(sin(bpa)/bmaj,2)+pow(cos(bpa)/bmin,2));


	i=0;
	j=0;

	for(k=0;k<imsize2;k++)
	{
		x=double(i-xorigin)*cellsize;
		y=double(j-yorigin)*cellsize;
		matrix[k]=exp(-(a*pow(x,2)+2.0*b*x*y+c*pow(y,2)));	// 2d gaussian general form
		j++;
		if(j==imsize)
		{
			j=0;
			i++;
		}
	}

	return(0);
}