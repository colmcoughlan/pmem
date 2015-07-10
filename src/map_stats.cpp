#include "map_stats.hpp"
	using namespace std;

/*
	average_region

	average_region takes the average of a given region in the map

	inputs:
		map = the map in question
		imsize = the size of the map (one side)
		blcx = the x coord of the bottom left corner
		blcy = the y coord of the bottom left corner
		trcx = the x coord of the top right corner
		trcy = the y coord of the top right corner

	outputs:
		on return = the average of the region
*/

double average_region(double* map, int blcx, int blcy, int trcx, int trcy, int imsize)
{
	int i,j;
	double sum = 0.0;

	trcx++;	// increment ends to make the for loop design a bit easier
	trcy++;

	#pragma omp parallel for collapse(2) reduction(+:sum)
	for(i=blcx;i<trcx;i++)
	{
		for(j=blcy;j<trcy;j++)
		{
			sum += map[j * imsize + i];	// add up over entire region, with openmp if possible
		}
	}

	i = trcx - blcx-1;	// find number of pixels in region
	j = trcy - blcy-1;

	i = i *j;

	if(i!=0)
	{
		return(sum/double(i));
	}
	else
	{
		return(map[blcy * imsize + blcx]);
	}
}

/*
	rms_region_region

	rms_region takes the root mean square of a given region in the map

	inputs:
		map = the map in question
		imsize = the size of the map (one side)
		blcx = the x coord of the bottom left corner
		blcy = the y coord of the bottom left corner
		trcx = the x coord of the top right corner
		trcy = the y coord of the top right corner

	outputs:
		on return = the rms of the region
*/

double rms_region(double* map, int blcx, int blcy, int trcx, int trcy, int imsize)
{
	int i,j;
	double sum = 0.0;
	double mean;
	
	mean = average_region( map,  blcx,  blcy,  trcx,  trcy,  imsize );

	trcx++;	// increment ends to make the for loop design a bit easier
	trcy++;

	#pragma omp parallel for collapse(2) reduction(+:sum)
	for(i=blcx;i<trcx;i++)
	{
		for(j=blcy;j<trcy;j++)
		{
			sum += pow(map[j * imsize + i] - mean , 2.0);
		}
	}

	i = trcx - blcx - 1;	// find number of pixels in region
	j = trcy - blcy - 1;

	i = i *j;

	if(i!=0)
	{
		return(sqrt(sum/double(i)));
	}
	else
	{
		return(map[blcy * imsize + blcx] - mean);
	}
}
