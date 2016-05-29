// This function reads the input from the mem driver file and passes it back to the main MEM program.

#include "read_driver.hpp"
	using namespace std;

int read_driver(string driver_filename, int &npol, string* filename_dirty_map, string &filename_dirty_beam, string &filename_default_map, string &filename_mask, double &zsf, bool &conserve_flux, bool& estimate_flux
, double* rms_theoretical, int &niter, double* beam, int* box, double &acc_factor, double &q_factor, double &pol_factor, string &output_name, int &ignore_edge_pixels, int &nr_solve, bool& do_bfgs, bool &debug)
{
	fstream fout;	// file stream
	string line;	// temporary line
	int i;			// for loop iteration

	fout.open(driver_filename.c_str(),ios::in);	// read in filename, number of frequencies to be used and number of 3x3 grids required
	
	if (!fout.is_open())
	{
		cout<<"\t read_driver: Error opening "<<driver_filename<<endl;
		return(1);
	}
	
	getline(fout,line);
	npol = atoi(line.c_str());		// read number of pols

	getline(fout,filename_dirty_map[0]);	// read imap filename
	for(i=1;i<npol;i++)
	{
		getline(fout,filename_dirty_map[i]);	// read in the filenames of the Q,U,V files
	}
	for(i=npol;i<4;i++)
	{
		getline(fout,line);	// read in any blank lines
	}

	getline(fout,filename_dirty_beam);	// filename of dirty beam
	
	getline(fout,filename_default_map);	// filename of default map (if any)
	
	getline(fout,filename_mask);	// filename of mask (if any)
	
	getline(fout,line);
	zsf = atof(line.c_str());	// read in estimated zero spacing flux

	getline(fout,line);
	i = atoi(line.c_str());
	switch(i)
	{
		case 0:
			conserve_flux = false;
			estimate_flux = false;
			break;

		case 1:
			conserve_flux = true;
			estimate_flux = false;
			break;

		case 2:
			conserve_flux = true;
			estimate_flux = true;
			break;

		default:
			cout<<"Error reading flux control parameter."<<endl;
			return(1);
	}

	for(i=0;i<npol;i++)
	{
		getline(fout,line);
		rms_theoretical[i] = atof(line.c_str());	// read in estimated fluxes
	}

	for(i=npol;i<4;i++)
	{
		getline(fout,line);	// read in any blank lines
	}

	getline(fout,line);
	niter = atoi(line.c_str());	// read in number of iterations

	getline(fout,line);
	beam[0] = atof(line.c_str());	// read in restoring beam (bmaj, bmin,bpa)

	getline(fout,line);
	beam[1] = atof(line.c_str());

	getline(fout,line);
	beam[2] = atof(line.c_str());

	getline(fout,line);
	box[0] = atoi(line.c_str());	// read in noise box info (blcx, blcy, trcx, trcy)

	getline(fout,line);
	box[1] = atoi(line.c_str());	// read in blc and trc. Assume to 0 base.

	getline(fout,line);
	box[2] = atoi(line.c_str());

	getline(fout,line);
	box[3] = atoi(line.c_str());

	getline(fout,line);	// read in acceleration factor name
	acc_factor = atof(line.c_str());

	getline(fout,line);	// read in q factor name
	q_factor = atof(line.c_str());

	getline(fout,line);	// read in polarisation upweight factor name
	pol_factor = atof(line.c_str());

	getline(fout,output_name);	// read in output name

	getline(fout,line);
	ignore_edge_pixels = atoi(line.c_str());	// read in number of edge pixels to be clipped
	
	getline(fout,line);
	nr_solve = atoi(line.c_str());	// read in solution mode option (0 = newton-raphson, 1 = BFGS, 2 = DFP, 3 = steepest-descent, 4 = conj_grad)
	if(nr_solve == 1)
	{
		do_bfgs = true;
	}
	else
	{
		do_bfgs = false;
		if(nr_solve == 2)
		{
			nr_solve = 1; 	// internally use the same nr_solve value for all quasi newton methods, with do_bfgs differentiating the two
		}
	}

	getline(fout,line);
	if(atoi(line.c_str())==1)	// read in debug mode option
	{
		debug=true;
	}
	else
	{
		debug=false;
	}


	fout.close();
	return(0);
}