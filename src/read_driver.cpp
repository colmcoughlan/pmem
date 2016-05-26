// This function reads the input from the mem driver file and passes it back to the main MEM program.

#include "read_driver.hpp"
	using namespace std;

int read_driver(string driver_filename, int &npol, string* filename_dirty_map, string &filename_dirty_beam, string &filename_default_map, string &filename_mask, double &zsf, bool &conserve_flux_mode
, double* rms_theoretical, int &niter, double* beam, int* box, double &acc_factor, double &q_factor, double &pol_factor, string &output_name, int &ignore_edge_pixels, bool &nr_solve, bool &debug)
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
	
	getline(fout,filename_mask);	// filename of default map (if any)
	
	getline(fout,line);
	zsf = atof(line.c_str());	// read in estimated zero spacing flux

	getline(fout,line);
	conserve_flux_mode = atoi(line.c_str());

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
	ignore_edge_pixels = atoi(line.c_str());
	
	getline(fout,line);
	if(atoi(line.c_str())==0)	// read in solution mode option (0 = newton-raphson, 1 = steepest-descent)
	{
		nr_solve=true;
	}
	else
	{
		nr_solve=false;
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