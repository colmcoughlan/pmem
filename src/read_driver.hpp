/*
    pmem_heads.hpp is a header file used by pmem and some programs it links with
    Copyright (C) 2012  Colm Coughlan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include <fstream>
#include <string.h>
#include <iostream>
#include <stdlib.h>

	using namespace std;

int read_driver(string driver_filename, int &npol, string* filename_dirty_map, string &filename_dirty_beam, string &filename_default_map, string &filename_mask, double &zsf, bool &conserve_flux_mode
, double* rms_theoretical, int &niter, double* beam, int* box, double &acc_factor, double &q_factor, double &pol_factor, string &output_name, int &ignore_edge_pixels, bool &debug);