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


#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <cmath>
#include <stdlib.h>
#include <sstream>
#include <fftw3.h>
#include <fitsio.h>
#include <stdio.h>
#include <omp.h>

#ifndef grad_def
	#include "gradient_structure.hpp"
#endif
#include "read_driver.hpp"
#include "update_Lagrangian.hpp"
#include "CE_functions.hpp"
#include "FT_convolution.hpp"
#include "map_stats.hpp"

extern "C" {
	#include "quickfits.h"
};