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

#include <cmath>
#include <iostream>
#ifndef grad_def
	#include "gradient_structure.hpp"
#endif
	using namespace std;

extern "C" {
	void dgetrf_( const int * , const int * , double * , const int * , int * , int * );
	void dgetrs_( const char * , const int * , const int * , double * , const int * , int * , double * , const int* , int * );
};

int new_ABG(gradient_structure grad , double delta_E , double delta_F , double delta_G , double& alpha , double& beta , double& gamma , bool conserve_flux , int npol, bool force_chi2_method);