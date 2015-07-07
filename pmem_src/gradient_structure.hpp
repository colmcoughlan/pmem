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


struct gradient_structure {
	double EE;
	double EF;
	double EG;
	double EH;
	double EJ;
	double FF;
	double FG;
	double FH;
	double FJ;
	double GG;
	double GH;
	double GJ;
	double HH;
	double JJ;
	double II;
};

#define grad_def