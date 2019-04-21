/*
Copyright (C) 2016  The AlterPCB team
Contact: Maarten Baert <maarten-baert@hotmail.com>

This file is part of AlterPCB.

AlterPCB is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

AlterPCB is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this AlterPCB.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include "Basics.h"
#include "SparseMatrix.h"

void FemMatrix_EMPot_Rect(SparseMatrixC<complex_t> matrix[3], size_t vars[12], real_t delta_x, real_t delta_y,
		complex_t permittivity_x, complex_t permittivity_y, complex_t permittivity_z,
		complex_t permeability_ref, real_t omega);

void FemMatrix_EMPot_XLine(SparseMatrixC<complex_t> matrix[3], size_t vars[5], real_t delta_x,
		real_t conductivity_x, real_t conductivity_z, real_t omega);

void FemMatrix_EMPot_YLine(SparseMatrixC<complex_t> matrix[3], size_t vars[5], real_t delta_y,
		real_t conductivity_y, real_t conductivity_z, real_t omega);
