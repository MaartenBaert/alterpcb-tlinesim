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

void FemMatrix_StaticEPot_Rect(SparseBlockMatrix<complex_t> &matrix, size_t vars[8], real_t delta_x, real_t delta_y,
		complex_t permittivity_x, complex_t permittivity_y);

void FemMatrix_StaticMPot_Rect(SparseBlockMatrix<complex_t> &matrix, size_t vars[8], real_t delta_x, real_t delta_y,
		complex_t permeability_x, complex_t permeability_y);
void FemMatrix_StaticMPot_XLine(SparseBlockMatrix<complex_t> &matrix, size_t vars[4], real_t delta_x, real_t omega, complex_t impedance);
void FemMatrix_StaticMPot_YLine(SparseBlockMatrix<complex_t> &matrix, size_t vars[4], real_t delta_y, real_t omega, complex_t impedance);

void FemMatrix_FullEMPot_Rect(SparseMatrix<complex_t> matrix[2], size_t vars[24], real_t delta_x, real_t delta_y, real_t omega,
		complex_t permittivity_x, complex_t permittivity_y, complex_t permittivity_z,
		complex_t permeability_x, complex_t permeability_y, complex_t permeability_z);
void FemMatrix_FullEMPot_XLine(SparseMatrix<complex_t> matrix[2], size_t vars[8], real_t delta_x, real_t omega, complex_t impedance);
void FemMatrix_FullEMPot_YLine(SparseMatrix<complex_t> matrix[2], size_t vars[8], real_t delta_y, real_t omega, complex_t impedance);
