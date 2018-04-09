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

#include <fstream>

template<typename F> void MatrixMarket_WriteType(std::ofstream &stream);
template<> void MatrixMarket_WriteType<float>(std::ofstream &stream) { stream << " real"; }
template<> void MatrixMarket_WriteType<double>(std::ofstream &stream) { stream << " real"; }

template<typename F, bool ROWMAJOR, bool SYMMETRIC, bool UPPER>
void SaveMatrixMarket(const SparseMatrixBase<F, ROWMAJOR, SYMMETRIC, UPPER> &matrix, const std::string &filename) {

	// convert to CSC
	std::unique_ptr<size_t[]> offsets(new size_t[matrix.GetCols() + 1]);
	std::unique_ptr<size_t[]> indices(new size_t[matrix.GetCoefficients()]);
	std::unique_ptr<F[]> values(new F[matrix.GetCoefficients()]);
	if(SYMMETRIC && UPPER) {
		matrix.ToCSR(offsets.get(), indices.get(), values.get()); // transpose
	} else {
		matrix.ToCSC(offsets.get(), indices.get(), values.get());
	}

	// open file
	std::ofstream stream;
	stream.open(filename, std::ios_base::out | std::ios_base::trunc);
	if(stream.fail())
		throw std::runtime_error("Could not open file '" + filename + "' for writing.");

	// write header
	stream << "%%MatrixMarket matrix coordinate";
	MatrixMarket_WriteType<F>(stream);
	stream << ((SYMMETRIC)? " symmetric" : " general");
	stream << std::endl;

	// write size
	stream << matrix.GetRows() << ' ' << matrix.GetCols() << ' ' << matrix.GetCoefficients() << std::endl;

	// write data (always lower triangular)
	for(size_t col = 0; col < matrix.GetCols(); ++col) {
		for(size_t i = offsets[col]; i < offsets[col + 1]; ++i) {
			stream << indices[i] << ' ' << col << ' ' << values[i] << std::endl;
		}
	}

}
