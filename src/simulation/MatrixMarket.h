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

template<typename F> constexpr bool MatrixMarket_IsComplex();
template<> constexpr bool MatrixMarket_IsComplex<float>() { return false; }
template<> constexpr bool MatrixMarket_IsComplex<double>() { return false; }
template<> constexpr bool MatrixMarket_IsComplex<std::complex<float>>() { return true; }
template<> constexpr bool MatrixMarket_IsComplex<std::complex<double>>() { return true; }

template<typename F> constexpr int MatrixMarket_GetPrecision();
template<> constexpr int MatrixMarket_GetPrecision<float>() { return 9; }
template<> constexpr int MatrixMarket_GetPrecision<double>() { return 17; }
template<> constexpr int MatrixMarket_GetPrecision<std::complex<float>>() { return 9; }
template<> constexpr int MatrixMarket_GetPrecision<std::complex<double>>() { return 17; }

template<typename F> void MatrixMarket_WriteValue(std::ostream &stream, F value);
template<> void MatrixMarket_WriteValue<float>(std::ostream &stream, float value) { stream << value; }
template<> void MatrixMarket_WriteValue<double>(std::ostream &stream, double value) { stream << value; }
template<> void MatrixMarket_WriteValue<std::complex<float>>(std::ostream &stream, std::complex<float> value) { stream << value.real() << ' ' << value.imag(); }
template<> void MatrixMarket_WriteValue<std::complex<double>>(std::ostream &stream, std::complex<double> value) { stream << value.real() << ' ' << value.imag(); }

template<typename F, bool ROWMAJOR, bool SYMMETRIC, bool UPPER>
void SaveMatrixMarket(const SparseMatrixBase<F, ROWMAJOR, SYMMETRIC, UPPER> &matrix, const std::string &filename) {
	constexpr bool COMPLEX = MatrixMarket_IsComplex<F>();

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
	stream.precision(MatrixMarket_GetPrecision<F>());

	// write header
	stream << "%%MatrixMarket matrix coordinate";
	stream << ' ' << ((COMPLEX)? "complex" : "real");
	stream << ' ' << ((SYMMETRIC)? "symmetric" : "general");
	stream << '\n';

	// write size
	stream << matrix.GetRows() << ' ' << matrix.GetCols() << ' ' << matrix.GetCoefficients() << '\n';

	// write data (always lower triangular)
	for(size_t col = 0; col < matrix.GetCols(); ++col) {
		for(size_t i = offsets[col]; i < offsets[col + 1]; ++i) {
			stream << indices[i] << ' ' << col << ' ';
			MatrixMarket_WriteValue<F>(stream, values[i]);
			stream << '\n';
		}
	}

}
