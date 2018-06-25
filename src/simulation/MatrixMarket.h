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
#include "Eigen.h"
#include "EigenSparse.h"

#include <fstream>

namespace MatrixMarket {

template<typename F> constexpr bool IsComplex();
template<> constexpr bool IsComplex<float>() { return false; }
template<> constexpr bool IsComplex<double>() { return false; }
template<> constexpr bool IsComplex<std::complex<float>>() { return true; }
template<> constexpr bool IsComplex<std::complex<double>>() { return true; }

template<typename F> constexpr int GetPrecision();
template<> constexpr int GetPrecision<float>() { return 9; }
template<> constexpr int GetPrecision<double>() { return 17; }
template<> constexpr int GetPrecision<std::complex<float>>() { return 9; }
template<> constexpr int GetPrecision<std::complex<double>>() { return 17; }

template<typename F> void WriteValue(std::ostream &stream, F value);
template<> void WriteValue<float>(std::ostream &stream, float value) { stream << value; }
template<> void WriteValue<double>(std::ostream &stream, double value) { stream << value; }
template<> void WriteValue<std::complex<float>>(std::ostream &stream, std::complex<float> value) { stream << value.real() << ' ' << value.imag(); }
template<> void WriteValue<std::complex<double>>(std::ostream &stream, std::complex<double> value) { stream << value.real() << ' ' << value.imag(); }

template<typename F>
void Save(const std::string &filename, const Eigen::SparseMatrix<F> &matrix, bool symmetric) {

	// open file
	std::ofstream stream;
	stream.open(filename, std::ios_base::out | std::ios_base::trunc);
	if(stream.fail())
		throw std::runtime_error("Could not open file '" + filename + "' for writing.");
	stream.precision(GetPrecision<F>());

	// write header
	stream << "%%MatrixMarket matrix coordinate";
	stream << ' ' << ((IsComplex<F>())? "complex" : "real");
	stream << ' ' << ((symmetric)? "symmetric" : "general");
	stream << '\n';

	// write size
	stream << matrix.rows() << ' ' << matrix.cols() << ' ' << matrix.nonZeros() << '\n';

	// write data (always lower triangular)
	auto offsets = matrix.outerIndexPtr();
	auto indices = matrix.innerIndexPtr();
	auto values = matrix.valuePtr();
	for(size_t outer = 0; outer < (size_t) matrix.outerSize(); ++outer) {
		for(size_t inner = (size_t) offsets[outer]; inner < (size_t) offsets[outer + 1]; ++inner) {
			if(matrix.IsRowMajor) {
				stream << outer << ' ' << indices[inner] << ' ';
			} else {
				stream << indices[inner] << ' ' << outer << ' ';
			}
			WriteValue<F>(stream, values[inner]);
			stream << '\n';
		}
	}

}

}
