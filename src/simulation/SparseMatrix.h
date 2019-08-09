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
#include "HashTable.h"
#include "Eigen.h"
#include "EigenSparse.h"

#include <type_traits>

// This file contains a sparse block matrix implementation optimized for the construction of FEM problems.
// The matrix is partitioned as follows:
//   M = [ A B ]
//       [ C D ]
// The blocks each have their own range of row and column indices:
//   A: row = 0 .. INDEX_OFFSET - 1, col = 0 .. INDEX_OFFSET - 1
//   B: row = 0 .. INDEX_OFFSET - 1, col = INDEX_OFFSET .. INDEX_NONE - 1
//   C: row = INDEX_OFFSET .. INDEX_NONE - 1, col = 0 .. INDEX_OFFSET - 1
//   D: row = INDEX_OFFSET .. INDEX_NONE - 1, col = INDEX_OFFSET .. INDEX_NONE - 1
// This partitioning scheme simplifies the construction of FEM problems. Usually most of the variables are unknown,
// but a few variables may have a predefined fixed value as a result of boundary conditions or excitations.
// These variables get an index starting at INDEX_OFFSET instead of 0, so all related coefficients will end up
// in the B, C or D blocks, separate from the main A block.
// The FEM problem is then solved as follows:
//   solution = A \ (rhs - B * fixed_values)
// Finally the resulting residual at the fixed nodes can be calculated as follows:
//   residual = C * solution + D * fixed_values
// Alternatively it is possible to solve for all excitation modes simultaneously as follows:
//   solution = A \ -B
//   residual = C * solution + D

// Symmetric sparse matrices are stored in lower or upper triangular format, with diagonal entry values halved.
// The lower and upper triangular parts are assumed to be identical (i.e. not complex conjugates).
// The 'size1' and 'size2' parameters determine the size of the ABCD partitions as explained above.

template<typename F>
class SparseMatrix {

private:
	struct TripletEntry {
		size_t row, col;
		F value;
		inline TripletEntry(size_t row, size_t col, F value) : row(row), col(col), value(value) {}
	};
	struct CompressedEntry {
		size_t index;
		F value;
	};

private:
	size_t m_rows, m_cols;
	std::vector<TripletEntry> m_triplets;

public:
	SparseMatrix() {
		m_rows = 0;
		m_cols = 0;
	}

	void Free() {
		m_cols = 0;
		m_rows = 0;
		m_triplets.clear();
		m_triplets.shrink_to_fit();
	}

	void Reset(size_t rows, size_t cols) {
		m_rows = rows;
		m_cols = cols;
		m_triplets.clear();
	}

	void Insert(size_t row, size_t col, F value) {
		assert(row < m_rows);
		assert(col < m_cols);
		m_triplets.emplace_back(row, col, value);
	}

	// Convert to compressed sparse row/column format.
	// For symmetric matrices, only the lower or upper part will be stored, but the diagonal values will be doubled.
	template<class EigenSparseMatrix>
	void ToEigen(EigenSparseMatrix &output) const {

		// count inner coefficients
		std::vector<size_t> inner_offsets((EigenSparseMatrix::IsRowMajor)? m_cols : m_rows, 0);
		for(const TripletEntry &entry : m_triplets) {
			++inner_offsets[(EigenSparseMatrix::IsRowMajor)? entry.col : entry.row];
		}
		size_t inner_total = 0;
		for(size_t inner = 0; inner < inner_offsets.size(); ++inner) {
			size_t count = inner_offsets[inner];
			inner_offsets[inner] = inner_total;
			inner_total += count;
		}

		// group by inner index
		std::vector<size_t> inner_indices(m_triplets.size());
		std::vector<F> inner_values(m_triplets.size());
		for(const TripletEntry &entry : m_triplets) {
			size_t k = inner_offsets[(EigenSparseMatrix::IsRowMajor)? entry.col : entry.row]++;
			inner_indices[k] = (EigenSparseMatrix::IsRowMajor)? entry.row : entry.col;
			inner_values[k] = entry.value;
		}

		// deduplicate and count outer coefficients
		std::vector<size_t> outer_offsets(((EigenSparseMatrix::IsRowMajor)? m_rows : m_cols) + 1, 0);
		std::vector<size_t> outer_last((EigenSparseMatrix::IsRowMajor)? m_rows : m_cols, INDEX_NONE);
		size_t pos = 0;
		size_t n1 = 0, n2 = 0;
		for(size_t inner = 0; inner < inner_offsets.size(); ++inner) {
			for(size_t end = inner_offsets[inner]; pos < end; ++pos) {
				size_t outer = inner_indices[pos];
				if(outer_last[outer] != inner) {
					assert(outer_last[outer] == INDEX_NONE || outer_last[outer] < inner);
					++outer_offsets[outer];
					outer_last[outer] = inner;
					++n1;
				} else {
					++n2;
				}
			}
		}
		size_t outer_total = 0;
		for(size_t outer = 0; outer < ((EigenSparseMatrix::IsRowMajor)? m_rows : m_cols); ++outer) {
			size_t count = outer_offsets[outer];
			outer_offsets[outer] = outer_total;
			outer_total += count;
		}

		// allocate output buffers
		output.resize((Eigen::Index) GetRows(), (Eigen::Index) GetCols());
		output.resizeNonZeros((Eigen::Index) outer_total);
		typename EigenSparseMatrix::StorageIndex *offsets = output.outerIndexPtr();
		typename EigenSparseMatrix::StorageIndex *indices = output.innerIndexPtr();
		typename EigenSparseMatrix::Scalar *values = output.valuePtr();

		// write offsets
		for(size_t outer = 0; outer < ((EigenSparseMatrix::IsRowMajor)? m_rows : m_cols); ++outer) {
			offsets[outer] = (typename EigenSparseMatrix::StorageIndex) outer_offsets[outer];
		}
		offsets[(EigenSparseMatrix::IsRowMajor)? m_rows : m_cols] = (typename EigenSparseMatrix::StorageIndex) outer_total;

		// deduplicate and group by outer index
		pos = 0;
		std::fill_n(outer_last.data(), outer_last.size(), INDEX_NONE);
		for(size_t inner = 0; inner < inner_offsets.size(); ++inner) {
			for(size_t end = inner_offsets[inner]; pos < end; ++pos) {
				size_t outer = inner_indices[pos];
				if(outer_last[outer] != inner) {
					size_t k = outer_offsets[outer]++;
					outer_last[outer] = inner;
					indices[k] = (typename EigenSparseMatrix::StorageIndex) inner;
					values[k] = (typename EigenSparseMatrix::Scalar) inner_values[pos];
				} else {
					size_t k = outer_offsets[outer] - 1;
					values[k] += (typename EigenSparseMatrix::Scalar) inner_values[pos];
				}
			}
		}

	}

public:
	inline size_t GetRows() const { return m_rows; }
	inline size_t GetCols() const { return m_cols; }

};

template<typename F>
class SparseBlockMatrix {

private:
	SparseMatrix<F> m_matrix_a, m_matrix_b, m_matrix_c, m_matrix_d;

public:
	void Free() {
		m_matrix_a.Free();
		m_matrix_b.Free();
		m_matrix_c.Free();
		m_matrix_d.Free();
	}

	void Reset(size_t rows1, size_t rows2, size_t cols1, size_t cols2) {
		assert(rows1 <= INDEX_OFFSET);
		assert(rows2 < INDEX_OFFSET);
		assert(cols1 <= INDEX_OFFSET);
		assert(cols2 < INDEX_OFFSET);
		m_matrix_a.Reset(rows1, cols1);
		m_matrix_b.Reset(rows1, cols2);
		m_matrix_c.Reset(rows2, cols1);
		m_matrix_d.Reset(rows2, cols2);
	}

	void Insert(size_t row, size_t col, F value) {
		assert(row < m_matrix_a.GetRows() || (row >= INDEX_OFFSET && row < INDEX_OFFSET + m_matrix_d.GetRows()));
		assert(col < m_matrix_a.GetCols() || (col >= INDEX_OFFSET && col < INDEX_OFFSET + m_matrix_d.GetCols()));
		if(row < INDEX_OFFSET) {
			if(col < INDEX_OFFSET) {
				m_matrix_a.Insert(row, col, value);
			} else {
				m_matrix_b.Insert(row, col - INDEX_OFFSET, value);
			}
		} else {
			if(col < INDEX_OFFSET) {
				m_matrix_c.Insert(row - INDEX_OFFSET, col, value);
			} else {
				m_matrix_d.Insert(row - INDEX_OFFSET, col - INDEX_OFFSET, value);
			}
		}
	}

public:
	inline size_t GetRows1() const { return m_matrix_a.GetRows(); }
	inline size_t GetRows2() const { return m_matrix_d.GetRows(); }
	inline size_t GetCols1() const { return m_matrix_a.GetCols(); }
	inline size_t GetCols2() const { return m_matrix_d.GetCols(); }

	inline const SparseMatrix<F>& GetMatrixA() const { return m_matrix_a; }
	inline const SparseMatrix<F>& GetMatrixB() const { return m_matrix_b; }
	inline const SparseMatrix<F>& GetMatrixC() const { return m_matrix_c; }
	inline const SparseMatrix<F>& GetMatrixD() const { return m_matrix_d; }

};
