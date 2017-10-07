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

#include <vector>

// A sparse matrix implementation optimized for the construction of FEM problems.
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
// in the B, C or D blocks, separate from the main A block. The FEM problem is then solved as follows:
//   x = inv(A) * (rhs - B * fixed_values)
// Finally the resulting residual at the fixed nodes can be calculated as follows:
//   residual = C * x + D * fixed_values
// Usually the number of fixed variables is very small, so the B, C and D blocks don't need any special optimizations.

template<typename T>
struct SparseMatrixBulk {
	size_t row;
	T value;
	inline SparseMatrixBulk() {}
	inline SparseMatrixBulk(size_t row, T value) : row(row), value(value) {}
};

template<typename T>
struct SparseMatrixEntry {
	size_t row, col;
	T value;
	inline SparseMatrixEntry(size_t row, size_t col, T value) : row(row), col(col), value(value) {}
};

template<typename T>
struct SparseMatrixHasher {
	inline static bool Equal(const SparseMatrixEntry<T> &a, const SparseMatrixEntry<T> &b) {
		return (a.row == b.row && a.col == b.col);
	}
	inline static hash_t Hash(hash_t hash, const SparseMatrixEntry<T> &value) {
		hash = MurmurHash::HashData(hash, value.row);
		hash = MurmurHash::HashData(hash, value.col);
		return hash;
	}
};

// A symmetric sparse matrix stored in upper triangular format, with diagonal entry values halved.
// The lower triangular part is assumed to be identical to the upper triangular part.
// The 'size' parameter determines only the size of the A block.
template<typename T>
class SymmetricSparseMatrix {

private:
	typedef HashTable<SparseMatrixEntry<T>, SparseMatrixHasher<T>> EntryTable;

private:
	size_t m_size1, m_size2, m_bulk_rows;
	std::vector<SparseMatrixBulk<T>> m_data_a_bulk;
	EntryTable m_data_a_extra;
	EntryTable m_data_bcd;
	size_t m_coefficients_bulk;

public:
	inline SymmetricSparseMatrix() : m_data_a_extra(1024), m_data_bcd(1024) {
		m_size1 = 0;
		m_size2 = 0;
		m_bulk_rows = 0;
		m_coefficients_bulk = 0;
	}

	inline void Reset(size_t size1, size_t size2, size_t bulk_rows) {
		assert(size1 != INDEX_NONE);
		assert(size2 != INDEX_NONE);
		m_size1 = size1;
		m_size2 = size2;
		m_bulk_rows = bulk_rows;
		size_t total_bulk = m_bulk_rows * m_size1;
		m_data_a_bulk.resize(total_bulk);
		for(size_t i = 0; i < total_bulk; ++i) {
			m_data_a_bulk[i].row = INDEX_NONE;
		}
		m_data_a_extra.Clear();
		m_data_bcd.Clear();
		m_coefficients_bulk = 0;
	}

	inline void Add(size_t row, size_t col, T value) {
		assert(row < m_size1 || (row >= INDEX_OFFSET && row < INDEX_OFFSET + m_size2));
		assert(col < m_size1 || (col >= INDEX_OFFSET && col < INDEX_OFFSET + m_size2));
		if(row > col) {
			std::swap(row, col); // make upper triangular
		}
		if(col < INDEX_OFFSET) {
			SparseMatrixBulk<T> *col_bulk = m_data_a_bulk.data() + m_bulk_rows * col;
			for(size_t i = 0; i < m_bulk_rows; ++i) {
				if(col_bulk[i].row == row) {
					col_bulk[i].value += value;
					return;
				}
				if(col_bulk[i].row == INDEX_NONE) {
					col_bulk[i].row = row;
					col_bulk[i].value = value;
					++m_coefficients_bulk;
					return;
				}
			}
		}
		EntryTable &table = (col < INDEX_OFFSET)? m_data_a_extra : m_data_bcd;
		SparseMatrixEntry<T> temp = {row, col, value};
		std::pair<size_t, bool> p = table.TryPushBack(temp);
		if(!p.second) {
			table[p.first].value += value;
		}
	}

	inline size_t GetCoefficientsA() {
		return m_coefficients_bulk + m_data_a_extra.GetSize();
	}

	inline size_t GetCoefficientsBCD() {
		return m_data_bcd.GetSize();
	}

	// Converts the A block to compressed sparse column format.
	template<typename I>
	void MakeCscA(I *offsets, I *indices, T *values) {

		// count coefficients per column
		std::fill_n(offsets, m_size1, 0);
		for(size_t col = 0; col < m_size1; ++col) {
			SparseMatrixBulk<T> *col_bulk = m_data_a_bulk.data() + m_bulk_rows * col;
			for(size_t i = 0; i < m_bulk_rows; ++i) {
				if(col_bulk[i].row == INDEX_NONE)
					break;
				++offsets[col];
			}
		}
		for(size_t i = 0; i < m_data_a_extra.GetSize(); ++i) {
			SparseMatrixEntry<T> &entry = m_data_a_extra[i];
			++offsets[entry.col];
		}

		// calculate offsets
		I total = 0;
		for(size_t i = 0; i < m_size1; ++i) {
			total += offsets[i];
			offsets[i] = total;
		}
		offsets[m_size1] = total;

		// calculate indices and values
		for(size_t col = 0; col < m_size1; ++col) {
			SparseMatrixBulk<T> *col_bulk = m_data_a_bulk.data() + m_bulk_rows * col;
			for(size_t i = 0; i < m_bulk_rows; ++i) {
				SparseMatrixBulk<T> &entry = col_bulk[i];
				if(entry.row == INDEX_NONE)
					break;
				size_t j = --offsets[col];
				indices[j] = (I) entry.row;
				values[j] = (entry.row == col)? 2.0 * entry.value : entry.value;
			}
		}
		for(size_t i = 0; i < m_data_a_extra.GetSize(); ++i) {
			SparseMatrixEntry<T> &entry = m_data_a_extra[i];
			size_t j = --offsets[entry.col];
			indices[j] = (I) entry.row;
			values[j] = (entry.row == entry.col)? 2.0 * entry.value : entry.value;
		}

	}

	// Calculates B * fixed_values and subtracts the result from the right-hand-side matrix.
	void SubtractFromRhs(T *rhs, const T *fixed_values) {
		for(size_t i = 0; i < m_data_bcd.GetSize(); ++i) {
			SparseMatrixEntry<T> &entry = m_data_bcd[i];
			if(entry.row < INDEX_OFFSET) {
				rhs[entry.row] -= entry.value * fixed_values[entry.col - INDEX_OFFSET];
			}
		}
	}

	// Calculates C * x + D * fixed_values.
	void CalculateResidual(T *residual, const T *solution, const T *fixed_values) {
		std::fill_n(residual, m_size2, T());
		for(size_t i = 0; i < m_data_bcd.GetSize(); ++i) {
			SparseMatrixEntry<T> &entry = m_data_bcd[i];
			if(entry.row < INDEX_OFFSET) {
				residual[entry.col - INDEX_OFFSET] += entry.value * solution[entry.row];
			} else {
				residual[entry.col - INDEX_OFFSET] += entry.value * fixed_values[entry.row - INDEX_OFFSET];
				residual[entry.row - INDEX_OFFSET] += entry.value * fixed_values[entry.col - INDEX_OFFSET];
			}
		}
	}

public:
	inline size_t GetSize1() { return m_size1; }
	inline size_t GetSize2() { return m_size2; }

};
