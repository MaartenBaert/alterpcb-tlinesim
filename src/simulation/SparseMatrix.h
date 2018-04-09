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

#include <type_traits>
#include <vector>

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

// Sparse matrices are usually constructed in triplet form, then converted to CSR/CSC by sorting and deduplicating.
// Unfortunately this is quite slow when there are a lot of duplicates, and can consume a lot of memory. In that case,
// a better strategy is to use a hash table to deduplicate the coefficients on the fly. The downside of this approach is
// that when the matrix becomes too large, the hash table buckets no longer fit in the cache and the hash table becomes
// very slow. In order to avoid this problem, the A block uses a 'bulk storage' optimization: since the typical number
// of coefficients per column is usually known for FEM problems, we can preallocate memory for these coefficients as a
// flat array instead of a hash table. This greatly improves the locality, assuming that the matrix is constructed in
// nearly sequential order. For the few columns that might need extra storage, we can fall back to a hash table.
// Usually the number of fixed variables is very small, so the B, C and D blocks don't need any special optimizations.

// Symmetric sparse matrices are stored in lower or upper triangular format, with diagonal entry values halved.
// The lower and upper triangular parts are assumed to be identical (i.e. not complex conjugates).
// The 'size1' and 'size2' parameters determine the size of the ABCD partitions as explained above.

// Insertion sort for data stored in two separate arrays. Useful for sorting coefficients in CSR/CSC format.
template<typename A, typename B>
void InsertionSortPairs(A *data_a, B *data_b, size_t size) {
	for(size_t i = 1; i < size; ++i) {
		if(data_a[i - 1] > data_a[i]) {
			A value_a = data_a[i];
			B value_b = data_b[i];
			data_a[i] = data_a[i - 1];
			data_b[i] = data_b[i - 1];
			size_t j = i - 1;
			while(j != 0 && data_a[j - 1] > value_a) {
				data_a[j] = data_a[j - 1];
				data_b[j] = data_b[j - 1];
				--j;
			}
			data_a[j] = value_a;
			data_b[j] = value_b;
		}
	}
}

template<typename F, bool ROWMAJOR>
class DenseViewBase {

private:
	F *m_data;
	size_t m_stride;

public:
	DenseViewBase(F *data, size_t stride) : m_data(data), m_stride(stride) {}

	void SetZero(size_t rows, size_t cols) {
		std::fill_n(m_data, (ROWMAJOR)? rows * m_stride : cols * m_stride, F());
	}

	void Negate(size_t rows, size_t cols) {
		size_t n = (ROWMAJOR)? rows * m_stride : cols * m_stride;
		for(size_t i = 0; i < n; ++i) {
			m_data[i] = -m_data[i];
		}
	}

	DenseViewBase<F, !ROWMAJOR> TransposedView() {
		return DenseViewBase<F, !ROWMAJOR>(m_data, m_stride);
	}

public:
	inline F* GetData() { return m_data; }
	inline size_t GetStride() { return m_stride; }
	inline size_t GetTotalSize(size_t rows, size_t cols) { return (ROWMAJOR)? rows * m_stride : cols * m_stride; }
	inline F& operator()(size_t row, size_t col) { return m_data[(ROWMAJOR)? col + row * m_stride : row + col * m_stride]; }

};

template<typename F>
inline DenseViewBase<F, false> DenseViewC(F *data, size_t stride) {
	return DenseViewBase<F, false>(data, stride);
}

template<typename F>
inline DenseViewBase<F, true> DenseViewR(F *data, size_t stride) {
	return DenseViewBase<F, true>(data, stride);
}

template<typename F, bool ROWMAJOR, bool SYMMETRIC, bool UPPER>
class SparseMatrixBase {

private:
	struct BulkEntry {
		size_t inner;
		F value;
		inline BulkEntry() : inner(INDEX_NONE) {}
	};

	struct TableEntry {
		size_t outer, inner;
		F value;
		inline TableEntry(size_t outer, size_t inner, F value) : outer(outer), inner(inner), value(value) {}
	};

	struct TableHasher {
		inline static bool Equal(const TableEntry &a, const TableEntry &b) {
			return (a.outer == b.outer && a.inner == b.inner);
		}
		inline static hash_t Hash(hash_t hash, const TableEntry &value) {
			hash = MurmurHash::HashData(hash, value.outer);
			hash = MurmurHash::HashData(hash, value.inner);
			return hash;
		}
	};

private:
	size_t m_outer_size, m_inner_size, m_bulk_size, m_bulk_coefficients;
	std::unique_ptr<BulkEntry[]> m_bulk_data;
	HashTable<TableEntry, TableHasher> m_table_data;

private:
	inline static size_t Row(size_t outer, size_t inner) { return (ROWMAJOR)? outer : inner; }
	inline static size_t Col(size_t outer, size_t inner) { return (ROWMAJOR)? inner : outer; }
	inline static size_t Outer(size_t row, size_t col) { return (ROWMAJOR)? row : col; }
	inline static size_t Inner(size_t row, size_t col) { return (ROWMAJOR)? col : row; }

public:
	SparseMatrixBase() {
		m_outer_size = 0;
		m_inner_size = 0;
		m_bulk_size = 0;
		m_bulk_coefficients = 0;
	}

	void Free() {
		m_outer_size = 0;
		m_inner_size = 0;
		m_bulk_size = 0;
		m_bulk_coefficients = 0;
		m_bulk_data.reset();
		m_table_data.Free();
	}

	void Reset(size_t rows, size_t cols, size_t bulk_size) {
		assert(!SYMMETRIC || rows == cols);
		m_outer_size = Outer(rows, cols);
		m_inner_size = Inner(rows, cols);
		m_bulk_size = bulk_size;
		m_bulk_coefficients = 0;
		m_bulk_data.reset(new BulkEntry[m_bulk_size * m_outer_size]());
		m_table_data.Clear();
	}

	// Insert a new coefficient. If it already exists, the new value will be added to the existing one.
	template<bool MAKE_SYMMETRIC=true>
	void Insert(size_t row, size_t col, F value) {
		size_t outer = Outer(row, col), inner = Inner(row, col);
		assert(outer < m_outer_size);
		assert(inner < m_inner_size);
		if(SYMMETRIC) {
			if(MAKE_SYMMETRIC) {
				if(UPPER) {
					if(row > col)
						std::swap(outer, inner);
				} else {
					if(row < col)
						std::swap(outer, inner);
				}
			} else {
				if(UPPER) {
					assert(row <= col);
				} else {
					assert(row >= col);
				}
			}
		}
		BulkEntry *bulk = m_bulk_data.get() + m_bulk_size * outer;
		for(size_t i = 0; i < m_bulk_size; ++i) {
			BulkEntry &entry = bulk[i];
			if(entry.inner == inner) {
				entry.value += value;
				return;
			}
			if(entry.inner == INDEX_NONE) {
				entry.inner = inner;
				entry.value = value;
				++m_bulk_coefficients;
				return;
			}
		}
		TableEntry temp = {outer, inner, value};
		std::pair<size_t, bool> p = m_table_data.TryPushBack(temp);
		if(!p.second) {
			m_table_data[p.first].value += value;
		}
	}

	// Add to dense matrix.
	template<typename F1, bool ROWMAJOR1>
	void DenseAdd(DenseViewBase<F1, ROWMAJOR1> output) const {
		for(size_t outer = 0; outer < m_outer_size; ++outer) {
			const BulkEntry *bulk = m_bulk_data.get() + m_bulk_size * outer;
			for(size_t i = 0; i < m_bulk_size; ++i) {
				const BulkEntry &entry = bulk[i];
				if(entry.inner == INDEX_NONE)
					break;
				size_t row = Row(outer, entry.inner), col = Col(outer, entry.inner);
				output(row, col) += entry.value;
				if(SYMMETRIC) {
					output(col, row) += entry.value;
				}
			}
		}
		for(size_t i = 0; i < m_table_data.GetSize(); ++i) {
			const TableEntry &entry = m_table_data[i];
			size_t row = Row(entry.outer, entry.inner), col = Col(entry.outer, entry.inner);
			output(row, col) += entry.value;
			if(SYMMETRIC) {
				output(col, row) += entry.value;
			}
		}
	}

	// Subtract from dense matrix.
	template<typename F1, bool ROWMAJOR1>
	void DenseSubtract(DenseViewBase<F1, ROWMAJOR1> output) const {
		for(size_t outer = 0; outer < m_outer_size; ++outer) {
			const BulkEntry *bulk = m_bulk_data.get() + m_bulk_size * outer;
			for(size_t i = 0; i < m_bulk_size; ++i) {
				const BulkEntry &entry = bulk[i];
				if(entry.inner == INDEX_NONE)
					break;
				size_t row = Row(outer, entry.inner), col = Col(outer, entry.inner);
				output(row, col) -= entry.value;
				if(SYMMETRIC) {
					output(col, row) -= entry.value;
				}
			}
		}
		for(size_t i = 0; i < m_table_data.GetSize(); ++i) {
			const TableEntry &entry = m_table_data[i];
			size_t row = Row(entry.outer, entry.inner), col = Col(entry.outer, entry.inner);
			output(row, col) -= entry.value;
			if(SYMMETRIC) {
				output(col, row) -= entry.value;
			}
		}
	}

	// Left or right matrix multiplication + addition.
	template<typename F1, typename F2, bool ROWMAJOR1, bool ROWMAJOR2, bool RIGHT>
	void MultiplyAddBase(DenseViewBase<F1, ROWMAJOR1> output, DenseViewBase<F2, ROWMAJOR2> input, size_t num) const {
		for(size_t outer = 0; outer < m_outer_size; ++outer) {
			const BulkEntry *bulk = m_bulk_data.get() + m_bulk_size * outer;
			for(size_t i = 0; i < m_bulk_size; ++i) {
				const BulkEntry &entry = bulk[i];
				if(entry.inner == INDEX_NONE)
					break;
				size_t row = Row(outer, entry.inner), col = Col(outer, entry.inner);
				F value = entry.value;
				if(RIGHT) {
					for(size_t j = 0; j < num; ++j) {
						output(j, col) += value * input(j, row);
						if(SYMMETRIC) {
							output(j, row) += value * input(j, col);
						}
					}
				} else {
					for(size_t j = 0; j < num; ++j) {
						output(row, j) += value * input(col, j);
						if(SYMMETRIC) {
							output(col, j) += value * input(row, j);
						}
					}
				}
			}
		}
		for(size_t i = 0; i < m_table_data.GetSize(); ++i) {
			const TableEntry &entry = m_table_data[i];
			size_t row = Row(entry.outer, entry.inner), col = Col(entry.outer, entry.inner);
			F value = entry.value;
			if(RIGHT) {
				for(size_t j = 0; j < num; ++j) {
					output(j, col) += value * input(j, row);
					if(SYMMETRIC) {
						output(j, row) += value * input(j, col);
					}
				}
			} else {
				for(size_t j = 0; j < num; ++j) {
					output(row, j) += value * input(col, j);
					if(SYMMETRIC) {
						output(col, j) += value * input(row, j);
					}
				}
			}
		}
	}

	// Left or right matrix multiplication.
	template<typename F1, typename F2, bool ROWMAJOR1, bool ROWMAJOR2, bool RIGHT>
	void MultiplyBase(DenseViewBase<F1, ROWMAJOR1> output, DenseViewBase<F2, ROWMAJOR2> input, size_t num) const {
		if(RIGHT) {
			output.SetZero(num, GetCols());
		} else {
			output.SetZero(GetRows(), num);
		}
		MultiplyAddBase<F1, F2, ROWMAJOR1, ROWMAJOR2, RIGHT>(output, input, num);
	}

	// Calculates output += A * input (convenience function).
	// size of output = (rows, num), size of input = (cols, num).
	template<typename F1, typename F2, bool ROWMAJOR1, bool ROWMAJOR2>
	void LeftMultiplyAdd(DenseViewBase<F1, ROWMAJOR1> output, DenseViewBase<F2, ROWMAJOR2> input, size_t num) const {
		MultiplyAddBase<F1, F2, ROWMAJOR1, ROWMAJOR2, false>(output, input, num);
	}

	// Calculates output += input * A (convenience function).
	// size of output = (num, cols), size of input = (num, rows).
	template<typename F1, typename F2, bool ROWMAJOR1, bool ROWMAJOR2>
	void RightMultiplyAdd(DenseViewBase<F1, ROWMAJOR1> output, DenseViewBase<F2, ROWMAJOR2> input, size_t num) const {
		MultiplyAddBase<F1, F2, ROWMAJOR1, ROWMAJOR2, true>(output, input, num);
	}

	// Calculates output = A * input (convenience function).
	// size of output = (rows, num), size of input = (cols, num).
	template<typename F1, typename F2, bool ROWMAJOR1, bool ROWMAJOR2>
	void LeftMultiply(DenseViewBase<F1, ROWMAJOR1> output, DenseViewBase<F2, ROWMAJOR2> input, size_t num) const {
		MultiplyBase<F1, F2, ROWMAJOR1, ROWMAJOR2, false>(output, input, num);
	}

	// Calculates output = input * A (convenience function).
	// size of output = (num, cols), size of input = (num, rows).
	template<typename F1, typename F2, bool ROWMAJOR1, bool ROWMAJOR2>
	void RightMultiply(DenseViewBase<F1, ROWMAJOR1> output, DenseViewBase<F2, ROWMAJOR2> input, size_t num) const {
		MultiplyBase<F1, F2, ROWMAJOR1, ROWMAJOR2, true>(output, input, num);
	}

	// Convert to dense matrix.
	template<typename F1, bool ROWMAJOR1>
	void ToDense(DenseViewBase<F1, ROWMAJOR1> output) const {
		output.SetZero(GetRows(), GetCols());
		DenseAdd(output);
	}

	// Convert to compressed sparse row/column format.
	// For symmetric matrices, only the lower or upper part will be stored, but the diagonal values will be doubled.
	template<typename I, typename F1, bool TRANSPOSE>
	void ToCompressedSparse(I *offsets, I *indices, F1 *values) const {

		// count coefficients per row/column
		std::fill_n(offsets, (TRANSPOSE)? m_inner_size : m_outer_size, 0);
		for(size_t outer = 0; outer < m_outer_size; ++outer) {
			const BulkEntry *bulk = m_bulk_data.get() + m_bulk_size * outer;
			for(size_t i = 0; i < m_bulk_size; ++i) {
				const BulkEntry &entry = bulk[i];
				if(entry.inner == INDEX_NONE)
					break;
				++offsets[(TRANSPOSE)? entry.inner : outer];
			}
		}
		for(size_t i = 0; i < m_table_data.GetSize(); ++i) {
			const TableEntry &entry = m_table_data[i];
			++offsets[(TRANSPOSE)? entry.inner : entry.outer];
		}

		// calculate offsets
		I total = 0;
		for(size_t i = 0; i < ((TRANSPOSE)? m_inner_size : m_outer_size); ++i) {
			total += offsets[i];
			offsets[i] = total;
		}
		offsets[(TRANSPOSE)? m_inner_size : m_outer_size] = total;
		assert((size_t) total == GetCoefficients());

		// calculate indices and values
		for(size_t outer = 0; outer < m_outer_size; ++outer) {
			const BulkEntry *bulk = m_bulk_data.get() + m_bulk_size * outer;
			for(size_t i = 0; i < m_bulk_size; ++i) {
				const BulkEntry &entry = bulk[i];
				if(entry.inner == INDEX_NONE)
					break;
				size_t j = --offsets[(TRANSPOSE)? entry.inner : outer];
				indices[j] = I((TRANSPOSE)? outer : entry.inner);
				values[j] = F1((outer == entry.inner)? entry.value + entry.value : entry.value);
			}
		}
		for(size_t i = 0; i < m_table_data.GetSize(); ++i) {
			const TableEntry &entry = m_table_data[i];
			size_t j = --offsets[(TRANSPOSE)? entry.inner : entry.outer];
			indices[j] = I((TRANSPOSE)? entry.outer : entry.inner);
			values[j] = F1((entry.outer == entry.inner)? entry.value + entry.value : entry.value);
		}

		// sort by index within each column
		for(size_t i = 0; i < ((TRANSPOSE)? m_inner_size : m_outer_size); ++i) {
			I begin = offsets[i], end = offsets[i + 1];
			InsertionSortPairs(indices + begin, values + begin, end - begin);
		}

	}

	// Convert to compressed sparse column format (convenience function).
	template<typename I, typename F1>
	void ToCSC(I *offsets, I *indices, F1 *values) const {
		ToCompressedSparse<I, F1, ROWMAJOR>(offsets, indices, values);
	}

	// Convert to compressed sparse row format (convenience function).
	template<typename I, typename F1>
	void ToCSR(I *offsets, I *indices, F1 *values) const {
		ToCompressedSparse<I, F1, !ROWMAJOR>(offsets, indices, values);
	}

public:
	inline size_t GetOuterSize() const { return m_outer_size; }
	inline size_t GetInnerSize() const { return m_inner_size; }
	inline size_t GetRows() const { return Row(m_outer_size, m_inner_size); }
	inline size_t GetCols() const { return Col(m_outer_size, m_inner_size); }
	inline size_t GetCoefficients() const { return m_bulk_coefficients + m_table_data.GetSize(); }

};

template<typename F> using SparseMatrixC   = SparseMatrixBase<F, false, false, false>;
template<typename F> using SparseMatrixCSL = SparseMatrixBase<F, false, true , false>;
template<typename F> using SparseMatrixCSU = SparseMatrixBase<F, false, true , true >;
template<typename F> using SparseMatrixR   = SparseMatrixBase<F, true , false, false>;
template<typename F> using SparseMatrixRSL = SparseMatrixBase<F, true , true , false>;
template<typename F> using SparseMatrixRSU = SparseMatrixBase<F, true , true , true >;

template<typename F, bool ROWMAJOR, bool SYMMETRIC, bool UPPER>
class SparseBlockMatrixBase {};

template<typename F, bool ROWMAJOR, bool UPPER>
class SparseBlockMatrixBase<F, ROWMAJOR, false, UPPER> {

private:
	SparseMatrixBase<F, ROWMAJOR, false, false> m_matrix_a, m_matrix_b, m_matrix_c, m_matrix_d;

public:
	void Free() {
		m_matrix_a.Free();
		m_matrix_b.Free();
		m_matrix_c.Free();
		m_matrix_d.Free();
	}

	void Reset(size_t rows1, size_t rows2, size_t cols1, size_t cols2, size_t bulk_size) {
		assert(rows1 < INDEX_OFFSET);
		assert(rows2 < INDEX_OFFSET - 1);
		assert(cols1 < INDEX_OFFSET);
		assert(cols2 < INDEX_OFFSET - 1);
		m_matrix_a.Reset(rows1, cols1, bulk_size);
		m_matrix_b.Reset(rows1, cols2, bulk_size);
		m_matrix_c.Reset(rows2, cols1, bulk_size);
		m_matrix_d.Reset(rows2, cols2, bulk_size);
	}

	// Insert a new coefficient. If it already exists, the new value will be added to the existing one.
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

	inline const SparseMatrixBase<F, ROWMAJOR, false, false>& GetMatrixA() const { return m_matrix_a; }
	inline const SparseMatrixBase<F, ROWMAJOR, false, false>& GetMatrixB() const { return m_matrix_b; }
	inline const SparseMatrixBase<F, ROWMAJOR, false, false>& GetMatrixC() const { return m_matrix_c; }
	inline const SparseMatrixBase<F, ROWMAJOR, false, false>& GetMatrixD() const { return m_matrix_d; }

};

template<typename F, bool ROWMAJOR, bool UPPER>
class SparseBlockMatrixBase<F, ROWMAJOR, true, UPPER> {

private:
	SparseMatrixBase<F, ROWMAJOR, true , UPPER> m_matrix_a;
	SparseMatrixBase<F, ROWMAJOR, false, false> m_matrix_bc;
	SparseMatrixBase<F, ROWMAJOR, true , UPPER> m_matrix_d;

public:
	void Free() {
		m_matrix_a.Free();
		m_matrix_bc.Free();
		m_matrix_d.Free();
	}

	void Reset(size_t rows1, size_t rows2, size_t cols1, size_t cols2, size_t bulk_size) {
		assert(rows1 < INDEX_OFFSET);
		assert(rows2 < INDEX_OFFSET - 1);
		assert(cols1 < INDEX_OFFSET);
		assert(cols2 < INDEX_OFFSET - 1);
		assert(rows1 == cols1);
		assert(rows2 == cols2);
		m_matrix_a.Reset(rows1, cols1, bulk_size);
		if(UPPER) {
			m_matrix_bc.Reset(rows1, cols2, bulk_size); // use B
		} else {
			m_matrix_bc.Reset(rows2, cols1, bulk_size); // use C
		}
		m_matrix_d.Reset(rows2, cols2, bulk_size);
	}

	// Insert a new coefficient. If it already exists, the new value will be added to the existing one.
	void Insert(size_t row, size_t col, F value) {
		assert(row < GetRows1() || (row >= INDEX_OFFSET && row < INDEX_OFFSET + GetRows2()));
		assert(col < GetCols1() || (col >= INDEX_OFFSET && col < INDEX_OFFSET + GetCols2()));
		if((UPPER)? (row > col) : (row < col)) {
			std::swap(row, col);
		}
		if(((UPPER)? col : row) < INDEX_OFFSET) {
			m_matrix_a.template Insert<false>(row, col, value);
		} else if(((UPPER)? row : col) < INDEX_OFFSET) {
			if(UPPER) {
				m_matrix_bc.template Insert<false>(row, col - INDEX_OFFSET, value);
			} else {
				m_matrix_bc.template Insert<false>(row - INDEX_OFFSET, col, value);
			}
		} else {
			m_matrix_d.template Insert<false>(row - INDEX_OFFSET, col - INDEX_OFFSET, value);
		}
	}

	// Calculates rhs = -B * fixed.
	// size of rhs = (rows1, num), size of fixed = (cols2, num).
	template<typename F1, typename F2, bool ROWMAJOR1, bool ROWMAJOR2>
	void CalculateRHS(DenseViewBase<F1, ROWMAJOR1> rhs, DenseViewBase<F2, ROWMAJOR2> fixed, size_t num) const {
		if(UPPER) {
			m_matrix_bc.LeftMultiply(rhs, fixed, num); // use B
		} else {
			m_matrix_bc.RightMultiply(rhs.TransposedView(), fixed.TransposedView(), num); // use C
		}
		rhs.Negate(GetRows1(), num);
	}

	// Calculates residual = fixed' * C * solution + fixed' * D * fixed.
	// size of residual = (num, num), size of solution = (cols1, num), size of fixed = (cols2, num).
	template<typename F1, typename F2, typename F3, bool ROWMAJOR1, bool ROWMAJOR2, bool ROWMAJOR3>
	void CalculateResidual(DenseViewBase<F1, ROWMAJOR1> residual, DenseViewBase<F2, ROWMAJOR2> solution, DenseViewBase<F3, ROWMAJOR3> fixed, size_t num) const {
		std::unique_ptr<F3[]> temp(new F3[GetRows2() * num]);
		DenseViewBase<F3, ROWMAJOR3> temp_view(temp.get(), (ROWMAJOR3)? num : GetRows2());
		if(UPPER) {
			m_matrix_bc.RightMultiply(temp_view.TransposedView(), solution.TransposedView(), num); // use B
		} else {
			m_matrix_bc.LeftMultiply(temp_view, solution, num); // use C
		}
		m_matrix_d.LeftMultiplyAdd(temp_view, fixed, num);
		for(size_t i = 0; i < num; ++i) {
			for(size_t j = 0; j < num; ++j) {
				F2 sum = F2();
				for(size_t k = 0; k < GetRows2(); ++k) {
					sum += fixed(k, i) * temp_view(k, j);
				}
				residual(i, j) = sum;
			}
		}
	}

	// Calculates output = solution' * A * solution + solution' * B * fixed + fixed' * C * solution + fixed' * D * fixed
	// size of output = (num, num), size of solution = (cols1, num), size of fixed = (cols2, num).
	template<typename F1, typename F2, typename F3, bool ROWMAJOR1, bool ROWMAJOR2, bool ROWMAJOR3>
	void CalculateQuadratic(DenseViewBase<F1, ROWMAJOR1> output, DenseViewBase<F2, ROWMAJOR2> solution, DenseViewBase<F3, ROWMAJOR3> fixed, size_t num) {
		std::unique_ptr<F2[]> temp1(new F2[GetRows1() * num]);
		std::unique_ptr<F3[]> temp2(new F3[GetRows2() * num]);
		DenseViewBase<F2, ROWMAJOR2> temp1_view(temp1.get(), (ROWMAJOR2)? num : GetRows1());
		DenseViewBase<F3, ROWMAJOR3> temp2_view(temp2.get(), (ROWMAJOR3)? num : GetRows2());
		m_matrix_a.LeftMultiply(temp1_view, solution, num);
		if(UPPER) {
			m_matrix_bc.LeftMultiplyAdd(temp1_view, fixed, num); // use B
			m_matrix_bc.RightMultiply(temp2_view.TransposedView(), solution.TransposedView(), num); // use B
		} else {
			m_matrix_bc.RightMultiplyAdd(temp1_view.TransposedView(), fixed.TransposedView(), num); // use C
			m_matrix_bc.LeftMultiply(temp2_view, solution, num); // use C
		}
		m_matrix_d.LeftMultiplyAdd(temp2_view, fixed, num);
		for(size_t i = 0; i < num; ++i) {
			for(size_t j = 0; j < num; ++j) {
				F2 sum = F2();
				for(size_t k = 0; k < GetRows1(); ++k) {
					sum += solution(k, i) * temp1_view(k, j);
				}
				for(size_t k = 0; k < GetRows2(); ++k) {
					sum += fixed(k, i) * temp2_view(k, j);
				}
				output(i, j) = sum;
			}
		}
	}

public:
	inline size_t GetRows1() const { return m_matrix_a.GetRows(); }
	inline size_t GetRows2() const { return m_matrix_d.GetRows(); }
	inline size_t GetCols1() const { return m_matrix_a.GetCols(); }
	inline size_t GetCols2() const { return m_matrix_d.GetCols(); }

	inline const SparseMatrixBase<F, ROWMAJOR, true , UPPER>& GetMatrixA()  const { return m_matrix_a;  }
	inline const SparseMatrixBase<F, ROWMAJOR, false, false>& GetMatrixBC() const { return m_matrix_bc; }
	inline const SparseMatrixBase<F, ROWMAJOR, true , UPPER>& GetMatrixD()  const { return m_matrix_d;  }

};

template<typename F> using SparseBlockMatrixC   = SparseBlockMatrixBase<F, false, false, false>;
template<typename F> using SparseBlockMatrixCSL = SparseBlockMatrixBase<F, false, true , false>;
template<typename F> using SparseBlockMatrixCSU = SparseBlockMatrixBase<F, false, true , true >;
template<typename F> using SparseBlockMatrixR   = SparseBlockMatrixBase<F, true , false, false>;
template<typename F> using SparseBlockMatrixRSL = SparseBlockMatrixBase<F, true , true , false>;
template<typename F> using SparseBlockMatrixRSU = SparseBlockMatrixBase<F, true , true , true >;
