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
				size_t j = (size_t) --offsets[(TRANSPOSE)? entry.inner : outer];
				indices[j] = I((TRANSPOSE)? outer : entry.inner);
				values[j] = (SYMMETRIC && outer == entry.inner)? entry.value + entry.value : entry.value;
			}
		}
		for(size_t i = 0; i < m_table_data.GetSize(); ++i) {
			const TableEntry &entry = m_table_data[i];
			size_t j = (size_t) --offsets[(TRANSPOSE)? entry.inner : entry.outer];
			indices[j] = I((TRANSPOSE)? entry.outer : entry.inner);
			values[j] = (SYMMETRIC && entry.outer == entry.inner)? entry.value + entry.value : entry.value;
		}

		// sort by index within each column
		for(size_t i = 0; i < ((TRANSPOSE)? m_inner_size : m_outer_size); ++i) {
			I begin = offsets[i], end = offsets[i + 1];
			InsertionSortPairs(indices + begin, values + begin, (size_t) (end - begin));
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

	// Convert to Eigen sparse matrix.
	template<class EigenSparseMatrix>
	void ToEigen(EigenSparseMatrix &output) const {
		output.resize((Eigen::Index) GetRows(), (Eigen::Index) GetCols());
		output.resizeNonZeros((Eigen::Index) GetCoefficients());
		if(output.IsRowMajor) {
			ToCSR(output.outerIndexPtr(), output.innerIndexPtr(), output.valuePtr());
		} else {
			ToCSC(output.outerIndexPtr(), output.innerIndexPtr(), output.valuePtr());
		}
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
