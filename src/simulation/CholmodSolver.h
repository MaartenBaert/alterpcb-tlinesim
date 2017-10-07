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

#include <complex>

#include <cholmod.h>

// Wrapper for cholmod_common.
class CholmodSolver {

private:
	static CholmodSolver *s_instance;

private:
	cholmod_common m_common;

public:
	CholmodSolver();
	~CholmodSolver();

public:
	inline static cholmod_common* GetCommon() { assert(s_instance != NULL); return &s_instance->m_common; }

};

// Wrapper for cholmod_sparse.
// Note: CHOLMOD stores sparse matrices in CSC format.
class CholmodSparseMatrix {

private:
	cholmod_sparse *m_sparse;

public:
	CholmodSparseMatrix();
	~CholmodSparseMatrix();

	void Reset();
	void Reset(size_t rows, size_t cols, size_t coefs, bool sorted, bool packed, int stype, int xtype);
	void Reset(SymmetricSparseMatrix<real_t> &matrix);
	void Reset(SymmetricSparseMatrix<complex_t> &matrix);

public:
	inline cholmod_sparse* GetSparse() { return m_sparse; }

};

// Wrapper for cholmod_dense.
// Note: CHOLMOD stores dense matrices in column-major (Fortran) order.
class CholmodDenseMatrix {

private:
	cholmod_dense *m_dense;

public:
	CholmodDenseMatrix();
	~CholmodDenseMatrix();

	void Reset();
	//void Reset(cholmod_dense *dense);
	void Reset(size_t rows, size_t cols, size_t stride, int xtype);
	void ResetReal(size_t rows, size_t cols);
	void ResetComplex(size_t rows, size_t cols);

public:
	inline cholmod_dense*& GetDense() { return m_dense; }
	inline bool IsValid() const { return m_dense != NULL; }
	inline size_t GetRows() const { assert(m_dense != NULL); return m_dense->nrow; }
	inline size_t GetCols() const { assert(m_dense != NULL); return m_dense->ncol; }
	inline size_t GetStride() const { assert(m_dense != NULL); return m_dense->d; }
	inline       real_t* GetRealData()       { assert(m_dense != NULL); assert(m_dense->xtype == CHOLMOD_REAL); return (real_t*) m_dense->x; }
	inline const real_t* GetRealData() const { assert(m_dense != NULL); assert(m_dense->xtype == CHOLMOD_REAL); return (real_t*) m_dense->x; }
	inline       complex_t* GetComplexData()       { assert(m_dense != NULL); assert(m_dense->xtype == CHOLMOD_COMPLEX); return (complex_t*) m_dense->x; }
	inline const complex_t* GetComplexData() const { assert(m_dense != NULL); assert(m_dense->xtype == CHOLMOD_COMPLEX); return (complex_t*) m_dense->x; }

};

class CholmodFactorization {

private:
	cholmod_factor *m_factor;

public:
	CholmodFactorization();
	~CholmodFactorization();

	void Reset();
	void Analyze(const CholmodSparseMatrix &matrix);
	void Factorize(const CholmodSparseMatrix &matrix);
	void Solve(CholmodDenseMatrix &solution, const CholmodDenseMatrix &rhs);
	void Solve(CholmodDenseMatrix &solution, const CholmodDenseMatrix &rhs, CholmodDenseMatrix &ws1, CholmodDenseMatrix &ws2);

public:
	inline cholmod_factor* GetFactor() { return m_factor; }

};
