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

#include "CholmodSolver.h"

static_assert(sizeof(real_t) == sizeof(double), "real_t has incorrect size!");
static_assert(sizeof(complex_t) == 2 * sizeof(double), "complex_t has incorrect size!");

CholmodSolver *CholmodSolver::s_instance = NULL;

CholmodSolver::CholmodSolver() {
	if(!cholmod_start(&m_common))
		throw std::runtime_error("Cholmod start failed.");
	m_common.nmethods = 2;
	s_instance = this;
}

CholmodSolver::~CholmodSolver() {
	s_instance = NULL;
	//cholmod_print_common("CholmodSolver", &m_common);
	cholmod_finish(&m_common);
}

CholmodSparseMatrix::CholmodSparseMatrix() {
	m_sparse = NULL;
}

CholmodSparseMatrix::~CholmodSparseMatrix() {
	cholmod_free_sparse(&m_sparse, CholmodSolver::GetCommon());
}

void CholmodSparseMatrix::Reset() {
	cholmod_free_sparse(&m_sparse, CholmodSolver::GetCommon());
}

void CholmodSparseMatrix::Reset(size_t rows, size_t cols, size_t coefs, bool sorted, bool packed, int stype, int xtype) {
	if(m_sparse == NULL || m_sparse->nrow != rows || m_sparse->ncol != cols || m_sparse->nzmax != coefs
			|| m_sparse->sorted != sorted || m_sparse->packed != packed || m_sparse->stype != stype || m_sparse->xtype != xtype) {
		cholmod_free_sparse(&m_sparse, CholmodSolver::GetCommon());
		m_sparse = cholmod_allocate_sparse(rows, cols, coefs, sorted, packed, stype, xtype, CholmodSolver::GetCommon());
		if(m_sparse == NULL)
			throw std::runtime_error("Cholmod sparse matrix allocation failed.");
	}
}

void CholmodSparseMatrix::Reset(SymmetricSparseMatrix<real_t> &matrix) {
	Reset(matrix.GetSize1(), matrix.GetSize1(), matrix.GetCoefficientsA(), false, true, 1, CHOLMOD_REAL);
	matrix.MakeCscA((int*) m_sparse->p, (int*) m_sparse->i, (real_t*) m_sparse->x);
}

void CholmodSparseMatrix::Reset(SymmetricSparseMatrix<complex_t> &matrix) {
	Reset(matrix.GetSize1(), matrix.GetSize1(), matrix.GetCoefficientsA(), false, true, 1, CHOLMOD_COMPLEX);
	matrix.MakeCscA((int*) m_sparse->p, (int*) m_sparse->i, (complex_t*) m_sparse->x);
}

CholmodDenseMatrix::CholmodDenseMatrix() {
	m_dense = NULL;
}

CholmodDenseMatrix::~CholmodDenseMatrix() {
	cholmod_free_dense(&m_dense, CholmodSolver::GetCommon());
}

void CholmodDenseMatrix::Reset() {
	cholmod_free_dense(&m_dense, CholmodSolver::GetCommon());
}

/*void CholmodDenseMatrix::Reset(cholmod_dense *dense) {
	cholmod_free_dense(&m_dense, CholmodSolver::GetCommon());
	m_dense = dense;
}*/

void CholmodDenseMatrix::Reset(size_t rows, size_t cols, size_t stride, int xtype) {
	if(m_dense == NULL || m_dense->nrow != rows || m_dense->ncol != cols || m_dense->d != stride || m_dense->xtype != xtype) {
		cholmod_free_dense(&m_dense, CholmodSolver::GetCommon());
		m_dense = cholmod_allocate_dense(rows, cols, stride, xtype, CholmodSolver::GetCommon());
		if(m_dense == NULL)
			throw std::runtime_error("Cholmod dense matrix allocation failed.");
	}
}

void CholmodDenseMatrix::ResetReal(size_t rows, size_t cols) {
	Reset(rows, cols, rows, CHOLMOD_REAL);
}

void CholmodDenseMatrix::ResetComplex(size_t rows, size_t cols) {
	Reset(rows, cols, rows, CHOLMOD_COMPLEX);
}

CholmodFactorization::CholmodFactorization() {
	m_factor = NULL;
}

CholmodFactorization::~CholmodFactorization() {
	cholmod_free_factor(&m_factor, CholmodSolver::GetCommon());
}

void CholmodFactorization::Reset() {
	cholmod_free_factor(&m_factor, CholmodSolver::GetCommon());
}

void CholmodFactorization::Analyze(const CholmodSparseMatrix &matrix) {
	assert(const_cast<CholmodSparseMatrix&>(matrix).GetSparse() != NULL);
	cholmod_free_factor(&m_factor, CholmodSolver::GetCommon());
	m_factor = cholmod_analyze(const_cast<CholmodSparseMatrix&>(matrix).GetSparse(), CholmodSolver::GetCommon());
	if(m_factor == NULL)
		throw std::runtime_error("Cholmod sparse matrix analysis failed.");
}

void CholmodFactorization::Factorize(const CholmodSparseMatrix &matrix) {
	assert(const_cast<CholmodSparseMatrix&>(matrix).GetSparse() != NULL);
	if(m_factor == NULL) {
		Analyze(matrix);
	}
	if(!cholmod_factorize(const_cast<CholmodSparseMatrix&>(matrix).GetSparse(), m_factor, CholmodSolver::GetCommon())) {
		throw std::runtime_error("Cholmod sparse matrix factorization failed.");
	}
	//cholmod_print_common("CholmodSolver", CholmodSolver::GetCommon());
}

void CholmodFactorization::Solve(CholmodDenseMatrix &solution, const CholmodDenseMatrix &rhs) {
	if(!cholmod_solve2(CHOLMOD_A, m_factor, const_cast<CholmodDenseMatrix&>(rhs).GetDense(), NULL, &solution.GetDense(), NULL, NULL, NULL, CholmodSolver::GetCommon()))
		throw std::runtime_error("Cholmod sparse matrix solve failed.");
}

void CholmodFactorization::Solve(CholmodDenseMatrix &solution, const CholmodDenseMatrix &rhs, CholmodDenseMatrix &ws1, CholmodDenseMatrix &ws2) {
	if(!cholmod_solve2(CHOLMOD_A, m_factor, const_cast<CholmodDenseMatrix&>(rhs).GetDense(), NULL, &solution.GetDense(), NULL, &ws1.GetDense(), &ws2.GetDense(), CholmodSolver::GetCommon()))
		throw std::runtime_error("Cholmod sparse matrix solve failed.");
}
