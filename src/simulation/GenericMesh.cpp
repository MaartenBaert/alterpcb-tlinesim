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

#include "GenericMesh.h"

#include "EigenSparse.h"
#include "StringHelper.h"

#include <iostream>

GenericMesh::GenericMesh() {
	m_initialized = false;
	m_solved = false;
	m_frequency = 0.0;
}

GenericMesh::~GenericMesh() {
	// nothing
}

void GenericMesh::Initialize() {
	if(m_initialized)
		throw std::runtime_error("GenericMesh error: The mesh has already been initialized.");
	DoInitialize();
	m_initialized = true;
}

void GenericMesh::Solve(const Eigen::MatrixXr &modes, real_t frequency) {
	if(!m_initialized)
		throw std::runtime_error("GenericMesh error: The mesh must be initialized first.");
	if((size_t) modes.rows() != GetFixedVariableCount())
		throw std::runtime_error(MakeString("GenericMesh error: Expected mode matrix with ", GetFixedVariableCount(), " rows, got ", modes.rows(), " instead."));
	if((size_t) modes.cols() == 0)
		throw std::runtime_error(MakeString("GenericMesh error: Expected mode matrix with at least one column, got ", modes.cols(), " instead."));
	if((size_t) modes.cols() >= GetFixedVariableCount())
		throw std::runtime_error(MakeString("GenericMesh error: Expected mode matrix with less than ", GetFixedVariableCount(), " columns, got ", modes.cols(), " instead."));
	m_solved = false;
	m_modes = modes;
	m_frequency = frequency;
	DoSolve();
	SolveEigenModes();
	m_solved = true;
}

void GenericMesh::Cleanup() {
	DoCleanup();
}

void GenericMesh::SolveEigenModes() {

	std::cerr << "inductance =\n" << m_inductance_matrix << std::endl;
	std::cerr << "capacitance =\n" << m_capacitance_matrix << std::endl;
	std::cerr << "resistance =\n" << m_resistance_matrix << std::endl;
	std::cerr << "conductance =\n" << m_conductance_matrix << std::endl;
	std::cerr << std::endl;

	// convert to impedance and admittance matrices
	real_t solution_omega = 2.0 * M_PI * m_frequency;
	Eigen::MatrixXc impedance(m_modes.cols(), m_modes.cols()), admittance(m_modes.cols(), m_modes.cols());
	impedance.real() = m_resistance_matrix;
	impedance.imag() = solution_omega * m_inductance_matrix;
	admittance.real() = m_conductance_matrix;
	admittance.imag() = solution_omega * m_capacitance_matrix;

	// calculate eigenmodes
	Eigen::MatrixXc matrix_zy = impedance * admittance;
	Eigen::ComplexEigenSolver<Eigen::MatrixXc> eigensolver(matrix_zy);
	if(eigensolver.info() != Eigen::Success)
		throw std::runtime_error("Eigenmode decomposition failed!");
	auto &eigval = eigensolver.eigenvalues();
	auto &eigvec = eigensolver.eigenvectors();

	std::cerr << "eigenvalues =\n" << eigval << std::endl;
	std::cerr << "eigenvectors =\n" << eigvec << std::endl;
	std::cerr << std::endl;

	// make sure we get the correct complex square root (positive imaginary part)
	Eigen::VectorXc eigval_sqrt = (-eigval).cwiseSqrt() * complex_t(0.0, 1.0);

	// calculate impedances, admittance and (approximated) characteristic impedance
	Eigen::MatrixXc matrix_zy_invsqrt = eigvec * eigval_sqrt.cwiseInverse().asDiagonal() * eigvec.inverse();
	m_characteristic_impedance_matrix = matrix_zy_invsqrt * impedance;
	Eigen::VectorXc diag_impedance = m_characteristic_impedance_matrix.diagonal();
	Eigen::VectorXc diag_admittance = m_characteristic_impedance_matrix.inverse().diagonal();
	m_characteristic_impedances = (diag_impedance.array() / diag_admittance.array()).sqrt().matrix();

	// calculate propagation constants
	m_propagation_constants = (-matrix_zy.diagonal()).cwiseSqrt() * complex_t(0.0, 1.0);

	std::cerr << "m_characteristic_impedance_matrix =\n" << m_characteristic_impedance_matrix << std::endl;
	std::cerr << "m_characteristic_impedances =\n" << m_characteristic_impedances << std::endl;
	std::cerr << "m_propagation_constants =\n" << m_propagation_constants << std::endl;
	std::cerr << std::endl;

	// reorder eigenmodes to match user-provided modes and calculate eigenmode propagation constants
	m_eigenmodes.resize(m_modes.cols(), m_modes.cols());
	m_eigenmode_propagation_constants.resize(m_modes.cols());
	std::vector<size_t> modemap((size_t) m_modes.cols());
	for(size_t i = 0; i < (size_t) m_modes.cols(); ++i) {
		modemap[i] = i;
	}
	for(size_t i = 0; i < (size_t) m_modes.cols(); ++i) {
		size_t best_index = i;
		complex_t best_value = complex_t(0.0, 0.0);
		for(size_t j = i; j < (size_t) m_modes.cols(); ++j) {
			complex_t value = eigvec((Eigen::Index) i, (Eigen::Index) modemap[j]);
			if(std::norm(value) > std::norm(best_value)) {
				best_index = j;
				best_value = value;
			}
		}
		std::swap(modemap[i], modemap[best_index]);
		m_eigenmodes.col((Eigen::Index) i) = eigvec.col((Eigen::Index) modemap[i]); // / best_value;
		m_eigenmode_propagation_constants[(Eigen::Index) i] = eigval_sqrt((Eigen::Index) modemap[i]);
	}

	std::cerr << "m_eigenmodes =\n" << m_eigenmodes << std::endl;
	std::cerr << "m_eigenmode_propagation_constants =\n" << m_eigenmode_propagation_constants << std::endl;
	std::cerr << std::endl;

}
