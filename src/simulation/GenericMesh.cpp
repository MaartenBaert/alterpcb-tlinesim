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
	if((size_t) modes.rows() != GetPortCount())
		throw std::runtime_error(MakeString("GenericMesh error: Expected mode matrix with ", GetPortCount(), " rows, got ", modes.rows(), " instead."));
	if((size_t) modes.cols() == 0)
		throw std::runtime_error(MakeString("GenericMesh error: Expected mode matrix with at least one column, got ", modes.cols(), " instead."));
	if((size_t) modes.cols() >= GetPortCount())
		throw std::runtime_error(MakeString("GenericMesh error: Expected mode matrix with less than ", GetPortCount(), " columns, got ", modes.cols(), " instead."));
	m_solved = false;
	m_modes = modes;
	m_frequency = frequency;
	DoSolve();
	m_solved = true;
}

void GenericMesh::Cleanup() {
	DoCleanup();
}
