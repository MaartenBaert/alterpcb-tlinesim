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
#include "Vector.h"

#include <vector>

enum MeshImageType {
	MESHIMAGETYPE_MESH,
	MESHIMAGETYPE_EFIELD,
	MESHIMAGETYPE_EFIELD_X,
	MESHIMAGETYPE_EFIELD_Y,
	MESHIMAGETYPE_EFIELD_Z,
	MESHIMAGETYPE_MFIELD,
	MESHIMAGETYPE_MFIELD_X,
	MESHIMAGETYPE_MFIELD_Y,
	MESHIMAGETYPE_MFIELD_Z,
	MESHIMAGETYPE_EPOT,
	MESHIMAGETYPE_MPOT,
	MESHIMAGETYPE_ENERGY,
	MESHIMAGETYPE_CURRENT,
};

enum SolverType {
	SOLVERTYPE_QUASISTATIC,
	SOLVERTYPE_FULLWAVE,
};

enum ElementType {
	ELEMENTTYPE_LINEAR,
	ELEMENTTYPE_QUADRATIC,
};

class GenericMesh {

private:
	bool m_initialized, m_solved;
	Eigen::MatrixXr m_modes;
	real_t m_frequency;

protected:
	Eigen::MatrixXr m_inductance_matrix, m_capacitance_matrix, m_resistance_matrix, m_conductance_matrix;
	Eigen::MatrixXc m_characteristic_impedance_matrix;
	Eigen::VectorXc m_characteristic_impedances, m_propagation_constants;
	Eigen::MatrixXc m_eigenmodes;
	Eigen::VectorXc m_eigenmode_propagation_constants;

public:
	GenericMesh();
	virtual ~GenericMesh();

	void Initialize();
	void Solve(const Eigen::MatrixXr &modes, real_t frequency);
	void Cleanup();

	virtual Box2D GetWorldBox2D() = 0;
	virtual Box2D GetWorldFocus2D() = 0;
	virtual void GetImage2D(std::vector<real_t> &image_value, std::vector<Vector2D> &image_gradient,
							size_t width, size_t height, const Box2D &view, MeshImageType type, size_t mode) = 0;

public:
	inline bool IsInitialized() { return m_initialized; }
	inline bool IsSolved() { return m_solved; }

	inline size_t GetModeCount() { return (size_t) m_modes.cols(); }
	inline const Eigen::MatrixXr& GetModes() { return m_modes; }
	inline real_t GetFrequency() { return m_frequency; }

	inline const Eigen::MatrixXr& GetInductanceMatrix() { return m_inductance_matrix; }
	inline const Eigen::MatrixXr& GetCapacitanceMatrix() { return m_capacitance_matrix; }
	inline const Eigen::MatrixXr& GetResistanceMatrix() { return m_resistance_matrix; }
	inline const Eigen::MatrixXr& GetConductanceMatrix() { return m_conductance_matrix; }

	inline const Eigen::MatrixXc& GetCharacteristicImpedanceMatrix() { return m_characteristic_impedance_matrix; }
	inline const Eigen::VectorXc& GetCharacteristicImpedances() { return m_characteristic_impedances; }
	inline const Eigen::VectorXc& GetPropagationConstants() { return m_propagation_constants; }

	inline const Eigen::MatrixXc& GetEigenmodes() { return m_eigenmodes; }
	inline const Eigen::VectorXc& GetEigenmodePropagationConstants() { return m_eigenmode_propagation_constants; }

protected:
	virtual void DoInitialize() = 0;
	virtual void DoSolve() = 0;
	virtual void DoCleanup() = 0;
	virtual size_t GetPortCount() = 0;

};
