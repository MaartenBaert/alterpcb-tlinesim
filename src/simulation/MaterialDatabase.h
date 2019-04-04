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

constexpr real_t VACUUM_PERMITTIVITY = 8.854187817620389e-12;
constexpr real_t VACUUM_PERMEABILITY = 1.2566370614359173e-6;
constexpr real_t SPEED_OF_LIGHT = 299792458.0;

struct MaterialConductor {
	std::string m_name;
	real_t m_conductivity;
	real_t m_permeability;
	real_t m_permeability_unity_frequency;
};

struct MaterialDielectric {
	std::string m_name;
	real_t m_permittivity_x, m_permittivity_y;
	real_t m_loss_tangent_x, m_loss_tangent_y;
	real_t m_test_frequency;
};

struct MaterialConductorProperties {
	real_t m_conductivity;
	real_t m_permeability;
	real_t m_surface_conductivity;
};

struct MaterialDielectricProperties {
	std::complex<real_t> m_permittivity_x, m_permittivity_y;
};

class MaterialDatabase {

private:
	std::vector<MaterialConductor> m_conductors;
	std::vector<MaterialDielectric> m_dielectrics;

public:
	MaterialDatabase();

	void LoadFile(const std::string &filename);
	void Finish();

	const MaterialConductor* FindConductor(const std::string &name);
	const MaterialDielectric* FindDielectric(const std::string &name);

public:
	inline const std::vector<MaterialConductor>& GetConductors() { return m_conductors; }
	inline const std::vector<MaterialDielectric>& GetDielectrics() { return m_dielectrics; }

};

void GetConductorProperties(MaterialConductorProperties &result, const MaterialConductor *source, real_t target_frequency);
void GetDielectricProperties(MaterialDielectricProperties &result, const MaterialDielectric *source, real_t target_frequency);
