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

enum MaterialModel {
	MATERIALMODEL_CONSTANT,
	MATERIALMODEL_DJORDJEVIC_SARKAR,
};

struct MaterialConductor {

	struct Properties {
		real_t conductivity;
		real_t permeability;
		complex_t impedance;
	};

	std::string name;
	real_t conductivity;
	real_t permeability;

	void GetProperties(real_t frequency, Properties &properties) const;

};

struct MaterialDielectric {

	struct Properties {
		complex_t permittivity_x, permittivity_y, permittivity_z;
	};

	std::string name;
	MaterialModel model;
	real_t permittivity_x, permittivity_y, permittivity_z;
	real_t loss_tangent_x, loss_tangent_y, loss_tangent_z;
	real_t reference_frequency;

	void GetProperties(real_t frequency, Properties &properties) const;

};
