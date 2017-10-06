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
#include "GenericMesh.h"
#include "VData.h"
#include "VDataPath.h"

#include <memory>
#include <string>
#include <vector>

class MaterialDatabase;
class MaterialConductor;
class MaterialDielectric;

enum TLineParameterType {
	TLINE_PARAMETERTYPE_BOOL,
	TLINE_PARAMETERTYPE_REAL,
	TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR,
	TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC,
	TLINE_PARAMETERTYPE_COUNT, // must be last
};

enum TLineResult {
	TLINERESULT_IMPEDANCE,
	TLINERESULT_VELOCITY,
	TLINERESULT_WAVELENGTH,
	TLINERESULT_LOSS,
	TLINERESULT_CAPACITANCE,
	TLINERESULT_INDUCTANCE,
	TLINERESULT_ALPHA,
	TLINERESULT_BETA,
	TLINERESULT_COUNT, // must be last
};

struct TLineParameter {
	std::string m_name;
	TLineParameterType m_type;
	VData m_default_value;
	bool m_unit_mm;
	bool m_separator;
};

struct TLineContext {
	MaterialDatabase *m_material_database;
	std::vector<real_t> m_frequencies;
	VData::Dict m_parameters;
	std::vector<real_t> m_results;
	std::unique_ptr<GenericMesh> m_output_mesh;
};

typedef void (*TLineSimulate)(TLineContext&);

struct TLineType {
	std::string m_name, m_description;
	std::vector<TLineParameter> m_parameters;
	std::vector<std::string> m_modes;
	TLineSimulate m_simulate;
};

extern const char *const TLINERESULT_NAMES[TLINERESULT_COUNT];
extern const char *const TLINERESULT_UNITS[TLINERESULT_COUNT];

extern std::vector<TLineType> g_tline_types;

void RegisterTLineTypes();

std::string CanonicalName(const std::string& name);
const MaterialConductor* FindConductor(VDataPathDict &root, const char *key, MaterialDatabase *material_database);
const MaterialDielectric* FindDielectric(VDataPathDict &root, const char *key, MaterialDatabase *material_database);
void TLineSolveModes(TLineContext &context, const std::vector<real_t> &modes, const std::vector<real_t> &mode_scale);
