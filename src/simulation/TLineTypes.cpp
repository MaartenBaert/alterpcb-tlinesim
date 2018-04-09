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

#include "TLineTypes.h"

#include "MaterialDatabase.h"

// TODO: remove
#include <iostream>

const char *const TLINERESULT_NAMES[TLINERESULT_COUNT] = {
	"Impedance",
	"Velocity",
	"Wavelength",
	"Loss",
	"Inductance",
	"Capacitance",
	"Resistance",
	"Conductance",
	"Alpha",
	"Beta",
};

const char *const TLINERESULT_UNITS[TLINERESULT_COUNT] = {
	"Ω",
	"m/s",
	"mm",
	"dB/mm",
	"nH/mm",
	"pF/mm",
	"Ω/mm",
	"S/mm",
	"Np/m",
	"rad/m",
};

std::vector<TLineType> g_tline_types;

void RegisterTLine_Microstrip();
void RegisterTLine_Stripline();
void RegisterTLine_CoplanarWaveguide();

void RegisterTLineTypes() {
	RegisterTLine_Microstrip();
	RegisterTLine_Stripline();
	RegisterTLine_CoplanarWaveguide();
}

const MaterialConductor* FindConductor(VDataDictReader &root, const char *key, MaterialDatabase *material_database) {
	std::string name = root.GetMember(key).AsString();
	const MaterialConductor *res = material_database->FindConductor(name);
	if(res == NULL)
		throw std::runtime_error(MakeString("Conductor '", name, "' not found in material database."));
	return res;
}

const MaterialDielectric* FindDielectric(VDataDictReader &root, const char *key, MaterialDatabase *material_database) {
	std::string name = root.GetMember(key).AsString();
	const MaterialDielectric *res = material_database->FindDielectric(name);
	if(res == NULL)
		throw std::runtime_error(MakeString("Dielectric '", name, "' not found in material database."));
	return res;
}

std::string CanonicalName(const std::string &name) {
	std::string out;
	for(size_t i = 0; i < name.size(); ++i) {
		char c = name[i];
		if(c >= 'a' && c <= 'z') {
			out += c;
		} else if(c >= 'A' && c <= 'Z') {
			out += (char) (c - 'A' + 'a');
		} else if(c >= '0' && c <= '9') {
			out += c;
		} else {
			out += '_';
		}
	}
	return out;
}

void TLineSolveModes(TLineContext &context, const Eigen::MatrixXr &modes) {

	// initialize
	context.m_output_mesh->Initialize();

	context.m_results.clear();
	context.m_results.resize(TLINERESULT_COUNT * modes.cols() * context.m_frequencies.size());
	for(size_t i = 0; i < context.m_frequencies.size(); ++i) {

		// solve
		context.m_output_mesh->Solve(modes, context.m_frequencies[i]);

		real_t omega = 2.0 * M_PI * context.m_frequencies[i];
		for(size_t j = 0; j < (size_t) modes.cols(); ++j) {

			// process results
			real_t ind = context.m_output_mesh->GetInductanceMatrix()(j, j);
			real_t cap = context.m_output_mesh->GetCapacitanceMatrix()(j, j);
			real_t res = context.m_output_mesh->GetResistanceMatrix()(j, j);
			real_t cond = context.m_output_mesh->GetConductanceMatrix()(j, j);
			complex_t z0 = context.m_output_mesh->GetCharacteristicImpedances()[j];
			complex_t gamma = context.m_output_mesh->GetPropagationConstants()[j];

			// generate outputs
			real_t *output_values = context.m_results.data() + TLINERESULT_COUNT * (modes.cols() * i + j);
			output_values[TLINERESULT_IMPEDANCE] = z0.real();
			output_values[TLINERESULT_VELOCITY] = omega / gamma.imag();
			output_values[TLINERESULT_WAVELENGTH] = 2.0 * M_PI / gamma.imag() * 1e3;
			output_values[TLINERESULT_LOSS] = gamma.real() * 20.0 / M_LN10 * 1e-3;
			output_values[TLINERESULT_INDUCTANCE] = ind * 1e6;
			output_values[TLINERESULT_CAPACITANCE] = cap * 1e9;
			output_values[TLINERESULT_RESISTANCE] = res * 1e-3;
			output_values[TLINERESULT_CONDUCTANCE] = cond * 1e-3;
			output_values[TLINERESULT_ALPHA] = gamma.real();
			output_values[TLINERESULT_BETA] = gamma.imag();

		}

		// update progress
		if(context.m_progress_callback) {
			context.m_progress_callback(i + 1);
		}

	}

	// cleanup
	context.m_output_mesh->Cleanup();

}
