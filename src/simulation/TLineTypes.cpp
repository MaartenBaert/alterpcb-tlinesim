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

const char *const TLINERESULT_NAMES[TLINERESULT_COUNT] = {
	"Impedance",
	"Velocity",
	"Wavelength",
	"Loss",
	"Capacitance",
	"Inductance",
	"Alpha",
	"Beta",
};

const char *const TLINERESULT_UNITS[TLINERESULT_COUNT] = {
	"Ohm",
	"m/s",
	"mm",
	"dB/mm",
	"pF/mm",
	"nH/mm",
	"Np/m",
	"rad/m",
};

std::vector<TLineType> g_tline_types;

void RegisterTLine_Microstrip();
void RegisterTLine_CoplanarWaveguide();

void RegisterTLineTypes() {
	RegisterTLine_Microstrip();
	RegisterTLine_CoplanarWaveguide();
}

const MaterialConductor* FindConductor(VDataPathDict &root, const char *key, MaterialDatabase *material_database) {
	std::string name = root.GetMember(key).AsString();
	const MaterialConductor *res = material_database->FindConductor(name);
	if(res == NULL)
		throw std::runtime_error(MakeString("Conductor '", name, "' not found in material database."));
	return res;
}

const MaterialDielectric* FindDielectric(VDataPathDict &root, const char *key, MaterialDatabase *material_database) {
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

void TLineSolveModes(TLineContext &context, const std::vector<real_t> &modes, const std::vector<real_t> &mode_scale) {
	size_t ports = modes.size() / mode_scale.size();

	// initialize
	context.m_output_mesh->Initialize();

	context.m_results.clear();
	context.m_results.resize(TLINERESULT_COUNT * mode_scale.size() * context.m_frequencies.size());
	for(size_t i = 0; i < context.m_frequencies.size(); ++i) {

		// solve
		std::vector<real_t> charges, currents;
		real_t freq = context.m_frequencies[i];
		context.m_output_mesh->Solve(charges, currents, modes, mode_scale.size(), freq);
		assert(charges.size() == modes.size());
		assert(currents.size() == modes.size());

		for(size_t j = 0; j < mode_scale.size(); ++j) {

			// process results
			real_t total_charge = 0.0, total_current = 0.0;
			for(size_t k = 0; k < ports; ++k) {
				total_charge += charges[ports * j + k] * modes[ports * j + k] / sqr(mode_scale[j]);
				total_current += currents[ports * j + k] * modes[ports * j + k] / sqr(mode_scale[j]);
			}
			real_t cap = VACUUM_PERMITTIVITY * total_charge;
			real_t ind = VACUUM_PERMEABILITY / total_current;

			// generate outputs
			real_t *output_values = context.m_results.data() + TLINERESULT_COUNT * (mode_scale.size() * i + j);
			output_values[TLINERESULT_IMPEDANCE] = sqrt(ind / cap);
			output_values[TLINERESULT_VELOCITY] = 1.0 / sqrt(ind * cap);
			output_values[TLINERESULT_WAVELENGTH] = 1.0 / (sqrt(ind * cap) * freq) * 1e3;
			output_values[TLINERESULT_LOSS] = 0.0;
			output_values[TLINERESULT_CAPACITANCE] = cap * 1e9;
			output_values[TLINERESULT_INDUCTANCE] = ind * 1e6;
			output_values[TLINERESULT_ALPHA] = 0.0;
			output_values[TLINERESULT_BETA] = 2.0 * M_PI * freq * sqrt(ind * cap);

		}

		// update progress
		if(context.m_progress_callback) {
			context.m_progress_callback(i + 1);
		}

	}

	// cleanup
	context.m_output_mesh->Cleanup();

}
