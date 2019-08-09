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

#include "Materials.h"

#include "EnumTranslator.h"

ENUMSTRINGS(MaterialModel) {
	{MATERIALMODEL_CONSTANT, "constant"},
	{MATERIALMODEL_DJORDJEVIC_SARKAR, "djordjevic-sarkar"},
};

complex_t DjordjevicSarkar(real_t permittivity, real_t loss_tangent, real_t reference_frequency, real_t target_frequency) {
	// The Djordjevic-Sarkar model generates permittivity values that satisfy the Kramers-Kronig relations,
	// which is necessary to ensure that frequency domain simulations will produce causal results.
	// The parameters f_min and f_max exist purely for numerical reasons, they have no physical meaning.
	constexpr real_t f_min = 1.0e3, f_max = 1.0e12;
	//real_t permittivity_slope = permittivity * loss_tangent * 2.0 / M_PI;
	//real_t permittivity_inf = permittivity - permittivity_slope * log(f_max / reference_frequency);

	complex_t reference_log = log(complex_t(f_max, reference_frequency) / complex_t(f_min, reference_frequency));
	complex_t target_log = log(complex_t(f_max, target_frequency) / complex_t(f_min, target_frequency));

	//complex_t(permittivity, -permittivity * loss_tangent) = permittivity_inf + permittivity_slope * reference_log;
	//permittivity = permittivity_inf + permittivity_slope * reference_log.real();
	//-permittivity * loss_tangent = permittivity_inf + permittivity_slope * reference_log.imag();

	complex_t permittivity_slope = permittivity * (loss_tangent + 1.0) / (reference_log.real() - reference_log.imag());
	complex_t permittivity_inf = permittivity - permittivity_slope * reference_log.real();

	//permittivity_inf = permittivity * (reference_log.real() * loss_tangent + reference_log.imag()) / (reference_log.imag() - reference_log.real());
	//permittivity_slope = (permittivity - permittivity_inf) / reference_log.real();

	return  permittivity_inf + permittivity_slope * target_log;
}

void MaterialConductor::GetProperties(real_t frequency, Properties &properties) const {
	properties.conductivity = conductivity;
	properties.permeability = VACUUM_PERMEABILITY * permeability;
	properties.impedance = std::sqrt(complex_t(0.0, 2.0 * M_PI * frequency) * properties.permeability / properties.conductivity);
}

void MaterialDielectric::GetProperties(real_t frequency, MaterialDielectric::Properties &properties) const {
	switch(model) {
		case MATERIALMODEL_CONSTANT: {
			properties.permittivity_x = VACUUM_PERMITTIVITY * permittivity_x * complex_t(1.0, -loss_tangent_x);
			properties.permittivity_y = VACUUM_PERMITTIVITY * permittivity_y * complex_t(1.0, -loss_tangent_y);
			properties.permittivity_z = VACUUM_PERMITTIVITY * permittivity_z * complex_t(1.0, -loss_tangent_z);
			break;
		}
		case MATERIALMODEL_DJORDJEVIC_SARKAR: {
			properties.permittivity_x = VACUUM_PERMITTIVITY * DjordjevicSarkar(permittivity_x, loss_tangent_x, reference_frequency, frequency);
			properties.permittivity_y = VACUUM_PERMITTIVITY * DjordjevicSarkar(permittivity_y, loss_tangent_y, reference_frequency, frequency);
			properties.permittivity_z = VACUUM_PERMITTIVITY * DjordjevicSarkar(permittivity_z, loss_tangent_z, reference_frequency, frequency);
			break;
		}
	}
}
