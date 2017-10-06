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

#include "MaterialDatabase.h"

#include "EnumTranslator.h"
#include "Json.h"
#include "NaturalSort.h"
#include "StringRegistry.h"
#include "VDataPath.h"

template<typename T>
struct NameCompare {
	inline bool operator()(const T &a, const T &b) const {
		return a.m_name < b.m_name;
	}
	inline bool operator()(const T &a, const std::string &b) const {
		return a.m_name < b;
	}
};

template<typename T>
struct NaturalNameCompare {
	inline bool operator()(const T &a, const T &b) const {
		return NaturalStringLess(a.m_name, b.m_name);
	}
	inline bool operator()(const T &a, const std::string &b) const {
		return NaturalStringLess(a.m_name, b);
	}
};

MaterialDatabase::MaterialDatabase() {
	// nothing
}

void MaterialDatabase::LoadFile(const std::string &filename) {

	VData default_permeability = Json::FromString("1.0");
	VData default_permeability_unity_frequency = Json::FromString("1.0e9");

	VData data;
	Json::FromFile(data, filename);
	VDataPath root(data);

	VDataPath conductors = root.GetMember("conductors");
	for(index_t i = 0; i < conductors.GetElementCount(); ++i) {
		VDataPath conductor = conductors.GetElement(i);
		m_conductors.push_back(MaterialConductor{
			conductor.GetMember("name").AsString(),
			FloatUnscale(conductor.GetMember("conductivity").AsFloat()),
			FloatUnscale(conductor.GetMemberDefault("permeability", default_permeability).AsFloat()),
			FloatUnscale(conductor.GetMemberDefault("permeability_unity_frequency", default_permeability_unity_frequency).AsFloat()),
		});
	}

	VDataPath dielectrics = root.GetMember("dielectrics");
	for(index_t i = 0; i < dielectrics.GetElementCount(); ++i) {
		VDataPath dielectric = dielectrics.GetElement(i);
		m_dielectrics.push_back(MaterialDielectric{
			dielectric.GetMember("name").AsString(),
			FloatUnscale(dielectric.GetMember("permittivity_x").AsFloat()),
			FloatUnscale(dielectric.GetMember("permittivity_y").AsFloat()),
			FloatUnscale(dielectric.GetMember("loss_tangent_x").AsFloat()),
			FloatUnscale(dielectric.GetMember("loss_tangent_y").AsFloat()),
			FloatUnscale(dielectric.GetMember("test_frequency").AsFloat()),
		});
	}

}

void MaterialDatabase::Finish() {
	std::stable_sort(m_conductors.begin(), m_conductors.end(), NaturalNameCompare<MaterialConductor>());
	std::stable_sort(m_dielectrics.begin(), m_dielectrics.end(), NaturalNameCompare<MaterialDielectric>());
}

const MaterialConductor *MaterialDatabase::FindConductor(const std::string &name) {
	auto it = std::lower_bound(m_conductors.begin(), m_conductors.end(), name, NameCompare<MaterialConductor>());
	if(it == m_conductors.end() || it->m_name != name)
		return NULL;
	return &(*it);
}

const MaterialDielectric *MaterialDatabase::FindDielectric(const std::string &name) {
	auto it = std::lower_bound(m_dielectrics.begin(), m_dielectrics.end(), name, NameCompare<MaterialDielectric>());
	if(it == m_dielectrics.end() || it->m_name != name)
		return NULL;
	return &(*it);
}

complex_t DjordjevicSarkar(real_t permittivity, real_t loss_tangent, real_t reference_frequency, real_t target_frequency) {
	// The Djordjevic-Sarkar model generates permittivity values that satisfy the Kramers-Kronig relations,
	// which is necessary to ensure that frequency domain simulations will produce causal results.
	// The parameters omega_min and omega_max exist purely for numerical reasons, they have no physical meaning.
	real_t omega_min = 2.0 * M_PI * 1.0e6;
	real_t omega_max = 2.0 * M_PI * 1.0e12;
	real_t reference_omega = 2.0 * M_PI * reference_frequency;
	real_t target_omega = 2.0 * M_PI * target_frequency;
	real_t permittivity_slope = permittivity * loss_tangent * 2.0 / M_PI;
	real_t permittivity_inf = permittivity - permittivity_slope * log(omega_max / reference_omega);
	return  permittivity_inf + permittivity_slope * log(complex_t(omega_max, target_omega) / complex_t(omega_min, target_omega));
}

void GetConductorProperties(MaterialConductorProperties &result, const MaterialConductor *source, real_t target_frequency) {
	result.m_conductivity = source->m_conductivity;
	result.m_permeability = fmin(source->m_permeability, fmax(1.0, source->m_permeability_unity_frequency / target_frequency));
}

void GetDielectricProperties(MaterialDielectricProperties &result, const MaterialDielectric *source, real_t target_frequency) {
	result.m_permittivity_x = DjordjevicSarkar(source->m_permittivity_x, source->m_loss_tangent_x, source->m_test_frequency, target_frequency);
	result.m_permittivity_y = DjordjevicSarkar(source->m_permittivity_y, source->m_loss_tangent_y, source->m_test_frequency, target_frequency);
}
