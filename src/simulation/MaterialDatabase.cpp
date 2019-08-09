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

#include "Json.h"
#include "NaturalSort.h"
#include "StringRegistry.h"
#include "VDataReader.h"

template<typename T>
struct NaturalNameCompare {
	inline bool operator()(const T &a, const T &b) const {
		return NaturalSort::Less(a.name, b.name);
	}
	inline bool operator()(const T &a, const std::string &b) const {
		return NaturalSort::Less(a.name, b);
	}
};

MaterialDatabase::MaterialDatabase() {
	// nothing
}

void MaterialDatabase::LoadFile(const std::string &filename) {

	VData data;
	Json::FromFile(data, filename);
	VDataReader root(data);

	VDataReader conductors = root.GetMember("conductors");
	for(size_t i = 0; i < conductors.GetElementCount(); ++i) {
		VDataReader conductor = conductors.GetElement(i);
		m_conductors.push_back(MaterialConductor{
			conductor.GetMember("name").AsString(),
			conductor.GetMember("conductivity").AsFloat(),
			conductor.GetMember("permeability").AsFloat(),
		});
	}

	VDataReader dielectrics = root.GetMember("dielectrics");
	for(size_t i = 0; i < dielectrics.GetElementCount(); ++i) {
		VDataReader dielectric = dielectrics.GetElement(i);
		m_dielectrics.emplace_back();
		MaterialDielectric &mat = m_dielectrics.back();
		mat.name = dielectric.GetMember("name").AsString();
		mat.model = dielectric.GetMember("model").AsEnum<MaterialModel>();
		mat.permittivity_x = dielectric.GetMember("permittivity_x").AsFloat();
		mat.permittivity_y = dielectric.GetMember("permittivity_y").AsFloat();
		mat.permittivity_z = dielectric.GetMember("permittivity_z").AsFloat();
		mat.loss_tangent_x = dielectric.GetMember("loss_tangent_x").AsFloat();
		mat.loss_tangent_y = dielectric.GetMember("loss_tangent_y").AsFloat();
		mat.loss_tangent_z = dielectric.GetMember("loss_tangent_z").AsFloat();
		switch(mat.model) {
			case MATERIALMODEL_CONSTANT: {
				break;
			}
			case MATERIALMODEL_DJORDJEVIC_SARKAR: {
				mat.reference_frequency = dielectric.GetMember("reference_frequency").AsFloat();
				break;
			}
		}
	}

}

void MaterialDatabase::Finish() {
	std::stable_sort(m_conductors.begin(), m_conductors.end(), NaturalNameCompare<MaterialConductor>());
	std::stable_sort(m_dielectrics.begin(), m_dielectrics.end(), NaturalNameCompare<MaterialDielectric>());
}

const MaterialConductor *MaterialDatabase::FindConductor(const std::string &name) {
	auto it = std::lower_bound(m_conductors.begin(), m_conductors.end(), name, NaturalNameCompare<MaterialConductor>());
	if(it == m_conductors.end() || it->name != name)
		return NULL;
	return &(*it);
}

const MaterialDielectric *MaterialDatabase::FindDielectric(const std::string &name) {
	auto it = std::lower_bound(m_dielectrics.begin(), m_dielectrics.end(), name, NaturalNameCompare<MaterialDielectric>());
	if(it == m_dielectrics.end() || it->name != name)
		return NULL;
	return &(*it);
}
