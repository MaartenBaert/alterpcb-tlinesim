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
#include "Materials.h"

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
