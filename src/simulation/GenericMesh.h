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
#include "Vector.h"

#include <vector>

enum MeshImageType {
	MESHIMAGETYPE_MESH,
	MESHIMAGETYPE_FIELD_E,
	MESHIMAGETYPE_FIELD_H,
};

class GenericMesh {

public:
	GenericMesh();
	virtual ~GenericMesh();

	virtual void Initialize() = 0;
	virtual void Solve(std::vector<real_t> &charges, std::vector<real_t> &currents, const std::vector<real_t> &modes, size_t mode_count, real_t frequency) = 0;
	virtual void Cleanup() = 0;

	virtual Box2D GetWorldBox2D() = 0;
	virtual Box2D GetWorldFocus2D() = 0;
	virtual bool IsInitialized() = 0;
	virtual bool IsSolved() = 0;
	virtual size_t GetModeCount() = 0;
	virtual bool GetImage2D(std::vector<real_t> &image_value, std::vector<Vector2D> &image_gradient, size_t width, size_t height, const Box2D &view, MeshImageType type, size_t mode) = 0;

};
