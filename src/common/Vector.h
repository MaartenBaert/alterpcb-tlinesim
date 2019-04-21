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
#include "MiscMath.h"
#include "Vector.h"

struct Vector2D {
	real_t x, y;
	inline Vector2D() {}
	inline Vector2D(real_t x, real_t y) : x(x), y(y) {}
};

struct Box2D {
	real_t x1, x2, y1, y2;
	inline Box2D() {}
	inline Box2D(real_t x1, real_t x2, real_t y1, real_t y2) : x1(x1), x2(x2), y1(y1), y2(y2) {}
	inline Box2D Normalize() const {
		return Box2D(std::min(x1, x2), std::max(x1, x2), std::min(y1, y2), std::max(y1, y2));
	}
	inline Box2D MirrorX(real_t ref = 0.0) const {
		return Box2D(2.0 * ref - x2, 2.0 * ref - x1, y1, y2);
	}
	inline Box2D MirrorY(real_t ref = 0.0) const {
		return Box2D(x1, x2, 2.0 * ref - y2, 2.0 * ref - y1);
	}
	inline Box2D Clip(const Box2D &clip) const {
		return Box2D(clamp(x1, clip.x1, clip.x2), clamp(x2, clip.x1, clip.x2), clamp(y1, clip.y1, clip.y2), clamp(y2, clip.y1, clip.y2));
	}
	inline real_t CenterX() const {
		return 0.5 * (x1 + x2);
	}
	inline real_t CenterY() const {
		return 0.5 * (y1 + y2);
	}
	inline Vector2D Center() const {
		return Vector2D(CenterX(), CenterY());
	}
};

struct GLvec2 {
	float x, y;
	inline GLvec2() {}
	inline GLvec2(float x, float y) : x(x), y(y) {}
	inline GLvec2(const Vector2D &v) : x((float) v.x), y((float) v.y) {}
	inline operator float*() { return &x; }
};

struct GLvec3 {
	float x, y, z;
	inline GLvec3() {}
	inline GLvec3(float x, float y, float z) : x(x), y(y), z(z) {}
	inline operator float*() { return &x; }
};

struct GLvec4 {
	float x, y, z, w;
	inline GLvec4() {}
	inline GLvec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
	inline operator float*() { return &x; }
};

static_assert(sizeof(GLvec2) == sizeof(float) * 2, "Invalid GLvec2 struct!");
static_assert(sizeof(GLvec3) == sizeof(float) * 3, "Invalid GLvec3 struct!");
static_assert(sizeof(GLvec4) == sizeof(float) * 4, "Invalid GLvec4 struct!");
