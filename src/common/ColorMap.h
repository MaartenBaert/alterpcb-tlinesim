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
#include "Color.h"

#include <initializer_list>
#include <vector>

class ColorMap {

private:
	std::vector<Color> m_colors;
	real_t m_scale;

public:
	inline ColorMap(std::initializer_list<Color> list) : m_colors(list), m_scale((real_t) (list.size() - 1)) {}

	inline Color operator()(real_t val) const {
		real_t scaled = clamp<real_t>(val, 0.0, 1.0) * m_scale;
		size_t index = clamp<ptrdiff_t>(rints(scaled - 0.5), 0, m_colors.size() - 2);
		real_t frac = scaled - (real_t) index;
		return ColorMix(m_colors[index], m_colors[index + 1], (float) frac);
	}

};

extern const ColorMap COLORMAP_GRAYSCALE;
extern const ColorMap COLORMAP_MAGMA;
