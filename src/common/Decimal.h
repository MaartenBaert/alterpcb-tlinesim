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

#include <cstddef>
#include <utility>

enum FloatType {
	FLOATTYPE_NORMAL,
	FLOATTYPE_INF,
	FLOATTYPE_NAN,
};

struct Decimal {
	FloatType type;
	bool negative;
	int32_t expo;
	uint64_t mant;
};

constexpr size_t DECIMAL_BUFFER_SIZE = 25;
constexpr uint64_t DECIMAL_MANT_MAX = UINT64_C(13446744073709551615);

Decimal ToDecimal(double value, uint32_t precision);
double FromDecimal(Decimal dec);
