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

int NaturalStringCompare(const char *a, size_t a_len, const char *b, size_t b_len);

inline int NaturalStringCompare(const std::string &a, const std::string &b) {
	return NaturalStringCompare(a.data(), a.length(), b.data(), b.length());
}
inline int NaturalStringCompare(const std::string &a, const char *b) {
	return NaturalStringCompare(a.data(), a.length(), b, strlen(b));
}
inline int NaturalStringCompare(const char *a, const std::string &b) {
	return NaturalStringCompare(a, strlen(a), b.data(), b.length());
}

template<typename A, typename B>
inline bool NaturalStringLess(A &&a, B &&b) {
	return NaturalStringCompare(std::forward<A>(a), std::forward<B>(b)) < 0;
}
