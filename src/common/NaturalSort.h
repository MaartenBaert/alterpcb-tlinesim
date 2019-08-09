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

namespace NaturalSort {

int Compare(const char *a, size_t a_len, const char *b, size_t b_len);

inline int Compare(const std::string &a, const std::string &b) {
	return Compare(a.data(), a.length(), b.data(), b.length());
}
inline int Compare(const std::string &a, const char *b) {
	return Compare(a.data(), a.length(), b, strlen(b));
}
inline int Compare(const char *a, const std::string &b) {
	return Compare(a, strlen(a), b.data(), b.length());
}

template<typename A, typename B>
inline bool Less(A &&a, B &&b) {
	return Compare(std::forward<A>(a), std::forward<B>(b)) < 0;
}

}
