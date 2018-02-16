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

#include <sstream>
#include <vector>

inline void ExtendString(std::ostringstream&) {}

template<typename V, typename... Args>
inline void ExtendString(std::ostringstream &ss, V &&value, Args&&... args) {
	ss << std::forward<V>(value);
	ExtendString(ss, std::forward<Args>(args)...);
}

template<typename... Args>
inline std::string MakeString(Args&&... args) {
	std::ostringstream ss;
	ExtendString(ss, std::forward<Args>(args)...);
	return ss.str();
}

template<typename T>
inline std::ostream& operator<<(std::ostream &stream, const std::vector<T> &data) {
	stream << '{';
	if(!data.empty()) {
		stream << data[0];
		for(size_t i = 1; i < data.size(); ++i) {
			stream << ',' << ' ' << data[i];
		}
	}
	stream << '}';
	return stream;
}
