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

#include <algorithm>

template<typename E>
class EnumTranslator {

public:
	static const std::string NAME;

private:
	static const EnumTranslator SINGLETON;

private:
	std::vector<std::pair<E, std::string>> m_enum_to_string, m_string_to_enum;

private:
	inline static bool CompareEnum1(const std::pair<E, std::string> &a, E b) {
		return (a.first < b);
	}
	inline static bool CompareEnum2(const std::pair<E, std::string> &a, const std::pair<E, std::string> &b) {
		return (a.first < b.first);
	}
	inline static bool CompareString1(const std::pair<E, std::string> &a, const std::string &b) {
		return (a.second < b);
	}
	inline static bool CompareString2(const std::pair<E, std::string> &a, const std::pair<E, std::string> &b) {
		return (a.second < b.second);
	}

public:
	EnumTranslator(std::initializer_list<std::pair<E, std::string>> list) {
		m_enum_to_string = list;
		m_string_to_enum = list;
		std::stable_sort(m_enum_to_string.begin(), m_enum_to_string.end(), CompareEnum2);
		std::stable_sort(m_string_to_enum.begin(), m_string_to_enum.end(), CompareString2);
	}
	inline static bool EnumToString(const E &e, std::string &s) {
		auto it = std::lower_bound(SINGLETON.m_enum_to_string.begin(), SINGLETON.m_enum_to_string.end(), e, CompareEnum1);
		if(it == SINGLETON.m_enum_to_string.end() || it->first != e)
			return false;
		s = it->second;
		return true;
	}
	inline static bool StringToEnum(const std::string &s, E &e) {
		auto it = std::lower_bound(SINGLETON.m_string_to_enum.begin(), SINGLETON.m_string_to_enum.end(), s, CompareString1);
		if(it == SINGLETON.m_string_to_enum.end() || it->second != s)
			return false;
		e = it->first;
		return true;
	}

};

// convenience functions without full error checking
template<typename E>
inline std::string EnumToString(E e) {
	std::string s;
	bool success = EnumTranslator<E>::EnumToString(e, s);
	assert(success); UNUSED(success);
	return s;
}
template<typename E>
inline E StringToEnum(const std::string &string, E fallback) {
	E e = fallback;
	EnumTranslator<E>::StringToEnum(e, string);
	return e;
}

#define ENUMSTRINGS(E) \
	template<> const std::string EnumTranslator<E>::NAME = #E; \
	template<> const EnumTranslator<E> EnumTranslator<E>::SINGLETON
