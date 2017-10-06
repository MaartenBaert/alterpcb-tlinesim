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

#include "StringHelper.h"
#include "VData.h"

#include <streambuf>
#include <string>

namespace Json {

class ParseError : public std::runtime_error {
private:
	size_t m_line, m_column, m_length; // these start at zero - add one when displaying to the user
public:
	template<typename... Args>
	inline ParseError(size_t line, size_t column, size_t length, Args&&... args)
		: std::runtime_error(MakeString("JSON parse error: ", std::forward<Args>(args)...,
										" (line ", line + 1, ", column ", column + 1, ", length ", length, ")")),
		m_line(line), m_column(column), m_length(length) {}
	inline size_t GetLine() { return m_line; }
	inline size_t GetColumn() { return m_column; }
	inline size_t GetLength() { return m_length; }
};

struct Format {
	bool multiline;
	uint32_t precision;
	bool engineering;
	inline Format(bool multiline = true, uint32_t precision = 17, bool engineering = false)
		: multiline(multiline), precision(precision), engineering(engineering) {}
};

extern const Format DEFAULT_SINGLELINE_FORMAT, DEFAULT_MULTILINE_FORMAT;

void FromStream(VData &data, std::streambuf *stream);
void FromString(VData &data, const std::string &str);
void FromFile(VData &data, const std::string &filename);

void ToStream(const VData &data, std::streambuf *stream, const Format &format = DEFAULT_SINGLELINE_FORMAT);
void ToString(const VData &data, std::string &str, const Format &format = DEFAULT_SINGLELINE_FORMAT);
void ToFile(const VData &data, const std::string &filename, const Format &format = DEFAULT_SINGLELINE_FORMAT);

// convenience functions, possibly slower
inline VData FromString(const std::string &str) {
	VData data;
	FromString(data, str);
	return data;
}
inline std::string ToString(const VData &data, const Format &format = DEFAULT_SINGLELINE_FORMAT) {
	std::string str;
	ToString(data, str, format);
	return str;
}

}
