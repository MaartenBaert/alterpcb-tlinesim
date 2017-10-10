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
#include "EnumTranslator.h"
#include "StringHelper.h"
#include "StringRegistry.h"
#include "VData.h"

class VDataReader {

private:
	const VData &m_data;
	VDataReader *m_parent;
	stringtag_t m_key;
	size_t m_index;

public:
	inline VDataReader(const VData &data, VDataReader *parent = NULL, stringtag_t key = INDEX_NONE, size_t index = INDEX_NONE)
		: m_data(data), m_parent(parent), m_key(key), m_index(index) {}

	inline bool AsBool() {
		if(m_data.GetType() == VDATA_BOOL)
			return m_data.AsBool();
		throw std::runtime_error(MakeString("Expected '", *this, "' to be bool, got ", EnumToString(m_data.GetType()), " instead."));
	}

	inline int64_t AsInt() {
		if(m_data.GetType() == VDATA_INT)
			return m_data.AsInt();
		throw std::runtime_error(MakeString("Expected '", *this, "' to be int, got ", EnumToString(m_data.GetType()), " instead."));
	}

	inline real_t AsFloat() {
		if(m_data.GetType() == VDATA_INT)
			return (real_t) m_data.AsInt(); // implicit int to float conversion
		if(m_data.GetType() == VDATA_FLOAT)
			return FloatUnscale(m_data.AsFloat());
		throw std::runtime_error(MakeString("Expected '", *this, "' to be float, got ", EnumToString(m_data.GetType()), " instead."));
	}

	inline std::string AsString() {
		if(m_data.GetType() == VDATA_STRING)
			return m_data.AsString();
		throw std::runtime_error(MakeString("Expected '", *this, "' to be string, got ", EnumToString(m_data.GetType()), " instead."));
	}

	inline size_t GetElementCount() {
		if(m_data.GetType() != VDATA_LIST)
			throw std::runtime_error(MakeString("Expected '", *this, "' to be list, got ", EnumToString(m_data.GetType()), " instead."));
		const VData::List &list = m_data.AsList();
		return list.size();
	}

	inline VDataReader GetElement(size_t index) {
		if(m_data.GetType() != VDATA_LIST)
			throw std::runtime_error(MakeString("Expected '", *this, "' to be list, got ", EnumToString(m_data.GetType()), " instead."));
		const VData::List &list = m_data.AsList();
		if(index >= list.size())
			throw std::runtime_error(MakeString("Index '", index, "' not found in '", *this, "'."));
		return VDataReader(list[index], this, INDEX_NONE, index);
	}

	inline VDataReader GetMember(stringtag_t key) {
		if(m_data.GetType() != VDATA_DICT)
			throw std::runtime_error(MakeString("Expected '", *this, "' to be dict, got ", EnumToString(m_data.GetType()), " instead."));
		const VData::Dict &dict = m_data.AsDict();
		size_t index = dict.Find(key);
		if(index == INDEX_NONE)
			throw std::runtime_error(MakeString("Key '", StringRegistry::GetString(key), "' not found in '", *this, "'."));
		return VDataReader(dict[index].Value(), this, key, INDEX_NONE);
	}
	inline VDataReader GetMember(const char *key) {
		return GetMember(StringRegistry::NewTag(key));
	}

	inline VDataReader GetMemberDefault(stringtag_t key, const VData &default_value) {
		if(m_data.GetType() != VDATA_DICT)
			throw std::runtime_error(MakeString("Expected '", *this, "' to be dict, got ", EnumToString(m_data.GetType()), " instead."));
		const VData::Dict &dict = m_data.AsDict();
		size_t index = dict.Find(key);
		if(index == INDEX_NONE)
			return VDataReader(default_value, this, key, INDEX_NONE);
		return VDataReader(dict[index].Value(), this, key, INDEX_NONE);
	}
	inline VDataReader GetMemberDefault(const char *key, const VData &default_value) {
		return GetMemberDefault(StringRegistry::NewTag(key), default_value);
	}

public:
	//inline const VData& GetData() { return m_data; }
	//inline VDataPath* GetParent() { return m_parent; }

public:
	inline friend std::ostream& operator<<(std::ostream &stream, const VDataReader &path) {
		if(path.m_parent != NULL)
			stream << *path.m_parent;
		if(path.m_key != INDEX_NONE) {
			if(path.m_parent != NULL)
				stream << '.';
			stream << StringRegistry::GetString(path.m_key);
		} else if(path.m_index != INDEX_NONE) {
			stream << '[' << path.m_index << ']';
		} else {
			stream << '$';
		}
		return stream;
	}

};

class VDataDictReader {

private:
	const VData::Dict &m_dict;

public:
	inline VDataDictReader(const VData::Dict &dict) : m_dict(dict) {}

	inline VDataReader GetMember(stringtag_t key) {
		size_t index = m_dict.Find(key);
		if(index == INDEX_NONE)
			throw std::runtime_error(MakeString("Key '", StringRegistry::GetString(key), "' not found."));
		return VDataReader(m_dict[index].Value(), NULL, key, INDEX_NONE);
	}
	inline VDataReader GetMember(const char *key) {
		return GetMember(StringRegistry::NewTag(key));
	}

	inline VDataReader GetMemberDefault(stringtag_t key, const VData &default_value) {
		size_t index = m_dict.Find(key);
		if(index == INDEX_NONE)
			return VDataReader(default_value, NULL, key, INDEX_NONE);
		return VDataReader(m_dict[index].Value(), NULL, key, INDEX_NONE);
	}
	inline VDataReader GetMemberDefault(const char *key, const VData &default_value) {
		return GetMemberDefault(StringRegistry::NewTag(key), default_value);
	}

};
