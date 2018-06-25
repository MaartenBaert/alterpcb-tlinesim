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
#include "Cow.h"
#include "HashTable.h"

#include <ostream>

// The internal representation of floating point numbers is scaled up by a factor of 1000000. This allows us to store
// common values such as 0.1 and 0.254 exactly, which avoids a lot of problems with rounding errors. The JSON parser is
// aware of this shift and is able to do the necessary string/float conversions without loss of accuracy.
constexpr int32_t VDATA_DECIMAL_SHIFT = 6;
constexpr real_t VDATA_DECIMAL_SHIFT_SCALE = 1000000.0;
constexpr real_t VDATA_DECIMAL_SHIFT_UNSCALE = 0.000001;

constexpr real_t FloatScale(real_t x) { return x * VDATA_DECIMAL_SHIFT_SCALE; }
constexpr real_t FloatUnscale(real_t x) { return x * VDATA_DECIMAL_SHIFT_UNSCALE; }

// These are the corresponding Python typenames. They don't match with the C++ types.
enum VDataType {
	VDATA_NULL,
	VDATA_BOOL,
	VDATA_INT,
	VDATA_FLOAT,
	VDATA_STRING,
	VDATA_LIST,
	VDATA_DICT,
	//VDATA_ARRAY_INT,
	//VDATA_ARRAY_FLOAT,
};

class VDataDictEntry;

struct VDataDictHasher {
	static bool Equal(const VDataDictEntry &a, const VDataDictEntry &b);
	static bool Equal(const VDataDictEntry &a, stringtag_t b);
	static hash_t Hash(hash_t hash, const VDataDictEntry &value);
	static hash_t Hash(hash_t hash, stringtag_t value);
};

// C++11 allows us to put non-POD types in unions, but we have to call the constructors and destructors manually.
// This is equivalent to using a char array and placement new for manual memory management, but cleaner.
class VData {

public:
	typedef std::vector<VData> List;
	typedef HashTable<VDataDictEntry, VDataDictHasher> Dict;

private:
	VDataType m_type;
	union {
		bool m_value_bool;
		int64_t m_value_int;
		real_t m_value_float;
		std::string m_value_string;
		Cow<List> m_value_list;
		Cow<Dict> m_value_dict;
	};

public:
	inline VData() : m_type(VDATA_NULL) {}
	inline VData(std::nullptr_t) : m_type(VDATA_NULL) {}
	inline VData(bool value) : m_type(VDATA_BOOL), m_value_bool(value) {}
	inline VData(int8_t value) : m_type(VDATA_INT), m_value_int(value) {}
	inline VData(uint8_t value) : m_type(VDATA_INT), m_value_int(value) {}
	inline VData(int16_t value) : m_type(VDATA_INT), m_value_int(value) {}
	inline VData(uint16_t value) : m_type(VDATA_INT), m_value_int(value) {}
	inline VData(int32_t value) : m_type(VDATA_INT), m_value_int(value) {}
	inline VData(uint32_t value) : m_type(VDATA_INT), m_value_int(value) {}
	inline VData(int64_t value) : m_type(VDATA_INT), m_value_int(value) {}
	inline VData(uint64_t value) : m_type(VDATA_INT), m_value_int((int64_t) value) {} // overflow possible
	inline VData(float value) : m_type(VDATA_FLOAT), m_value_float(value) {}
	inline VData(double value) : m_type(VDATA_FLOAT), m_value_float(value) {}
	inline VData(const std::string &value) : m_type(VDATA_STRING), m_value_string(value) {}
	inline VData(std::string &&value) : m_type(VDATA_STRING), m_value_string(std::move(value)) {}
	inline VData(const char *value) : m_type(VDATA_STRING), m_value_string(value) { assert(value != NULL); }
	inline VData(const Cow<List> &value) : m_type(VDATA_LIST), m_value_list(value) {}
	inline VData(Cow<List> &&value) : m_type(VDATA_LIST), m_value_list(std::move(value)) {}
	inline VData(const List &value) : m_type(VDATA_LIST), m_value_list(std::make_shared<List>(value)) {}
	inline VData(List &&value) : m_type(VDATA_LIST), m_value_list(std::make_shared<List>(std::move(value))) {}
	inline VData(const Cow<Dict> &value) : m_type(VDATA_DICT), m_value_dict(value) {}
	inline VData(Cow<Dict> &&value) : m_type(VDATA_DICT), m_value_dict(std::move(value)) {}
	inline VData(const Dict &value) : m_type(VDATA_DICT), m_value_dict(std::make_shared<Dict>(value)) {}
	inline VData(Dict &&value) : m_type(VDATA_DICT), m_value_dict(std::make_shared<Dict>(std::move(value))) {}

	inline VData(VDataType type) {
		switch(type) {
			case VDATA_NULL: ConstructNull(); break;
			case VDATA_BOOL: ConstructBool(false); break;
			case VDATA_INT: ConstructInt(0); break;
			case VDATA_FLOAT: ConstructFloat(0.0); break;
			case VDATA_STRING: ConstructString(); break;
			case VDATA_LIST: ConstructList(); break;
			case VDATA_DICT: ConstructDict(); break;
		}
	}

	inline VData(const VData &other) {
		switch(other.m_type) {
			case VDATA_NULL: ConstructNull(); break;
			case VDATA_BOOL: ConstructBool(other.m_value_bool); break;
			case VDATA_INT: ConstructInt(other.m_value_int); break;
			case VDATA_FLOAT: ConstructFloat(other.m_value_float); break;
			case VDATA_STRING: ConstructString(other.m_value_string); break;
			case VDATA_LIST: ConstructList(other.m_value_list); break;
			case VDATA_DICT: ConstructDict(other.m_value_dict); break;
		}
	}

	inline VData(VData &&other) {
		switch(other.m_type) {
			case VDATA_NULL: ConstructNull(); break;
			case VDATA_BOOL: ConstructBool(other.m_value_bool); break;
			case VDATA_INT: ConstructInt(other.m_value_int); break;
			case VDATA_FLOAT: ConstructFloat(other.m_value_float); break;
			case VDATA_STRING: ConstructString(std::move(other.m_value_string)); break;
			case VDATA_LIST: ConstructList(std::move(other.m_value_list)); other.Destruct(); break;
			case VDATA_DICT: ConstructDict(std::move(other.m_value_dict)); other.Destruct(); break;
		}
	}

	inline VData& operator=(const VData &other) {
		if(this != &other) {
			Destruct();
			switch(other.m_type) {
				case VDATA_NULL: ConstructNull(); break;
				case VDATA_BOOL: ConstructBool(other.m_value_bool); break;
				case VDATA_INT: ConstructInt(other.m_value_int); break;
				case VDATA_FLOAT: ConstructFloat(other.m_value_float); break;
				case VDATA_STRING: ConstructString(other.m_value_string); break;
				case VDATA_LIST: ConstructList(other.m_value_list); break;
				case VDATA_DICT: ConstructDict(other.m_value_dict); break;
			}
		}
		return *this;
	}

	inline VData& operator=(VData &&other) {
		if(this != &other) {
			Destruct();
			switch(other.m_type) {
				case VDATA_NULL: ConstructNull(); break;
				case VDATA_BOOL: ConstructBool(other.m_value_bool); break;
				case VDATA_INT: ConstructInt(other.m_value_int); break;
				case VDATA_FLOAT: ConstructFloat(other.m_value_float); break;
				case VDATA_STRING: ConstructString(std::move(other.m_value_string)); break;
				case VDATA_LIST: ConstructList(std::move(other.m_value_list)); other.Destruct(); break;
				case VDATA_DICT: ConstructDict(std::move(other.m_value_dict)); other.Destruct(); break;
			}
		}
		return *this;
	}

	inline ~VData() {
		Destruct();
	}

	inline VData& operator=(std::nullptr_t) { Destruct(); ConstructNull(); return *this; }
	inline VData& operator=(bool value) { Destruct(); ConstructBool(value); return *this; }
	inline VData& operator=(int8_t value) { Destruct(); ConstructInt(value); return *this; }
	inline VData& operator=(uint8_t value) { Destruct(); ConstructInt(value); return *this; }
	inline VData& operator=(int16_t value) { Destruct(); ConstructInt(value); return *this; }
	inline VData& operator=(uint16_t value) { Destruct(); ConstructInt(value); return *this; }
	inline VData& operator=(int32_t value) { Destruct(); ConstructInt(value); return *this; }
	inline VData& operator=(uint32_t value) { Destruct(); ConstructInt(value); return *this; }
	inline VData& operator=(int64_t value) { Destruct(); ConstructInt(value); return *this; }
	inline VData& operator=(uint64_t value) { Destruct(); ConstructInt((int64_t) value); return *this; } // overflow possible
	inline VData& operator=(float value) { Destruct(); ConstructFloat(value); return *this; }
	inline VData& operator=(double value) { Destruct(); ConstructFloat(value); return *this; }
	inline VData& operator=(const std::string &value) { Destruct(); ConstructString(value); return *this; }
	inline VData& operator=(std::string &&value) { Destruct(); ConstructString(value); return *this; }
	inline VData& operator=(const char *value) { Destruct(); ConstructString(value); return *this; }
	inline VData& operator=(const Cow<List> &value) { Destruct(); ConstructList(value); return *this; }
	inline VData& operator=(Cow<List> &&value) { Destruct(); ConstructList(std::move(value)); return *this; }
	inline VData& operator=(const List &value) { Destruct(); ConstructList(value); return *this; }
	inline VData& operator=(List &&value) { Destruct(); ConstructList(std::move(value)); return *this; }
	inline VData& operator=(const Cow<Dict> &value) { Destruct(); ConstructDict(value); return *this; }
	inline VData& operator=(Cow<Dict> &&value) { Destruct(); ConstructDict(std::move(value)); return *this; }
	inline VData& operator=(const Dict &value) { Destruct(); ConstructDict(value); return *this; }
	inline VData& operator=(Dict &&value) { Destruct(); ConstructDict(std::move(value)); return *this; }

	template<typename... Args>
	inline std::string& NewString(Args&&... args) {
		Destruct();
		ConstructString(std::forward<Args>(args)...);
		return m_value_string;
	}

	template<typename... Args>
	inline List& NewList(Args&&... args) {
		Destruct();
		ConstructList(std::forward<Args>(args)...);
		return m_value_list.UniqueBypass();
	}

	template<typename... Args>
	inline Dict& NewDict(Args&&... args) {
		Destruct();
		ConstructDict(std::forward<Args>(args)...);
		return m_value_dict.UniqueBypass();
	}

	inline VDataType GetType() const { return m_type; }
	inline       bool& AsBool()       { assert(m_type == VDATA_BOOL); return m_value_bool; }
	inline const bool& AsBool() const { assert(m_type == VDATA_BOOL); return m_value_bool; }
	inline       int64_t& AsInt()       { assert(m_type == VDATA_INT); return m_value_int; }
	inline const int64_t& AsInt() const { assert(m_type == VDATA_INT); return m_value_int; }
	inline       real_t& AsFloat()       { assert(m_type == VDATA_FLOAT); return m_value_float; }
	inline const real_t& AsFloat() const { assert(m_type == VDATA_FLOAT); return m_value_float; }
	inline       std::string& AsString()       { assert(m_type == VDATA_STRING); return m_value_string; }
	inline const std::string& AsString() const { assert(m_type == VDATA_STRING); return m_value_string; }
	inline       List& AsListUnique() { assert(m_type == VDATA_LIST); return m_value_list.Unique(); }
	inline const List& AsList() const { assert(m_type == VDATA_LIST); return m_value_list.Ref(); }
	inline       Dict& AsDictUnique() { assert(m_type == VDATA_DICT); return m_value_dict.Unique(); }
	inline const Dict& AsDict() const { assert(m_type == VDATA_DICT); return m_value_dict.Ref(); }

private:
	inline void ConstructNull() {
		m_type = VDATA_NULL;
	}

	inline void ConstructBool(bool value) {
		m_value_bool = value;
		m_type = VDATA_BOOL;
	}

	inline void ConstructInt(int64_t value) {
		m_value_int = value;
		m_type = VDATA_INT;
	}

	inline void ConstructFloat(real_t value) {
		m_value_float = value;
		m_type = VDATA_FLOAT;
	}

	template<typename... Args>
	inline void ConstructString(Args&&... args) {
		new(&m_value_string) std::string(std::forward<Args>(args)...);
		m_type = VDATA_STRING;
	}

	inline void ConstructList(const Cow<List> &other) {
		new(&m_value_list) Cow<List>(other);
		m_type = VDATA_LIST;
	}

	inline void ConstructList(Cow<List> &&other) {
		new(&m_value_list) Cow<List>(std::move(other));
		m_type = VDATA_LIST;
	}

	template<typename... Args>
	inline void ConstructList(Args&&... args) {
		new(&m_value_list) Cow<List>(std::make_shared<List>(std::forward<Args>(args)...));
		m_type = VDATA_LIST;
	}

	inline void ConstructDict(const Cow<Dict> &other) {
		new(&m_value_list) Cow<Dict>(other);
		m_type = VDATA_DICT;
	}

	inline void ConstructDict(Cow<Dict> &&other) {
		new(&m_value_list) Cow<Dict>(std::move(other));
		m_type = VDATA_DICT;
	}

	template<typename... Args>
	inline void ConstructDict(Args&&... args) {
		new(&m_value_list) Cow<Dict>(std::make_shared<Dict>(std::forward<Args>(args)...));
		m_type = VDATA_DICT;
	}

	inline void Destruct() {
		switch(m_type) {
			case VDATA_NULL: break;
			case VDATA_BOOL: break;
			case VDATA_INT: break;
			case VDATA_FLOAT: break;
			case VDATA_STRING: m_value_string.~basic_string(); break;
			case VDATA_LIST: m_value_list.~Cow(); break;
			case VDATA_DICT: m_value_dict.~Cow(); break;
		}
		m_type = VDATA_NULL;
	}

};

class VDataDictEntry {

private:
	stringtag_t m_key;
	VData m_value;

public:
	inline explicit VDataDictEntry(stringtag_t key) : m_key(key) {}
	//inline VDataDictEntry(stringtag_t key, const VData& value) : m_key(key), m_value(value) {}
	//inline VDataDictEntry(stringtag_t key, VData&& value) : m_key(key), m_value(std::move(value)) {}
	template<typename... Args>
	inline explicit VDataDictEntry(stringtag_t key, Args&&... args) : m_key(key), m_value(std::forward<Args>(args)...) {}

	// default copy and assignment
	VDataDictEntry(const VDataDictEntry&) = default;
	VDataDictEntry(VDataDictEntry&&) = default;
	VDataDictEntry& operator=(const VDataDictEntry&) = default;
	VDataDictEntry& operator=(VDataDictEntry&&) = default;

	inline stringtag_t Key() const { return m_key; }
	inline VData& Value() { return m_value; }
	inline const VData& Value() const { return m_value; }

};

inline bool VDataDictHasher::Equal(const VDataDictEntry &a, const VDataDictEntry &b) {
	return a.Key() == b.Key();
}

inline bool VDataDictHasher::Equal(const VDataDictEntry &a, stringtag_t b) {
	return a.Key() == b;
}

inline hash_t VDataDictHasher::Hash(hash_t hash, const VDataDictEntry &value) {
	return MurmurHash::HashData(hash, value.Key());
}

inline hash_t VDataDictHasher::Hash(hash_t hash, stringtag_t value) {
	return MurmurHash::HashData(hash, value);
}

std::ostream& operator<<(std::ostream &stream, const VData &data);

int VDataCompare(const VData &a, const VData &b);
bool operator==(const VData &a, const VData &b);

inline bool operator!=(const VData &a, const VData &b) { return !(a == b); }
inline bool operator<(const VData &a, const VData &b) { return (VDataCompare(a, b) < 0); }
inline bool operator>(const VData &a, const VData &b) { return (VDataCompare(a, b) > 0); }
inline bool operator<=(const VData &a, const VData &b) { return (VDataCompare(a, b) <= 0); }
inline bool operator>=(const VData &a, const VData &b) { return (VDataCompare(a, b) >= 0); }

inline void ExtendVList(VData::List&) {}

template<typename V, typename... Args>
inline void ExtendVList(VData::List &ref, V &&value, Args&&... args) {
	ref.emplace_back(std::forward<V>(value));
	ExtendVList(ref, std::forward<Args>(args)...);
}

template<typename... Args>
inline VData MakeVList(Args&&... args) {
	VData data;
	VData::List &ref = data.NewList();
	ref.reserve(sizeof...(Args));
	ExtendVList(ref, std::forward<Args>(args)...);
	return data;
}

inline void ExtendVDict(VData::Dict&) {}

template<typename K, typename V, typename... Args>
inline void ExtendVDict(VData::Dict &ref, K &&key, V &&value, Args&&... args) {
	ref.EmplaceBack(std::forward<K>(key), std::forward<V>(value));
	ExtendVDict(ref, std::forward<Args>(args)...);
}

template<typename... Args>
inline VData MakeVDict(Args&&... args) {
	VData data;
	VData::Dict &ref = data.NewDict(sizeof...(Args));
	ExtendVDict(ref, std::forward<Args>(args)...);
	return data;
}
