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

#include "VData.h"

#include "EnumTranslator.h"
#include "StringRegistry.h"

ENUMSTRINGS(VDataType) {
	{VDATA_NULL, "null"},
	{VDATA_BOOL, "bool"},
	{VDATA_INT, "int"},
	{VDATA_FLOAT, "float"},
	{VDATA_STRING, "string"},
	{VDATA_LIST, "list"},
	{VDATA_DICT, "dict"},
};

// This is purely for convenience (e.g. for debugging), it doesn't produce 100% correct results in all cases.
// For example, strings are not escaped, floats may be printed with minor rounding errors, ...
// Use Json::ToStream or Json::ToString when accurate results are required.
std::ostream& operator<<(std::ostream &stream, const VData &data) {
	switch(data.GetType()) {
		case VDATA_NULL: {
			stream << "null";
			break;
		}
		case VDATA_BOOL: {
			stream << ((data.AsBool())? "true" : "false");
			break;
		}
		case VDATA_INT: {
			stream << data.AsInt();
			break;
		}
		case VDATA_FLOAT: {
			stream << FloatUnscale(data.AsFloat());
			break;
		}
		case VDATA_STRING: {
			stream << '"' << data.AsString() << '"';
			break;
		}
		case VDATA_LIST: {
			stream << '[';
			const VData::List &ref = data.AsList();
			if(ref.size() != 0) {
				stream << ref[0];
				for(size_t i = 1; i < ref.size(); ++i) {
					stream << ", " << ref[i];
				}
			}
			stream << ']';
			break;
		}
		case VDATA_DICT: {
			stream << '{';
			const VData::Dict &ref = data.AsDict();
			if(ref.GetSize() != 0) {
				stream << '"' << StringRegistry::GetString(ref[0].Key()) << "\": " << ref[0].Value();
				for(size_t i = 1; i < ref.GetSize(); ++i) {
					stream << ", \"" << StringRegistry::GetString(ref[i].Key()) << "\": " << ref[i].Value();
				}
			}
			stream << '}';
			break;
		}
	}
	return stream;
}

inline int IntCompare(int64_t  a, int64_t  b) {
	return (a == b)? 0 : (a < b)? -1 : 1;
}

inline int FloatCompare(real_t a, real_t b) {
	// this is necessary to produce consistent results for NaN and infinity
	int64_t ai = MemCast<int64_t>(a), bi = MemCast<int64_t>(b);
	ai ^= (ai >> 63) & INT64_MAX;
	bi ^= (bi >> 63) & INT64_MAX;
	return (ai == bi)? 0 : (ai < bi)? -1 : 1;
}

int VDataCompare(const VData &a, const VData &b) {
	if(a.GetType() != b.GetType())
		return (int) a.GetType() - (int) b.GetType();
	switch(a.GetType()) {
		case VDATA_NULL: return 0;
		case VDATA_BOOL: return (int) a.AsBool() - (int) b.AsBool();
		case VDATA_INT: return IntCompare(a.AsInt(), b.AsInt());
		case VDATA_FLOAT: return FloatCompare(a.AsFloat(), b.AsFloat());
		case VDATA_STRING: return a.AsString().compare(b.AsString());
		case VDATA_LIST: {
			const VData::List &ref1 = a.AsList(), &ref2 = b.AsList();
			if(ref1.size() != ref2.size())
				return (ref1.size() < ref2.size())? -1 : 1;
			for(size_t i = 0; i < ref1.size(); ++i) {
				int res = VDataCompare(ref1[i], ref2[i]);
				if(res != 0)
					return res;
			}
			return 0;
		}
		case VDATA_DICT: {
			const VData::Dict &ref1 = a.AsDict(), &ref2 = b.AsDict();
			if(ref1.GetSize() != ref2.GetSize())
				return (ref1.GetSize() < ref2.GetSize())? -1 : 1;
			for(size_t i = 0; i < ref1.GetSize(); ++i) {
				if(ref1[i].Key() != ref2[i].Key())
					return (ref1[i].Key() < ref2[i].Key())? -1 : 1;
				int res = VDataCompare(ref1[i].Value(), ref2[i].Value());
				if(res != 0)
					return res;
			}
			return 0;
		}
	}
	// this should never be reached
	assert(false);
	return 0;
}

bool operator==(const VData &a, const VData &b) {
	if(a.GetType() != b.GetType())
		return false;
	switch(a.GetType()) {
		case VDATA_NULL: return true;
		case VDATA_BOOL: return (a.AsBool() == b.AsBool());
		case VDATA_INT: return (a.AsInt() == b.AsInt());
		case VDATA_FLOAT: return (MemCast<int64_t>(a.AsFloat()) == MemCast<int64_t>(b.AsFloat()));
		case VDATA_STRING: return (a.AsString() == b.AsString());
		case VDATA_LIST: {
			const VData::List &ref1 = a.AsList(), &ref2 = b.AsList();
			if(ref1.size() != ref2.size())
				return false;
			for(size_t i = 0; i < ref1.size(); ++i) {
				if(ref1[i] != ref2[i])
					return false;
			}
			return true;
		}
		case VDATA_DICT: {
			const VData::Dict &ref1 = a.AsDict(), &ref2 = b.AsDict();
			if(ref1.GetSize() != ref2.GetSize())
				return false;
			for(size_t i = 0; i < ref1.GetSize(); ++i) {
				if(ref1[i].Key() != ref2[i].Key())
					return false;
				if(ref1[i].Value() != ref2[i].Value())
					return false;
			}
			return true;
		}
	}
	// this should never be reached
	assert(false);
	return false;
}
