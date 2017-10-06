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

#include "Json.h"

#include "Decimal.h"
#include "StringRegistry.h"

#include <fstream>
#include <limits>
#include <sstream>
#include <streambuf>

namespace Json {

const Format DEFAULT_SINGLELINE_FORMAT(false), DEFAULT_MULTILINE_FORMAT(true);

struct ReadContext {
	std::streambuf *stream;
	size_t line, column;
	template<typename... Args>
	inline void Error(size_t length, Args&&... args) {
		throw ParseError(line, (column > length)? column - length : 0, length, std::forward<Args>(args)...);
	}
};

struct WriteContext {
	std::streambuf *stream;
	Format format;
	uint32_t indent;
};

// ==================== Lookup tables ====================

constexpr char LUT_HEX[16] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};

// ==================== I/O helpers ====================

inline int ReadCharRaw(ReadContext &context) {
	int c = context.stream->sbumpc();
	if(c == '\n') {
		++context.line;
		context.column = 0;
	} else if(c != EOF) {
		++context.column;
	}
	return c;
}

inline int PeekCharRaw(ReadContext &context) {
	return context.stream->sgetc();
}

inline char ReadChar(ReadContext &context) {
	int c = ReadCharRaw(context);
	if(c == EOF)
		context.Error(0, "Unexpected end of input");
	return (char) c;
}

inline char PeekChar(ReadContext &context) {
	int c = PeekCharRaw(context);
	if(c == EOF)
		context.Error(0, "Unexpected end of input");
	return (char) c;
}

inline void SkipWhitespace(ReadContext &context) {
	for( ; ; ) {
		int c = PeekCharRaw(context);
		if(c != ' ' && c != '\t' && c != '\n' && c != '\r')
			break;
		ReadCharRaw(context);
	}
}

inline void WriteChar(WriteContext &context, char c) {
	if(context.stream->sputc(c) == EOF)
		throw std::runtime_error("Could not write to output");
}

inline void WriteData(WriteContext &context, const char *data, std::streamsize len) {
	if(context.stream->sputn(data, len) != len)
		throw std::runtime_error("Could not write to output");
}

inline void WriteLine(WriteContext &context, bool space) {
	if(context.format.multiline) {
		WriteChar(context, '\n');
		for(uint32_t i = 0; i < context.indent; ++i) {
			WriteChar(context, '\t');
		}
	} else if(space) {
		WriteChar(context, ' ');
	}
}

// ==================== Hexadecimal helper ====================

inline uint32_t ReadHex(ReadContext &context) {
	char c = ReadChar(context);
	if(c >= '0' && c <= '9')
		return c - '0';
	if(c >= 'a' && c <= 'f')
		return c - 'a' + 10;
	if(c >= 'A' && c <= 'F')
		return c - 'A' + 10;
	context.Error(1, "Expected hexadecimal value");
	return 0; // this is never reached
}

// ==================== String helpers ====================

inline void ReadString(ReadContext &context, std::string &str) {
	for( ; ; ) {
		char c = ReadChar(context);
		if(c == '"')
			break;
		if(c == '\\') {
			c = ReadChar(context);
			if(c == '"') {
				str += '"';
			} else if(c == '\\') {
				str += '\\';
			} else if(c == '/') {
				str += '/';
			} else if(c == 'b') {
				str += '\b';
			} else if(c == 'f') {
				str += '\f';
			} else if(c == 'n') {
				str += '\n';
			} else if(c == 'r') {
				str += '\r';
			} else if(c == 't') {
				str += '\t';
			} else if(c == 'u') {
				uint32_t codepoint = 0;
				codepoint |= ReadHex(context) << 12;
				codepoint |= ReadHex(context) << 8;
				codepoint |= ReadHex(context) << 4;
				codepoint |= ReadHex(context);
				if(codepoint < 0x0080) {
					str += (char) codepoint;
				} else if(codepoint < 0x0800) {
					str += (char) (0xc0 | (codepoint >> 6));
					str += (char) (0x80 | (codepoint & 0x3f));
				} else if(codepoint >= 0xd800 && codepoint <= 0xdbff) { // high surrogate
					if(ReadChar(context) != '\\')
						context.Error(1, "Expected low surrogate after high surrogate");
					if(ReadChar(context) != 'u')
						context.Error(2, "Expected low surrogate after high surrogate");
					uint32_t codepoint2 = 0;
					codepoint2 |= ReadHex(context) << 12;
					codepoint2 |= ReadHex(context) << 8;
					codepoint2 |= ReadHex(context) << 4;
					codepoint2 |= ReadHex(context);
					if(codepoint2 >= 0xdc00 && codepoint2 <= 0xdfff) { // low surrogate
						codepoint = 0x10000 + (((codepoint - 0xd800) << 10) | (codepoint2 - 0xdc00));
						str += (char) (0xf0 | (codepoint >> 18));
						str += (char) (0x80 | ((codepoint >> 12) & 0x3f));
						str += (char) (0x80 | ((codepoint >> 6) & 0x3f));
						str += (char) (0x80 | (codepoint & 0x3f));
					} else {
						context.Error(6, "Expected low surrogate after high surrogate");
					}
				} else if(codepoint >= 0xdc00 && codepoint <= 0xdfff) { // low surrogate
					context.Error(6, "Unexpected low surrogate");
				} else {
					str += (char) (0xe0 | (codepoint >> 12));
					str += (char) (0x80 | ((codepoint >> 6) & 0x3f));
					str += (char) (0x80 | (codepoint & 0x3f));
				}
			} else {
				context.Error(1, "Invalid escape sequence");
			}
		} else if(c < 32 || c == 127) {
			context.Error(1, "Invalid character with code ", c);
		} else {
			str += c;
		}
	}
}

inline void WriteString(WriteContext &context, const std::string &str) {
	WriteChar(context, '"');
	for(unsigned char c : str) {
		if(c == '"') {
			WriteChar(context, '\\');
			WriteChar(context, '"');
		} else if(c == '\\') {
			WriteChar(context, '\\');
			WriteChar(context, '\\');
		} else if(c == '\b') {
			WriteChar(context, '\\');
			WriteChar(context, '\b');
		} else if(c == '\f') {
			WriteChar(context, '\\');
			WriteChar(context, '\f');
		} else if(c == '\n') {
			WriteChar(context, '\\');
			WriteChar(context, '\n');
		} else if(c == '\r') {
			WriteChar(context, '\\');
			WriteChar(context, '\r');
		} else if(c == '\t') {
			WriteChar(context, '\\');
			WriteChar(context, '\t');
		} else if(c < 32 || c == 127) {
			char code[6] = {'\\', 'u', '0', '0', LUT_HEX[c >> 4], LUT_HEX[c & 0xf]};
			WriteData(context, code, 6);
		} else {
			WriteChar(context, c);
		}
	}
	WriteChar(context, '"');
}

// ==================== Integer helpers ====================

inline void ReadName(ReadContext &context, std::string& name) {
	for( ; ; ) {
		int c = PeekCharRaw(context);
		if(c >= 'a' && c <= 'z') {
			ReadChar(context);
			name += (char) c;
		} else {
			return;
		}
	}
}

inline void ReadNameOrNumber(ReadContext &context, VData &data) {

	// read sign if present
	char prefix = PeekChar(context);
	if(prefix == '+' || prefix == '-') {
		ReadCharRaw(context);
	}

	SkipWhitespace(context);

	// is it a name or a number?
	char c = PeekChar(context);
	if(c >= 'a' && c <= 'z') {

		std::string name;
		ReadName(context, name);

		if(name == "inf") {
			Decimal dec;
			dec.negative = (prefix == '-');
			dec.type = FLOATTYPE_INF;
			data = FromDecimal(dec);
		} else if(name == "nan") {
			Decimal dec;
			dec.negative = (prefix == '-');
			dec.type = FLOATTYPE_NAN;
			data = FromDecimal(dec);
		} else if(name == "true") {
			if(prefix == '+' || prefix == '-')
				context.Error(name.size(), "Boolean values can't have a sign");
			data = true;
		} else if(name == "false") {
			if(prefix == '+' || prefix == '-')
				context.Error(name.size(), "Boolean values can't have a sign");
			data = false;
		} else if(name == "null") {
			if(prefix == '+' || prefix == '-')
				context.Error(name.size(), "Null can't have a sign");
			data = nullptr;
		} else {
			context.Error(name.size(), "Invalid keyword");
		}

	} else if((c >= '0' && c <= '9') || c == '.') {

		Decimal dec;
		dec.negative = (prefix == '-');
		dec.type = FLOATTYPE_NORMAL;
		dec.expo = VDATA_DECIMAL_SHIFT;
		dec.mant = 0;

		bool starts_with_point = (c == '.');
		bool point = false, limit = false;
		for( ; ; ) {
			int c2 = PeekCharRaw(context);
			if(c2 >= '0' && c2 <= '9') {
				ReadCharRaw(context);
				if(limit) {
					if(!point && dec.expo < INT32_MAX)
						++dec.expo;
				} else if(dec.mant > (DECIMAL_MANT_MAX - 10) / 10) {
					if(c2 >= '5')
						++dec.mant;
					if(!point && dec.expo < INT32_MAX)
						++dec.expo;
					limit = true;
				} else {
					if(point && dec.expo > INT32_MIN)
						--dec.expo;
					dec.mant = dec.mant * 10 + (uint64_t) (c2 - '0');
				}
			} else if(c2 == '.') {
				ReadCharRaw(context);
				if(point)
					context.Error(1, "Unexpected second decimal point.");
				point = true;
			} else {
				break;
			}
		}

		// make sure that the number isn't just a decimal point
		if(starts_with_point && dec.expo == VDATA_DECIMAL_SHIFT) {
			context.Error(1, "Unexpected decimal point.");
		}

		// try to read exponent
		int c2 = PeekCharRaw(context);
		if(c2 == 'e' || c2 == 'E') {
			ReadCharRaw(context);

			// read sign if present
			c2 = PeekChar(context);
			bool expo_negative = (c2 == '-');
			if(c2 == '+' || c2 == '-') {
				ReadCharRaw(context);
			}

			// make sure that there is at least one digit
			char c = PeekChar(context);
			if(!(c >= '0' && c <= '9')) {
				context.Error(1, "Expected exponent");
			}

			// read exponent
			uint32_t expo = 0;
			for( ; ; ) {
				int c = PeekCharRaw(context);
				if(c >= '0' && c <= '9') {
					ReadCharRaw(context);
					if(expo > (INT32_MAX - 10) / 10) {
						expo = INT32_MAX;
						// some people just like to break things
						for( ; ; ) {
							int c = PeekCharRaw(context);
							if(c >= '0' && c <= '9') {
								ReadCharRaw(context);
							} else {
								break;
							}
						}
						break;
					}
					expo = expo * 10 + (uint32_t) (c - '0');
				} else {
					break;
				}
			}
			if(expo_negative) {
				if(expo > (uint32_t) dec.expo - (uint32_t) INT32_MIN)
					dec.expo = INT32_MIN;
				else
					dec.expo -= expo;
			} else {
				if(expo > (uint32_t) INT32_MAX - (uint32_t) dec.expo)
					dec.expo = INT32_MAX;
				else
					dec.expo += expo;
			}

			// convert to floating point
			data = FromDecimal(dec);

		} else {
			uint64_t lim = (dec.negative)? (uint64_t) INT64_MAX + 1 : (uint64_t) INT64_MAX;
			if(!point && dec.expo == VDATA_DECIMAL_SHIFT && dec.mant <= lim) {
				data = (dec.negative)? -(int64_t) dec.mant : (int64_t) dec.mant;
			} else {
				data = FromDecimal(dec);
			}
		}

	} else {
		context.Error(0, "Expected keyword or number");
	}

}

inline void WriteInt(WriteContext &context, int64_t value) {

	// check the sign
	bool negative = (value < 0);
	uint64_t value2 = (negative)? -value : value;

	// convert to string
	char buf[20], *buf2 = buf + sizeof(buf);
	do {
		*(--buf2) = (char) ('0' + (uint32_t) (value2 % 10));
		value2 /= 10;
	} while(value2 != 0);
	if(negative)
		*(--buf2) = '-';
	WriteData(context, buf2, buf + sizeof(buf) - buf2);

}

// ==================== Floating point helpers ====================

inline void WriteFloat(WriteContext &context, real_t value) {
	Decimal dec = ToDecimal(value, context.format.precision);
	switch(dec.type) {
		case FLOATTYPE_NORMAL: {
			uint64_t mant = dec.mant;
			if(mant == 0) {
				if(dec.negative) {
					WriteData(context, "-0.0", 4);
				} else {
					WriteData(context, "0.0", 3);
				}
				return;
			}
			int32_t expo = dec.expo - VDATA_DECIMAL_SHIFT + (int32_t) (context.format.precision - 1);
			assert(expo >= -360 && expo <= 360);
			int32_t decimals;
			char buffer[DECIMAL_BUFFER_SIZE];
			char *ptr = buffer + DECIMAL_BUFFER_SIZE;
			if(expo >= -3 && expo <= 8) {
				decimals = (int32_t) context.format.precision - 1 - expo;
			} else {
				int32_t engshift;
				engshift = (context.format.engineering)? (expo + 360) % 3 : 0; // mod does not work properly for negative numbers
				expo -= engshift;
				bool expo_negative = (expo < 0);
				uint32_t expo_abs = (expo_negative)? -expo : expo;
				do {
					*(--ptr) = (char) ('0' + (uint32_t) (expo_abs % 10));
					expo_abs /= 10;
				} while(expo_abs != 0);
				*(--ptr) = (expo_negative)? '-' : '+';
				*(--ptr) = 'e';
				decimals = (int32_t) context.format.precision - 1 - engshift;
			}
			if(decimals > 0) {
				while(decimals != 1) {
					--decimals;
					uint32_t digit = (uint32_t) (mant % 10);
					mant /= 10;
					if(digit != 0) {
						*(--ptr) = (char) ('0' + digit);
						break;
					}
				}
				while(decimals != 0) {
					--decimals;
					*(--ptr) = (char) ('0' + (uint32_t) (mant % 10));
					mant /= 10;
				}
			} else {
				*(--ptr) = '0';
			}
			*(--ptr) = '.';
			while(decimals != 0) {
				++decimals;
				*(--ptr) = '0';
			}
			do {
				*(--ptr) = (char) ('0' + (uint32_t) (mant % 10));
				mant /= 10;
			} while(mant != 0);
			if(dec.negative)
				*(--ptr) = '-';
			assert(ptr >= buffer);
			WriteData(context, ptr, buffer + DECIMAL_BUFFER_SIZE - ptr);
			break;
		}
		case FLOATTYPE_INF: {
			if(dec.negative) {
				WriteData(context, "-inf", 4);
			} else {
				WriteData(context, "inf", 3);
			}
			break;
		}
		case FLOATTYPE_NAN: {
			WriteData(context, "nan", 3);
			break;
		}
	}
}

// ==================== JSON Reading ====================

void ReadVData(ReadContext &context, VData &data) {

	// read first character
	char c = PeekChar(context);

	// what type is it?
	if(c == '{') { // dict

		ReadCharRaw(context);

		// create dict
		VData::Dict& ref = data.NewDict();
		SkipWhitespace(context);

		// read contents
		for( ; ; ) {

			// read the key
			c = ReadChar(context);
			if(c == '}')
				break;
			if(c != '"')
				context.Error((c == EOF)? 0 : 1, "Expected '\"' or '}'");
			std::string key;
			ReadString(context, key);
			SkipWhitespace(context);

			// read the colon
			c = ReadChar(context);
			if(c != ':')
				context.Error((c == EOF)? 0 : 1, "Expected ':'");
			SkipWhitespace(context);

			// read one value
			ref.EmplaceBack(StringRegistry::NewTag(key), nullptr);
			ReadVData(context, ref.Back().Value());

			// is there more?
			c = ReadChar(context);
			if(c == '}')
				break;
			if(c != ',')
				context.Error((c == EOF)? 0 : 1, "Expected ',' or '}'");
			SkipWhitespace(context);

		}

	} else if(c == '[') { // list

		ReadCharRaw(context);

		// create list
		VData::List &ref = data.NewList();
		SkipWhitespace(context);

		// read contents
		for( ; ; ) {

			// is this the end?
			if(PeekChar(context) == ']') {
				ReadChar(context);
				break;
			}

			// read one value
			ref.emplace_back(nullptr);
			ReadVData(context, ref.back());

			// is there more?
			c = ReadChar(context);
			if(c == ']')
				break;
			if(c != ',')
				context.Error((c == EOF)? 0 : 1, "Expected ',' or ']'");
			SkipWhitespace(context);

		}

	} else if(c == '"') { // string

		ReadCharRaw(context);
		ReadString(context, data.NewString());

	} else if(c == '+' || c == '-' || c == '.' || (c >= '0' && c <= '9') || (c >= 'a' && c <= 'z')) { // number or name

		ReadNameOrNumber(context, data);

	} else {
		ReadCharRaw(context);
		context.Error(1, "Unexpected character");
	}

	SkipWhitespace(context);

}

// ==================== JSON Writing ====================

void WriteVData(WriteContext &context, const VData &data) {
	switch(data.GetType()) {
		case VDATA_NULL: {
			WriteData(context, "null", 4);
			break;
		}
		case VDATA_BOOL: {
			if(data.AsBool())
				WriteData(context, "true", 4);
			else
				WriteData(context, "false", 5);
			break;
		}
		case VDATA_INT: {
			WriteInt(context, data.AsInt());
			break;
		}
		case VDATA_FLOAT: {
			WriteFloat(context, data.AsFloat());
			break;
		}
		case VDATA_STRING: {
			WriteString(context, data.AsString());
			break;
		}
		case VDATA_LIST: {
			WriteChar(context, '[');
			const VData::List &ref = data.AsList();
			if(ref.size() != 0) {
				++context.indent;
				WriteLine(context, false);
				WriteVData(context, ref[0]);
				for(size_t i = 1; i < ref.size(); ++i) {
					WriteChar(context, ',');
					WriteLine(context, true);
					WriteVData(context, ref[i]);
				}
				--context.indent;
				WriteLine(context, false);
			}
			WriteChar(context, ']');
			break;
		}
		case VDATA_DICT: {
			WriteChar(context, '{');
			const VData::Dict &ref = data.AsDict();
			if(ref.GetSize() != 0) {
				++context.indent;
				WriteLine(context, false);
				WriteString(context, StringRegistry::GetString(ref[0].Key()));
				WriteData(context, ": ", 2);
				WriteVData(context, ref[0].Value());
				for(index_t i = 1; i < ref.GetSize(); ++i) {
					WriteChar(context, ',');
					WriteLine(context, true);
					WriteString(context, StringRegistry::GetString(ref[i].Key()));
					WriteData(context, ": ", 2);
					WriteVData(context, ref[i].Value());
				}
				--context.indent;
				WriteLine(context, false);
			}
			WriteChar(context, '}');
			break;
		}
	}
}

// ==================== Simple public interface ====================

void FromStream(VData &data, std::streambuf *stream) {
	ReadContext context;
	context.stream = stream;
	context.line = 0;
	context.column = 0;
	ReadVData(context, data);
	if(context.stream->sgetc() != EOF)
		context.Error(0, "Extra data at end of input.");
}

void FromString(VData &data, const std::string &str) {
	std::stringbuf stream(str);
	FromStream(data, &stream);
}

void FromFile(VData &data, const std::string &filename) {
	std::filebuf stream;
	if(stream.open(filename, std::ios_base::in | std::ios_base::binary) == NULL)
		throw std::runtime_error("Could not open file '" + filename + "' for reading.");
	FromStream(data, &stream);
}

void ToStream(const VData &data, std::streambuf *stream, const Format &format) {
	WriteContext context;
	context.stream = stream;
	context.format = format;
	context.format.precision = clamp<uint32_t>(context.format.precision, 2, 18);
	context.indent = 0;
	WriteVData(context, data);
}

void ToString(const VData &data, std::string &str, const Format &format) {
	std::stringbuf stream;
	ToStream(data, &stream, format);
	str = stream.str();
}

void ToFile(const VData &data, const std::string &filename, const Format &format) {
	std::filebuf stream;
	if(stream.open(filename, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc) == NULL)
		throw std::runtime_error("Could not open file '" + filename + "' for writing.");
	ToStream(data, &stream, format);
}

}
