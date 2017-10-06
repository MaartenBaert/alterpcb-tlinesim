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

#include "NaturalSort.h"

// locale-independent alternative to std::isdigit()
inline bool IsDigit(char c) {
	return c >= '0' && c <= '9';
}

// This comparison function is used for 'natural' string sorting. It detects numbers in strings and sorts them in a way
// that makes sense to humans, by comparing the value rather than the binary representation. It does this by padding the
// numbers with zeros until they have the same length. Normally numbers are padded from the left, however if a number is
// followed by a decimal dot followed by a second number, the second number is considered a fractional part, and will be
// padded from the right. For example, when comparing the strings "X 1.2" and "X 12.345", the comparison will behave as
// if the first string was actually "X 01.200". If two numbers have the same numerical value, the one with the highest
// number of leading/trailing zeros is considered to be greater. For example, "12" < "012" and "1.2" < "1.20".
int NaturalStringCompare(const char *a, size_t a_len, const char *b, size_t b_len) {
	size_t ia = 0, ib = 0;
	while(ia < a_len && ib < b_len) {
		if(IsDigit(a[ia]) && IsDigit(b[ib])) { // both strings contain a number
			if(ia >= 2 && a[ia - 1] == '.' && IsDigit(a[ia - 2])) { // fractional part of previous number
				do { // compare the numbers
					if(a[ia] != b[ib]) {
						return (int)(unsigned char) a[ia] - (int)(unsigned char) b[ib];
					}
					++ia; ++ib;
				} while(ia < a_len && ib < b_len && IsDigit(a[ia]) && IsDigit(b[ib]));
				if(ia < a_len && IsDigit(a[ia])) // number in a is longer
					return 1;
				if(ib < b_len && IsDigit(b[ib])) // number in b is longer
					return -1;
			} else { // regular number
				size_t ja = ia + 1, jb = ib + 1;
				while(ja < a_len && IsDigit(a[ja])) { // find the first non-digit in a
					++ja;
				}
				while(jb < b_len && IsDigit(b[jb])) { // find the first non-digit in b
					++jb;
				}
				size_t na = ja - ia, nb = jb - ib;
				while(ja - ia > jb - ib) { // number in a is longer
					if(a[ia++] != '0') {
						return 1; // number in a has the highest value
					}
				}
				while(ja - ia < jb - ib) { // number in b is longer
					if(b[ib++] != '0') {
						return -1; // number in b has the highest value
					}
				}
				while(ia < ja) { // compare the numbers
					if(a[ia] != b[ib]) {
						return (int)(unsigned char) a[ia] - (int)(unsigned char) b[ib];
					}
					++ia; ++ib;
				}
				if(na != nb) { // if both numbers have the same value, the longest one is considered higher
					return (na > nb)? 1 : -1;
				}
			}
		} else {
			if(a[ia] != b[ib]) {
				return (int)(unsigned char) a[ia] - (int)(unsigned char) b[ib];
			}
			++ia; ++ib;
		}
	}
	return (ia < a_len)? 1 : (ib < b_len)? -1 : 0;
}
