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

#include "Color.h"

void LumaInvert2(float &r, float &g, float &b) {
	if(r == b) {
		r = 1.0f - r;
		g = 1.0f - g;
		b = 1.0f - b;
	} else {
		float sum = r + g + b, invsum = 3.0f - sum;
		float r3 = 3.0f * r - sum;
		float g3 = 3.0f * g - sum;
		float b3 = 3.0f * b - sum;
		float num, den;
		if(sum * (b - r) < b3) {
			if(sum * (r - b) < r3) {
				num = sum;
				den = invsum;
			} else {
				num = r3;
				den = -b3;
			}
		} else {
			if(sum * (r - b) < r3) {
				num = -b3;
				den = r3;
			} else {
				num = invsum;
				den = sum;
			}
		}
		r = (invsum * den + r3 * num) / (3.0f * den);
		g = (invsum * den + g3 * num) / (3.0f * den);
		b = (invsum * den + b3 * num) / (3.0f * den);
	}
}

Color LumaInvert(const Color &color, bool invert) {
	if(!invert)
		return color;
	float r = color.r, g = color.g, b = color.b;
	if(r > g && r > b) {
		if(g > b) { // r > g > b
			LumaInvert2(r, g, b);
		} else { // r > b > g
			LumaInvert2(r, b, g);
		}
	} else if(g > b) {
		if(r > b) { // g > r > b
			LumaInvert2(g, r, b);
		} else { // g > b > r
			LumaInvert2(g, b, r);
		}
	} else {
		if(r > g) { // b > r > g
			LumaInvert2(b, r, g);
		} else { // b > g > r
			LumaInvert2(b, g, r);
		}
	}
	return Color(r, g, b, color.a);
}
