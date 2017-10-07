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
#include "MiscMath.h"

struct Color {
	float r, g, b, a;
	inline Color() {}
	inline Color(float r, float g, float b, float a = 1.0f) : r(r), g(g), b(b), a(a) {}
	//inline Color(const QColor &color) : r(color.redF()), g(color.greenF()), b(color.blueF()), a(color.alphaF()) {}
	//inline operator QColor() const { return QColor::fromRgbF(clamp(r, 0.0f, 1.0f), clamp(g, 0.0f, 1.0f), clamp(b, 0.0f, 1.0f), clamp(a, 0.0f, 1.0f)); }
	inline operator float*() { return &r; }
	inline uint32_t ToUint32() {
		uint32_t rr = rint32(clamp<float>(r * 255.0f, 0.0f, 255.0f));
		uint32_t gg = rint32(clamp<float>(g * 255.0f, 0.0f, 255.0f));
		uint32_t bb = rint32(clamp<float>(b * 255.0f, 0.0f, 255.0f));
		uint32_t aa = rint32(clamp<float>(a * 255.0f, 0.0f, 255.0f));
		return (aa << 24) | (rr << 16) | (gg << 8) | bb;
	}
};

static_assert(sizeof(Color) == sizeof(float) * 4, "Invalid Color struct!");

inline Color ColorPremultiply(const Color &color) {
	return Color(color.r * color.a, color.g * color.a, color.b * color.a, color.a);
}

inline Color ColorMix(const Color &col1, const Color &col2, float frac) {
	return Color(
		col1.r + (col2.r - col1.r) * frac,
		col1.g + (col2.g - col1.g) * frac,
		col1.b + (col2.b - col1.b) * frac,
		col1.a + (col2.a - col1.a) * frac
	);
}

inline Color ColorBlend(const Color &col1, const Color &col2) {
	return Color(
		col1.r + (col2.r - col1.r) * col2.a,
		col1.g + (col2.g - col1.g) * col2.a,
		col1.b + (col2.b - col1.b) * col2.a,
		col1.a * (1.0f - col2.a) + col2.a
	);
}

inline Color ColorBlendPremultiplied(const Color &col1, const Color &col2) {
	return Color(
		col1.r * (1.0f - col2.a) + col2.r,
		col1.g * (1.0f - col2.a) + col2.g,
		col1.b * (1.0f - col2.a) + col2.b,
		col1.a * (1.0f - col2.a) + col2.a
	);
}

inline float ToSRGB(float x) {
	return (x < 0.04045f / 12.92f)? 12.92f * x : 1.055f * powf(x, 1.0f / 2.4f) - 0.055f;
}

inline float FromSRGB(float x) {
	return (x < 0.04045f)? (1.0f / 12.92f) * x : powf((x + 0.055f) / 1.055f, 2.4f);
}

inline Color ToSRGB(const Color &color, bool dark = true) {
	return Color(ToSRGB(color.r), ToSRGB(color.g), ToSRGB(color.b), (dark)? ToSRGB(color.a) : 1.0f - ToSRGB(1.0f - color.a));
}

inline Color FromSRGB(const Color &color, bool dark = true) {
	return Color(FromSRGB(color.r), FromSRGB(color.g), FromSRGB(color.b), (dark)? FromSRGB(color.a) : 1.0f - FromSRGB(1.0f - color.a));
}

/*inline GLvec4 ColorToVec(const QColor &color) {
	return GLvec4(color.redF(), color.greenF(), color.blueF(), color.alphaF());
}*/

/*inline GLvec4 ColorToVec(const Color &color) {
	return GLvec4(color.r, color.g, color.b, color.a);
}*/

Color LumaInvert(const Color &color, bool invert = true);
