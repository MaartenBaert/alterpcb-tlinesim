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

template<typename T>
inline T vmin(T v) {
	return v;
}
template<typename T, typename... Args>
inline T vmin(T v1, T v2, Args&&... args) {
	return vmin(std::min(v1, v2), std::forward<Args>(args)...);
}

template<typename T>
inline T vmax(T v) {
	return v;
}
template<typename T, typename... Args>
inline T vmax(T v1, T v2, Args&&... args) {
	return vmax(std::max(v1, v2), std::forward<Args>(args)...);
}

template<typename T, typename U>
inline T clamp(U v, T lo, T hi) {
	assert(lo <= hi);
	if(v < (U) lo)
		return lo;
	if(v > (U) hi)
		return hi;
	return (T) v;
}
template<> inline float clamp<float, float>(float v, float lo, float hi) {
	assert(!(lo > hi)); // nan ok
	return std::min(std::max(v, lo), hi);
}
template<> inline double clamp<double, double>(double v, double lo, double hi) {
	assert(!(lo > hi)); // nan ok
	return std::min(std::max(v, lo), hi);
}

template<typename F>
inline F square(F x) {
	return x * x;
}
template<typename F>
inline F cube(F x) {
	return x * x * x;
}

template<typename F, typename G>
F lerp(F n0, F n1, G fx) {
	return n0 + (n1 - n0) * fx;
}

template<typename F>
inline int32_t rint32(F x) {
	return (sizeof(long int) >= sizeof(int32_t))? (int32_t) lrint(x) : (int32_t) llrint(x);
}
template<typename F>
inline int64_t rint64(F x) {
	return (sizeof(long int) >= sizeof(int64_t))? (int64_t) lrint(x) : (int64_t) llrint(x);
}
template<typename F>
inline int rinti(F x) {
	return (sizeof(long int) >= sizeof(int))? (int) lrint(x) : (int) llrint(x);
}
template<typename F>
inline ptrdiff_t rintp(F x) {
	return (sizeof(long int) >= sizeof(ptrdiff_t))? (ptrdiff_t) lrint(x) : (ptrdiff_t) llrint(x);
}

#ifdef _WIN32
template<typename F>
inline F exp10(F x) {
	return exp(x * (F) M_LN10);
}
#endif

template<typename F>
inline F GetEpsilon(F x) {
	return std::max(std::numeric_limits<real_t>::min(), std::numeric_limits<real_t>::epsilon() * fabs(x));
}
template<typename F>
inline F GetEpsilon(F x1, F x2) {
	return std::max(std::numeric_limits<real_t>::min(), std::numeric_limits<real_t>::epsilon() * std::max(fabs(x1), fabs(x2)));
}

template<typename F>
inline bool FinitePositive(F x) {
	return std::isfinite(x) && x > 0.0;
}
template<typename F>
inline bool FiniteNegative(F x) {
	return std::isfinite(x) && x < 0.0;
}
template<typename F>
inline bool FiniteNonPositive(F x) {
	return std::isfinite(x) && x <= 0.0;
}
template<typename F>
inline bool FiniteNonNegative(F x) {
	return std::isfinite(x) && x >= 0.0;
}
template<typename F>
inline bool FiniteGreater(F x, F y) {
	return std::isfinite(x) && std::isfinite(y) && x > y;
}
template<typename F>
inline bool FiniteLess(F x, F y) {
	return std::isfinite(x) && std::isfinite(y) && x < y;
}
template<typename F>
inline bool FiniteGreaterEqual(F x, F y) {
	return std::isfinite(x) && std::isfinite(y) && x >= y;
}
template<typename F>
inline bool FiniteLessEqual(F x, F y) {
	return std::isfinite(x) && std::isfinite(y) && x <= y;
}

template<typename F>
inline F ToRadians(F x) {
	return x * (F) (M_PI / 180.0);
}
template<typename F>
inline F ToDegrees(F x) {
	return x * (F) (180.0 / M_PI);
}

// IsPow2(0) = true
// IsPow2(1) = true
inline bool IsPow2(uint32_t x) {
	return ((x & (x - 1)) == 0);
}

// NextPow2(0) = 0
// NextPow2(1) = 1
inline uint32_t NextPow2(uint32_t x) {
	--x;
	x |= (x >>  1);
	x |= (x >>  2);
	x |= (x >>  4);
	x |= (x >>  8);
	x |= (x >> 16);
	return x + 1;
}

// FloorLog2(0) = 0
// FloorLog2(1) = 0
inline uint32_t FloorLog2(uint32_t x) {
	uint32_t i = 0;
	if(x >= ((uint32_t) 1 << 16)) { i += 16; x >>= 16; }
	if(x >= ((uint32_t) 1 <<  8)) { i +=  8; x >>=  8; }
	if(x >= ((uint32_t) 1 <<  4)) { i +=  4; x >>=  4; }
	if(x >= ((uint32_t) 1 <<  2)) { i +=  2; x >>=  2; }
	if(x >= ((uint32_t) 1 <<  1)) { i +=  1; }
	return i;
}
inline uint32_t FloorLog2(uint64_t x) {
	uint32_t i = 0;
	if(x >= ((uint64_t) 1 << 32)) { i += 32; x >>= 32; }
	if(x >= ((uint64_t) 1 << 16)) { i += 16; x >>= 16; }
	if(x >= ((uint64_t) 1 <<  8)) { i +=  8; x >>=  8; }
	if(x >= ((uint64_t) 1 <<  4)) { i +=  4; x >>=  4; }
	if(x >= ((uint64_t) 1 <<  2)) { i +=  2; x >>=  2; }
	if(x >= ((uint64_t) 1 <<  1)) { i +=  1; }
	return i;
}

// CeilLog2(0) = 0
// CeilLog2(1) = 0
inline uint32_t CeilLog2(uint32_t x) {
	return (x < 2)? 0 : FloorLog2(x - 1) + 1;
}
inline uint32_t CeilLog2(uint64_t x) {
	return (x < 2)? 0 : FloorLog2(x - 1) + 1;
}

template<typename T, typename U>
inline T MemCast(U value) {
	static_assert(sizeof(T) == sizeof(U), "Size of types does not match!");
	T res;
	memcpy(&res, &value, sizeof(T));
	return res;
}

// calculates (a * b) >> 64 without overflow
inline uint64_t FixMul64F(uint64_t a, uint64_t b) {
#if defined(__GNUC__) && defined(__SIZEOF_INT128__)
	typedef unsigned __int128 UINT128;
	static_assert(__SIZEOF_INT128__ == 16, "Size of 128-bit integer should be 16 bytes!");
	static_assert(sizeof(UINT128) == 16, "Size of 128-bit integer should be 16 bytes!");
	UINT128 c = (UINT128) a * (UINT128) b;
	return (uint64_t) (c >> 64);
#else
	uint32_t a0 = (uint32_t) a, a1 = (uint32_t) (a >> 32);
	uint32_t b0 = (uint32_t) b, b1 = (uint32_t) (b >> 32);
	uint64_t c00 = (uint64_t) a0 * (uint64_t) b0;
	uint64_t c01 = (uint64_t) a0 * (uint64_t) b1;
	uint64_t c10 = (uint64_t) a1 * (uint64_t) b0;
	uint64_t c11 = (uint64_t) a1 * (uint64_t) b1;
	uint64_t mid = (c00 >> 32) + (c01 & 0xffffffff) + (c10 & 0xffffffff);
	uint64_t hi = c11 + (c01 >> 32) + (c10 >> 32);
	return hi + (mid >> 32);
#endif
}

// calculates (a * b + (1 << 63)) >> 64 without overflow
inline uint64_t FixMul64R(uint64_t a, uint64_t b) {
#if defined(__GNUC__) && defined(__SIZEOF_INT128__)
	typedef unsigned __int128 UINT128;
	static_assert(__SIZEOF_INT128__ == 16, "Size of 128-bit integer should be 16 bytes!");
	static_assert(sizeof(UINT128) == 16, "Size of 128-bit integer should be 16 bytes!");
	UINT128 c = (UINT128) a * (UINT128) b;
	uint64_t lo = (uint64_t) c, hi = (uint64_t) (c >> 64);
	return hi + (lo >> 63);
#else
	uint32_t a0 = (uint32_t) a, a1 = (uint32_t) (a >> 32);
	uint32_t b0 = (uint32_t) b, b1 = (uint32_t) (b >> 32);
	uint64_t c00 = (uint64_t) a0 * (uint64_t) b0;
	uint64_t c01 = (uint64_t) a0 * (uint64_t) b1;
	uint64_t c10 = (uint64_t) a1 * (uint64_t) b0;
	uint64_t c11 = (uint64_t) a1 * (uint64_t) b1;
	uint64_t mid = (c00 >> 32) + (c01 & 0xffffffff) + (c10 & 0xffffffff);
	uint64_t hi = c11 + (c01 >> 32) + (c10 >> 32);
	return hi + ((mid + (UINT64_C(1) << 31)) >> 32);
#endif
}
