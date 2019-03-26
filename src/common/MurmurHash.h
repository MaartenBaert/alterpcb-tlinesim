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

namespace MurmurHash {

inline hash_t RotateLeft(hash_t x, uint32_t r) {
	return (x << r) | (x >> (32 - r));
}

inline hash_t HashData(hash_t hash, uint32_t data) {
	data *= 0xcc9e2d51;
	data = RotateLeft(data, 15);
	data *= 0x1b873593;
	hash ^= data;
	hash = RotateLeft(hash, 13);
	hash = hash * 5 + 0xe6546b64;
	return hash;
}

inline hash_t HashData(hash_t hash, uint64_t data) {
	hash = HashData(hash, (uint32_t) data);
	hash = HashData(hash, (uint32_t) (data >> 32));
	return hash;
}

inline hash_t HashData(hash_t hash, const void *data, size_t size) {
	size_t word_count = size / 4;
	for(size_t i = 0; i < word_count; ++i) {
		uint32_t word; // memcpy is required due to strict aliasing, the compiler will optimize this
		memcpy(&word, (const uint8_t*) data + i * 4, 4);
		hash = HashData(hash, word);
	}
	const uint8_t *remainder = (const uint8_t*) data + word_count * 4;
	uint32_t value = 0;
	switch(size & 3) {
		case 3: value |= (uint32_t) remainder[2] << 16; // fallthrough
		case 2: value |= (uint32_t) remainder[1] << 8; // fallthrough
		case 1: value |= (uint32_t) remainder[0];
		hash = HashData(hash, value);
	}
	return hash;
}

inline hash_t HashFinish(hash_t hash) {
	hash ^= hash >> 16;
	hash *= 0x85ebca6b;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35;
	hash ^= hash >> 16;
	return hash;
}

inline hash_t HashTruncate(hash_t hash, uint32_t bits) {
	assert(bits < 32);
	return hash & (((uint32_t) 1 << bits) - 1);
}

}
