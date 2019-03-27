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
#include "MurmurHash.h"

#include <utility>
#include <vector>

template<typename T, class Hasher>
class HashTable {

private:
	struct Bucket {
		size_t m_hashlink;
		size_t m_size;
		inline Bucket() {
			m_hashlink = INDEX_NONE;
			m_size = 0;
		}
	};
	struct Entry : public T {
		size_t m_hashlink;
		template<typename... Args>
		inline Entry(Args&&... args)
			: T(std::forward<Args>(args)...) {}
	};

private:
	constexpr static hash_t HASH_SEED1 = 0x5d02e5ec, HASH_SEED2 = 0xbfcb8a4e;
	static Bucket SENTINEL_BUCKETS[2];

private:
	Hasher m_hasher;
	Bucket *m_buckets;
	size_t m_num_buckets;
	std::vector<Entry> m_data;
	uint32_t m_hash_bits;

public:
	HashTable(const Hasher &hasher = Hasher())
		: m_hasher(hasher) {
		m_buckets = SENTINEL_BUCKETS;
		m_num_buckets = 0;
		m_hash_bits = 0;
	}
	~HashTable() {
		if(m_buckets != SENTINEL_BUCKETS)
			delete[] m_buckets;
	}

	// default copy and assignment
	HashTable(const HashTable&) = default;
	HashTable(HashTable&&) noexcept = default;
	HashTable& operator=(const HashTable&) = default;
	HashTable& operator=(HashTable&&) noexcept = default;

	void Free() noexcept {
		if(m_buckets != SENTINEL_BUCKETS) {
			delete[] m_buckets;
			m_buckets = SENTINEL_BUCKETS;
		}
		m_data = std::vector<Entry>();
		m_num_buckets = 0;
		m_hash_bits = 0;
	}

	void Clear() noexcept {
		m_data.clear();
		for(size_t i = 0; i < m_num_buckets; ++i) {
			m_buckets[i].m_hashlink = INDEX_NONE;
			m_buckets[i].m_size = 0;
		}
	}

	void Reserve(size_t capacity) {
		m_data.reserve(capacity);
		if(capacity > m_num_buckets) {
			m_hash_bits = std::max<uint32_t>(2, CeilLog2(capacity)) - 1;
			Rehash();
		}
	}

	std::pair<size_t, bool> TryPushBack(const T &value) {
		hash_t hash1 = Hash1(value);
		size_t res1 = FindAt(hash1, value);
		if(res1 != INDEX_NONE)
			return std::make_pair(res1, false);
		hash_t hash2 = Hash2(value);
		size_t res2 = FindAt(hash2, value);
		if(res2 != INDEX_NONE)
			return std::make_pair(res2, false);
		if(MaybeRehash()) {
			hash1 = Hash1(value);
			hash2 = Hash2(value);
		}
		m_data.push_back(value);
		size_t res = m_data.size() - 1;
		AddHash(hash1, hash2, res);
		return std::make_pair(res, true);
	}

	std::pair<size_t, bool> TryPushBack(T &&value) {
		hash_t hash1 = Hash1(value);
		size_t res1 = FindAt(hash1, value);
		if(res1 != INDEX_NONE)
			return std::make_pair(res1, false);
		hash_t hash2 = Hash2(value);
		size_t res2 = FindAt(hash2, value);
		if(res2 != INDEX_NONE)
			return std::make_pair(res2, false);
		if(MaybeRehash()) {
			hash1 = Hash1(value);
			hash2 = Hash2(value);
		}
		m_data.push_back(std::move(value));
		size_t res = m_data.size() - 1;
		AddHash(hash1, hash2, res);
		return std::make_pair(res, true);
	}

	template<typename K, typename... Args>
	std::pair<size_t, bool> TryEmplaceBack(const K &key, Args&&... args) {
		hash_t hash1 = Hash1(key);
		size_t res1 = FindAt(hash1, key);
		if(res1 != INDEX_NONE)
			return std::make_pair(res1, false);
		hash_t hash2 = Hash2(key);
		size_t res2 = FindAt(hash2, key);
		if(res2 != INDEX_NONE)
			return std::make_pair(res2, false);
		if(MaybeRehash()) {
			hash1 = Hash1(key);
			hash2 = Hash2(key);
		}
		m_data.emplace_back(std::forward<Args>(args)...);
		size_t res = m_data.size() - 1;
		AddHash(hash1, hash2, res);
		return std::make_pair(res, true);
	}

	size_t PushBack(const T &value) {
		MaybeRehash();
		m_data.push_back(value);
		size_t res = m_data.size() - 1;
		AddHash(Hash1(m_data[res]), Hash2(m_data[res]), res);
		return res;
	}

	size_t PushBack(T &&value) {
		MaybeRehash();
		m_data.push_back(std::move(value));
		size_t res = m_data.size() - 1;
		AddHash(Hash1(m_data[res]), Hash2(m_data[res]), res);
		return res;
	}

	template<typename... Args>
	size_t EmplaceBack(Args&&... args) {
		MaybeRehash();
		m_data.emplace_back(std::forward<Args>(args)...);
		size_t res = m_data.size() - 1;
		AddHash(Hash1(m_data[res]), Hash2(m_data[res]), res);
		return res;
	}

	template<typename K>
	size_t Find(const K &key) const noexcept {
		size_t res1 = FindAt(Hash1(key), key);
		if(res1 != INDEX_NONE)
			return res1;
		size_t res2 = FindAt(Hash2(key), key);
		return res2;
	}

private:
	template<typename K>
	hash_t Hash1(const K &key) const noexcept {
		hash_t hash = MurmurHash::HashFinish(m_hasher.Hash(HASH_SEED1, key));
		return MurmurHash::HashTruncate(hash, m_hash_bits);
	}

	template<typename K>
	hash_t Hash2(const K &key) const noexcept {
		hash_t hash = MurmurHash::HashFinish(m_hasher.Hash(HASH_SEED2, key));
		return MurmurHash::HashTruncate(hash, m_hash_bits) | ((hash_t) 1 << m_hash_bits);
	}

	template<typename K>
	size_t FindAt(hash_t hash, const K &key) const noexcept {
		for(size_t i = m_buckets[hash].m_hashlink; i != INDEX_NONE; i = m_data[i].m_hashlink) {
			if(m_hasher.Equal(m_data[i], key))
				return i;
		}
		return INDEX_NONE;
	}

	void AddHash(hash_t hash1, hash_t hash2, size_t index) noexcept {
		hash_t hash = (m_buckets[hash1].m_size <= m_buckets[hash2].m_size)? hash1 : hash2;
		m_data[index].m_hashlink = m_buckets[hash].m_hashlink;
		m_buckets[hash].m_hashlink = index;
		++m_buckets[hash].m_size;
	}

	bool MaybeRehash() {
		if(m_data.size() < m_num_buckets)
			return false;
		++m_hash_bits;
		Rehash();
		return true;
	}

	void Rehash() {
		size_t new_num_buckets = (size_t) 1 << (m_hash_bits + 1);
		Bucket *new_buckets = new Bucket[new_num_buckets];
		if(new_buckets == NULL)
			throw std::bad_alloc();
		if(m_buckets != SENTINEL_BUCKETS)
			delete[] m_buckets;
		m_buckets = new_buckets;
		m_num_buckets = new_num_buckets;
		for(size_t i = 0; i < m_data.size(); ++i) {
			hash_t hash1 = Hash1(m_data[i]);
			hash_t hash2 = Hash2(m_data[i]);
			AddHash(hash1, hash2, i);
		}
	}

public:
	inline size_t GetSize() const {
		return m_data.size();
	}
	inline T& operator[](size_t i) {
		assert(i < m_data.size());
		return m_data[i];
	}
	inline const T& operator[](size_t i) const {
		assert(i < m_data.size());
		return m_data[i];
	}
	inline T& Front() {
		assert(!m_data.empty());
		return m_data.front();
	}
	inline const T& Front() const {
		assert(!m_data.empty());
		return m_data.front();
	}
	inline T& Back() {
		assert(!m_data.empty());
		return m_data.back();
	}
	inline const T& Back() const {
		assert(!m_data.empty());
		return m_data.back();
	}

};

template<typename T, class Hasher>
typename HashTable<T, Hasher>::Bucket HashTable<T, Hasher>::SENTINEL_BUCKETS[2];
