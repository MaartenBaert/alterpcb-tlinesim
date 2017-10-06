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

#include "HashTable.h"

#include <string>

/*
AlterPCB uses short strings in lots of places, most importantly for parameters. String lookup is slow, so instead all
these strings are replaced with a string tag (stringtag_t). This is really just a number that identifies, nothing more.
These tags are tracked by StringRegistry.

Currently this is not thread-safe.
*/

class StringRegistryEntry {

private:
	std::string m_string;

public:
	inline StringRegistryEntry(const std::string &str) : m_string(str) {}
	inline StringRegistryEntry(std::string &&str) : m_string(std::move(str)) {}

	// default copy and assignment
	StringRegistryEntry(const StringRegistryEntry&) = default;
	StringRegistryEntry(StringRegistryEntry&&) = default;
	StringRegistryEntry& operator=(const StringRegistryEntry&) = default;
	StringRegistryEntry& operator=(StringRegistryEntry&&) = default;

	inline const std::string& GetString() const {
		return m_string;
	}

};

struct StringRegistryHasher {
	inline static bool Equal(const StringRegistryEntry &a, const StringRegistryEntry &b) {
		return a.GetString() == b.GetString();
	}
	inline static bool Equal(const StringRegistryEntry &a, const std::string &b) {
		return a.GetString() == b;
	}
	inline static bool Equal(const StringRegistryEntry &a, const char *b) {
		return a.GetString() == b;
	}
	inline static hash_t Hash(hash_t hash, const StringRegistryEntry &value) {
		return MurmurHash::HashData(hash, value.GetString().data(), value.GetString().size());
	}
	inline static hash_t Hash(hash_t hash, const std::string &value) {
		return MurmurHash::HashData(hash, value.data(), value.size());
	}
	inline static hash_t Hash(hash_t hash, const char *value) {
		return MurmurHash::HashData(hash, value, strlen(value));
	}
};

class StringRegistry {

private:
	static StringRegistry *s_instance;

private:
	HashTable<StringRegistryEntry, StringRegistryHasher> m_entries;

public:
	inline StringRegistry() : m_entries(1024) {
		s_instance = this;
	}
	inline ~StringRegistry() {
		s_instance = NULL;
	}

	// noncopyable
	StringRegistry(const StringRegistry&) = delete;
	StringRegistry& operator=(const StringRegistry&) = delete;

	// Converts a string to a tag, adding the string to the registry if necessary.
	inline static stringtag_t NewTag(const std::string &str) {
		assert(s_instance != NULL);
		return s_instance->m_entries.TryEmplaceBack(str, str).first;
	}
	inline static stringtag_t NewTag(std::string &&str) {
		assert(s_instance != NULL);
		return s_instance->m_entries.TryEmplaceBack(str, std::move(str)).first;
	}
	inline static stringtag_t NewTag(const char *str) {
		assert(s_instance != NULL);
		return s_instance->m_entries.TryEmplaceBack(str, str).first;
	}

	// Converts a string to a tag, but if the string is not already registered, the string is not added and instead the
	// function returns STRINGTAG_NONE. Useful for things like key lookup.
	inline static stringtag_t FindTag(const std::string &str) {
		assert(s_instance != NULL);
		return s_instance->m_entries.Find(str);
	}
	inline static stringtag_t FindTag(const char *str) {
		assert(s_instance != NULL);
		return s_instance->m_entries.Find(str);
	}

	// Converts a tag back to a string. the tag *must* be valid, there is no error checking.
	inline static const std::string& GetString(stringtag_t tag) {
		assert(s_instance != NULL);
		return s_instance->m_entries[tag].GetString();
	}

};

// convenience functions (shorter name)
inline stringtag_t SRNewTag(const std::string &str) { return StringRegistry::NewTag(str); }
inline stringtag_t SRNewTag(std::string &&str) { return StringRegistry::NewTag(str); }
inline stringtag_t SRNewTag(const char *str) { return StringRegistry::NewTag(str); }
inline stringtag_t SRFindTag(const std::string &str) { return StringRegistry::FindTag(str); }
inline stringtag_t SRFindTag(const char *str) { return StringRegistry::FindTag(str); }
inline const std::string& SRGetString(stringtag_t tag) { return StringRegistry::GetString(tag); }
