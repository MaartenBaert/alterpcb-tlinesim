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

#include <memory>
#include <utility>

/*
The Cow (Copy-On-Write) class provides an easy way to enforce copy-on-write semantics. Cow<T> manages a shared_ptr<T>
such that you can normally only get a const reference to the object using the Ref() function. In order to get a
non-const reference, you have to call Unique() which makes a copy of the object if necessary.

Currently polymorphic objects are not supported, but support could be added with a virtual Clone() function.
*/

template<typename T>
class Cow {

	template<typename U>
	friend class Cow;

private:
	std::shared_ptr<T> m_data;

public:
	inline Cow() {}

	// These are intended only for use with std::make_shared(), because we can't have a variadic constructor here.
	// Don't use these to store existing pointers unless you know what you are doing.
	inline explicit Cow(const std::shared_ptr<T> &data) : m_data(data) {}
	inline explicit Cow(std::shared_ptr<T> &&data) : m_data(std::move(data)) {}

	// default copy and assignment
	Cow(const Cow&) = default;
	Cow(Cow&&) = default;
	Cow& operator=(const Cow&) = default;
	Cow& operator=(Cow&&) = default;

	// null pointer creation and assignment
	inline Cow(std::nullptr_t) {}
	inline Cow& operator=(std::nullptr_t) {
		m_data.reset();
	}

	// The preferred way to create new objects if you already have the pointer. It behaves like calling MakeCow()
	// followed by UniqueBypass() to get a reference to the object you just created.
	template<typename... Args>
	inline T& New(Args&&... args) {
		m_data = std::make_shared<T>(std::forward<Args>(args)...);
		return *m_data;
	}

	// Sets the pointer to NULL, possibly deleting the old object.
	inline void Reset() {
		m_data.reset();
	}

	// Returns whether the pointer is NULL.
	inline bool IsNull() const {
		return m_data == NULL;
	}

	// Returns a constant reference to the object. The object can be read but not modified.
	inline const T& Ref() const {
		assert(!IsNull());
		return *m_data;
	}

	// Returns a modifyable reference to the object. If the pointer is not unique, this creates a copy of the object.
	// Copying can be costly, so only use this function when it is actually needed.
	inline T& Unique() {
		assert(!IsNull());
		if(!m_data.unique())
			m_data = std::make_shared<T>(*m_data);
		return *m_data;
	}

	// This function bypasses the uniqueness check for performance. This should only be used in situations where you
	// already know that the pointer will be unique, e.g. right after you have created a new object.
	inline T& UniqueBypass() {
		assert(!IsNull());
		assert(m_data.unique());
		return *m_data;
	}

};

// This is a convenience function. If you need to modify the object after creation, you should use New() instead.
template<typename T, typename... Args>
inline Cow<T> MakeCow(Args&&... args) {
	return Cow<T>(std::make_shared<T>(std::forward<Args>(args)...));
}
