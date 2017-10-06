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

#include "Qt.h"

inline void GroupEnabled(std::initializer_list<QAction*> actions, bool enabled) {
	for(QAction *a : actions) {
		a->setEnabled(enabled);
	}
}

inline void GroupEnabled(std::initializer_list<QWidget*> widgets, bool enabled) {
	for(QWidget *w : widgets) {
		w->setEnabled(enabled);
	}
}

inline void GroupVisible(std::initializer_list<QWidget*> widgets, bool visible) {
	for(QWidget *w : widgets) {
		w->setVisible(visible);
	}
}

inline void MultiGroupVisible(std::initializer_list<std::pair<std::initializer_list<QWidget*>, bool> > conditions) {
	// Qt updates the layout every time something is made visible or invisible, so the order is important.
	// First hide everything that needs to be hidden, then show everything that needs to be shown.
	// If the order is wrong, Qt will make the widget larger than necessary.
	for(auto &c : conditions) {
		if(!c.second)
			GroupVisible(c.first, false);
	}
	for(auto &c : conditions) {
		if(c.second)
			GroupVisible(c.first, true);
	}
}
