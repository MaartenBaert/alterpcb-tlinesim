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

#include "QLineEditSmall.h"

QLineEditSmall::QLineEditSmall(QWidget *parent, int size_factor)
	: QLineEdit(parent) {
	m_size_factor = size_factor;
}

QLineEditSmall::QLineEditSmall(const QString &contents, QWidget *parent, int size_factor)
	: QLineEdit(contents, parent) {
	m_size_factor = size_factor;
}

QLineEditSmall::~QLineEditSmall() {
	// nothing
}

QSize QLineEditSmall::sizeHint() const {
	QSize size = QLineEdit::sizeHint();
	return QSize(size.width() * m_size_factor / 100, size.height());
}
