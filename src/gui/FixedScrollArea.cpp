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

#include "FixedScrollArea.h"

FixedScrollArea::FixedScrollArea(QWidget* parent) : QScrollArea(parent) {
	// nothing
}

FixedScrollArea::~FixedScrollArea() {
	// nothing
}

QSize FixedScrollArea::sizeHint() const {
	int f = 2 * frameWidth();
	QSize sz(f, f);
	int h = fontMetrics().height();
	if(widget()) {
		sz += widget()->sizeHint();
	} else {
		sz += QSize(12 * h, 8 * h);
	}
	if(verticalScrollBarPolicy() == Qt::ScrollBarAlwaysOn)
		sz.setWidth(sz.width() + verticalScrollBar()->sizeHint().width());
	if(horizontalScrollBarPolicy() == Qt::ScrollBarAlwaysOn)
		sz.setHeight(sz.height() + horizontalScrollBar()->sizeHint().height());
	if(verticalScrollBarPolicy() != Qt::ScrollBarAlwaysOff)
		sz.setWidth(qMin(sz.width(), 36 * h));
	if(horizontalScrollBarPolicy() != Qt::ScrollBarAlwaysOff)
		sz.setHeight(qMin(sz.height(), 24 * h));
	return sz;
}

QSize FixedScrollArea::minimumSizeHint() const {
	int f = 2 * frameWidth();
	QSize sz(f, f);
	int h = fontMetrics().height();
	if(widget()) {
		sz += widget()->minimumSizeHint();
	} else {
		sz += QSize(12 * h, 8 * h);
	}
	if(verticalScrollBarPolicy() == Qt::ScrollBarAlwaysOn)
		sz.setWidth(sz.width() + verticalScrollBar()->sizeHint().width());
	if(horizontalScrollBarPolicy() == Qt::ScrollBarAlwaysOn)
		sz.setHeight(sz.height() + horizontalScrollBar()->sizeHint().height());
	if(verticalScrollBarPolicy() != Qt::ScrollBarAlwaysOff)
		sz.setWidth(qMin(sz.width(), 36 * h));
	if(horizontalScrollBarPolicy() != Qt::ScrollBarAlwaysOff)
		sz.setHeight(qMin(sz.height(), 24 * h));
	return sz;
}
