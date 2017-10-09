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
#include "Qt.h"

#include <atomic>

class QProgressDialogThreaded : public QProgressDialog {
	Q_OBJECT

private:
	QTimer m_update_timer;
	std::atomic<int> m_task_progress;
	std::atomic<bool> m_task_canceled, m_task_stopped;

public:
	QProgressDialogThreaded(QWidget *parent = 0, Qt::WindowFlags f = 0);
	QProgressDialogThreaded(const QString &labelText, const QString &cancelButtonText, int minimum, int maximum, QWidget *parent = 0, Qt::WindowFlags f = 0);
	~QProgressDialogThreaded();

	void execThreaded(std::function<void(std::atomic<int>&, std::atomic<bool>&)> task);

private slots:
	void OnUpdateTimer();
	void OnCancel();

};
