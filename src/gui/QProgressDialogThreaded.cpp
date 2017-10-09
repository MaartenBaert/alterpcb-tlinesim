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

#include "QProgressDialogThreaded.h"

#include <thread>

QProgressDialogThreaded::QProgressDialogThreaded(QWidget *parent, Qt::WindowFlags f)
	: QProgressDialog(parent, f) {
	connect(&m_update_timer, SIGNAL(timeout()), this, SLOT(OnUpdateTimer()));
	connect(this, SIGNAL(canceled()), this, SLOT(OnCancel()));
}

QProgressDialogThreaded::QProgressDialogThreaded(const QString &labelText, const QString &cancelButtonText, int minimum, int maximum, QWidget *parent, Qt::WindowFlags f)
	: QProgressDialog(labelText, cancelButtonText, minimum, maximum, parent, f) {
	connect(&m_update_timer, SIGNAL(timeout()), this, SLOT(OnUpdateTimer()));
	connect(this, SIGNAL(canceled()), this, SLOT(OnCancel()));
}

QProgressDialogThreaded::~QProgressDialogThreaded() {
	// nothing
}

void QProgressDialogThreaded::execThreaded(std::function<void(std::atomic<int>&, std::atomic<bool>&)> task) {

	// initialize state
	setWindowModality(Qt::WindowModal);
	setValue(0);
	m_task_progress = 0;
	m_task_canceled = false;
	m_task_stopped = false;

	// start task thread
	std::exception_ptr task_exception;
	std::thread task_thread([&]{

		// execute the provided task
		try {
			task(m_task_progress, m_task_canceled);
		} catch(...) {
			task_exception = std::current_exception();
		}

		// tell the main thread that we have stopped
		m_task_stopped = true;
		QMetaObject::invokeMethod(this, "OnUpdateTimer", Qt::QueuedConnection);

	});

	try {

		// process events while the task is running
		m_update_timer.start(50);
		while(!m_task_canceled && !m_task_stopped) {
			// There's a subtle unfixable bug here: QApplication::processEvents may continue to process several more
			// events after the dialog has been hidden by pressing the cancel() button, which means that the parent
			// window may start receiving events again before we have left this loop. As far as I can tell, all modal
			// Qt dialogs have this problem, because there is just no way to avoid it.
			QApplication::processEvents(QEventLoop::WaitForMoreEvents);
		}
		m_update_timer.stop();

	} catch(...) {

		// join with task thread and re-throw the exception
		m_task_canceled = true;
		task_thread.join();
		throw;

	}

	// join with task thread and propagate the exception if necessary
	task_thread.join();
	if(task_exception) {
		std::rethrow_exception(task_exception);
	}

}

void QProgressDialogThreaded::OnUpdateTimer() {
	if(!m_task_canceled) {
		setValue(m_task_progress);
	}
}

void QProgressDialogThreaded::OnCancel() {
	m_task_canceled = true;
}
