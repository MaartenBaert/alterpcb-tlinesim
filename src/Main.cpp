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

#include "ApplicationDirs.h"
#include "Basics.h"
#include "Icons.h"
#include "MainWindow.h"
#include "StringRegistry.h"
#include "TLineTypes.h"

int main(int argc, char* argv[]) {

	// create singletons
	StringRegistry string_registry;
	UNUSED(string_registry);

	// create application
	QApplication app(argc, argv);
	QCoreApplication::setOrganizationName("AlterPCB");
	QCoreApplication::setApplicationName("AlterPCB Transmission Line Simulator");
	QCoreApplication::setAttribute(Qt::AA_UseHighDpiPixmaps);

	InitApplicationDirs();
	LoadIcons();
	RegisterTLineTypes();

	MainWindow window;
	Q_UNUSED(window);
	int exit_code = app.exec();

	return exit_code;
}
