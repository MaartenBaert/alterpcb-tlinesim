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

#include "Icons.h"

#include "ApplicationDirs.h"

#include <iostream>

QIcon g_icon_simulation;

void LoadIcons() {

	if(g_application_data_dir.isEmpty()) {
		std::cerr << "Error: Could not load icons, data directory is missing." << std::endl;
		return;
	}
	QString icondir = g_application_data_dir + "/icons";

	g_icon_simulation = QIcon(icondir + "/icon-simulation.svg");

}
