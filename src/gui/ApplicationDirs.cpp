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

#include <iostream>

QString g_application_data_dir;

void InitApplicationDirs() {

	std::vector<QString> global_data_dirs = {
		QCoreApplication::applicationDirPath() + "/data",
		QCoreApplication::applicationDirPath() + "/../data",
	};
	for(auto &dir : global_data_dirs) {
		if(QDir(dir).exists()) {
			g_application_data_dir = dir;
			break;
		}
	}
	if(g_application_data_dir.isEmpty()) {
		std::cerr << "Warning: Could not find application data directory." << std::endl;
	}

}
