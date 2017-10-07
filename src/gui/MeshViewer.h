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
#include "GenericMesh.h"
#include "Qt.h"

#include <memory>

class GenericMesh;

class MeshViewer : public QWidget {
	Q_OBJECT

private:
	std::unique_ptr<GenericMesh> m_mesh;
	real_t m_zoom;
	MeshImageType m_image_type;
	bool m_mesh_overlay;
	size_t m_mode;

public:
	MeshViewer(QWidget* parent);
	~MeshViewer();

	void SetMesh(std::unique_ptr<GenericMesh> mesh);
	void SetZoom(real_t zoom);
	void SetImageType(MeshImageType image_type);
	void SetMeshOverlay(bool mesh_overlay);
	void SetMode(size_t mode);

	virtual QSize minimumSizeHint() const override;
	virtual QSize sizeHint() const override;

public:
	inline GenericMesh* GetMesh() { return m_mesh.get(); }

protected:
	virtual void paintEvent(QPaintEvent* event) override;

};
