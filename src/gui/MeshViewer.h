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

#include <thread>

class GenericMesh;
class MeshViewer;

struct MeshRequestInfo {
	std::shared_ptr<GenericMesh> m_mesh;
	MeshImageType m_image_type;
	bool m_mesh_overlay, m_contour_lines;
	size_t m_mode;
	size_t m_image_w, m_image_h;
	Box2D m_image_view;
};

// object responsible for rendering, uses the worker thread event loop
class MeshRenderer : public QThread {
	Q_OBJECT

private:
	MeshViewer *m_meshviewer;

public:
	MeshRenderer(MeshViewer *meshviewer);

private:
	void GenerateImage(MeshRequestInfo &request_info, QImage &image);

public slots:
	void ProcessRequest();

signals:
	void CompletedRequest();

};

// worker thread object, uses the main event loop
class MeshWorker : public QThread {
	Q_OBJECT

private:
	MeshViewer *m_meshviewer;

public:
	MeshWorker(MeshViewer *meshviewer);

	virtual void run() override;

};

class MeshViewer : public QWidget {
	Q_OBJECT

	friend class MeshRenderer;

private:
	std::shared_ptr<GenericMesh> m_mesh;
	real_t m_zoom;
	MeshImageType m_image_type;
	bool m_mesh_overlay, m_contour_lines;
	size_t m_mode;

	std::mutex m_request_mutex;
	size_t m_request_id, m_response_id;
	MeshRequestInfo m_request_info, m_response_info;
	QImage m_response_image;

	MeshWorker *m_worker;

public:
	MeshViewer(QWidget* parent);
	~MeshViewer();

	void SetMesh(std::shared_ptr<GenericMesh> mesh);
	void SetZoom(real_t zoom);
	void SetImageType(MeshImageType image_type);
	void SetMeshOverlay(bool mesh_overlay);
	void SetContourLines(bool contour_lines);
	void SetMode(size_t mode);

public:
	virtual QSize minimumSizeHint() const override;
	virtual QSize sizeHint() const override;

protected:
	virtual void paintEvent(QPaintEvent *event) override;

public:
	inline GenericMesh* GetMesh() { return m_mesh.get(); }

signals:
	void NewRequest();

};
