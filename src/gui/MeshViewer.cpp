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

#include "MeshViewer.h"

#include "Basics.h"
#include "ColorMap.h"
#include "GenericMesh.h"
#include "Vector.h"

#include <fstream> // TODO remove

inline Box2D AdjustAspectRatio(const Box2D &view, real_t w, real_t h) {
	if(view.Width() * h > view.Height() * w) {
		real_t center = (view.y1 + view.y2) * 0.5;
		real_t offset = (view.x2 - view.x1) * 0.5 * h / w;
		return Box2D(view.x1, view.x2, center - offset, center + offset);
	} else {
		real_t center = (view.x1 + view.x2) * 0.5;
		real_t offset = (view.y2 - view.y1) * 0.5 * w / h;
		return Box2D(center - offset, center + offset, view.y1, view.y2);
	}
}

inline Vector2D HaltonSequence(uint32_t index) {
	uint32_t i2 = index;
	i2 = ((i2 & 0xffff0000) >> 16) | ((i2 & 0x0000ffff) << 16);
	i2 = ((i2 & 0xff00ff00) >>  8) | ((i2 & 0x00ff00ff) <<  8);
	i2 = ((i2 & 0xf0f0f0f0) >>  4) | ((i2 & 0x0f0f0f0f) <<  4);
	i2 = ((i2 & 0xcccccccc) >>  2) | ((i2 & 0x33333333) <<  2);
	i2 = ((i2 & 0xaaaaaaaa) >>  1) | ((i2 & 0x55555555) <<  1);
	uint64_t i3 = index, i3a, i3b;
	i3a = (i3 * 2681326939) >> 44;
	i3b = i3 - i3a * 6561;
	i3 = i3a | (i3b << 32);
	i3a = ((i3 * 25891) >> 21) & 0x000007ff000007ff;
	i3b = i3 - i3a * 81;
	i3 = i3a | (i3b << 16);
	i3a = ((i3 * 57) >> 9) & 0x007f007f007f007f;
	i3b = i3 - i3a * 9;
	i3 = i3a | (i3b << 8);
	i3a = ((i3 * 11) >> 5) & 0x0707070707070707;
	i3b = i3 - i3a * 3;
	i3 = i3a + i3b * 3;
	i3 = (i3 & 0x00ff00ff00ff00ff) + ((i3 & 0xff00ff00ff00ff00) >>  8) * 9;
	i3 = (i3 & 0x0000ffff0000ffff) + ((i3 & 0xffff0000ffff0000) >> 16) * 81;
	i3 = (i3 & 0x00000000ffffffff) + ((i3 & 0xffffffff00000000) >> 32) * 6561;
	return Vector2D((real_t) i2 / 4294967296.0, (real_t) i3 / 43046721.0);
}

MeshRenderer::MeshRenderer(MeshViewer *meshviewer) {
	m_meshviewer = meshviewer;
	connect(meshviewer, SIGNAL(NewRequest()), this, SLOT(ProcessRequest()), Qt::QueuedConnection);
	connect(this, SIGNAL(CompletedRequest()), meshviewer, SLOT(update()), Qt::QueuedConnection);
}

void MeshRenderer::GenerateImage(MeshRequestInfo &request_info, QImage &image) {

	GenericMesh *mesh = request_info.m_mesh.get();
	auto image_type = request_info.m_image_type;
	auto mesh_overlay = request_info.m_mesh_overlay;
	auto contour_lines = request_info.m_contour_lines;
	auto mode = request_info.m_mode;
	auto image_w = request_info.m_image_w;
	auto image_h = request_info.m_image_h;
	auto image_view = request_info.m_image_view;
	std::swap(image_view.y1, image_view.y2); // Y-axis is upside down

	// check whether we have a mesh
	if(mesh == NULL)
		return;

	// check whether we can get an image
	bool has_image = (image_type == MESHIMAGETYPE_MESH)? mesh->IsInitialized() : mesh->IsSolved();
	if(!has_image)
		return;

	// get image data
	std::vector<real_t> image_value;
	mesh->GetImage2D(image_value, image_w, image_h, image_view, image_type, mode);

	// generate gradient
	std::vector<Vector2D> image_gradient(image_value.size());
	{
		uint32_t *row = (uint32_t*) image.scanLine(0);
		real_t *row1_value = image_value.data() + 0 * image_w;
		real_t *row2_value = image_value.data() + 1 * image_w;
		Vector2D *row_gradient = image_gradient.data() + 0 * image_w;
		row_gradient[0].x = row1_value[1] - row1_value[0];
		row_gradient[0].y = row2_value[0] - row1_value[0];
		for(size_t i = 1; i < image_w - 1; ++i) {
			row_gradient[i].x = 0.5 * (row1_value[i + 1] - row1_value[i - 1]);
			row_gradient[i].y = row2_value[i] - row1_value[i];
		}
		row_gradient[image_w - 1].x = row1_value[image_w - 1] - row1_value[image_w - 2];
		row_gradient[image_w - 1].y = row2_value[image_w - 1] - row1_value[image_w - 1];
	}
	for(size_t j = 1; j < image_h - 1; ++j) {
		uint32_t *row = (uint32_t*) image.scanLine((int) j);
		real_t *row0_value = image_value.data() + (j - 1) * image_w;
		real_t *row1_value = image_value.data() + j * image_w;
		real_t *row2_value = image_value.data() + (j + 1) * image_w;
		Vector2D *row_gradient = image_gradient.data() + j * image_w;
		row_gradient[0].x = row1_value[1] - row1_value[0];
		row_gradient[0].y = 0.5 * (row2_value[0] - row0_value[0]);
		for(size_t i = 1; i < image_w - 1; ++i) {
			row_gradient[i].x = 0.5 * (row1_value[i + 1] - row1_value[i - 1]);
			row_gradient[i].y = 0.5 * (row2_value[i] - row0_value[i]);
		}
		row_gradient[image_w - 1].x = row1_value[image_w - 1] - row1_value[image_w - 2];
		row_gradient[image_w - 1].y = 0.5 * (row2_value[image_w - 1] - row0_value[image_w - 1]);
	}
	{
		uint32_t *row = (uint32_t*) image.scanLine((int) image_h - 1);
		real_t *row0_value = image_value.data() + (image_h - 2) * image_w;
		real_t *row1_value = image_value.data() + (image_h - 1) * image_w;
		Vector2D *row_gradient = image_gradient.data() + (image_h - 1) * image_w;
		row_gradient[0].x = row1_value[1] - row1_value[0];
		row_gradient[0].y = row1_value[0] - row0_value[0];
		for(size_t i = 1; i < image_w - 1; ++i) {
			row_gradient[i].x = 0.5 * (row1_value[i + 1] - row1_value[i - 1]);
			row_gradient[i].y = row1_value[i] - row0_value[i];
		}
		row_gradient[image_w - 1].x = row1_value[image_w - 1] - row1_value[image_w - 2];
		row_gradient[image_w - 1].y = row1_value[image_w - 1] - row0_value[image_w - 1];
	}

	// create image
	image = QImage((int) image_w, (int) image_h, QImage::Format_RGB32);

	// generate image
	if(image_type == MESHIMAGETYPE_MESH) {

		// color plot
		const ColorMap &cmap = COLORMAP_GRAYSCALE;
		for(size_t j = 0; j < image_h; ++j) {
			uint32_t *row = (uint32_t*) image.scanLine((int) j);
			real_t *row_value = image_value.data() + j * image_w;
			for(size_t i = 0; i < image_w; ++i) {
				row[i] = cmap(row_value[i]).ToUint32();
			}
		}

	} else {

		bool logarithmic;
		real_t log_scale;
		if(image_type == MESHIMAGETYPE_EPOT || image_type == MESHIMAGETYPE_MPOT) {
			logarithmic = false;
			log_scale = 1.0;
		} else {
			logarithmic = true;
			log_scale = 1.0 / log((image_type == MESHIMAGETYPE_POYNTING)? 1e4 : 1e2);
		}

		// get mesh data
		std::vector<real_t> image_mesh_value;
		if(mesh_overlay) {
			mesh->GetImage2D(image_mesh_value, image_w, image_h, image_view, MESHIMAGETYPE_MESH, mode);
		}

		// field lines
		bool field_lines = (image_type == MESHIMAGETYPE_EFIELD || image_type == MESHIMAGETYPE_MFIELD);
		std::vector<real_t> image_field_value;
		if(field_lines) {
			image_field_value.resize((image_h + 2) * (image_w + 2), 0.0);
			size_t num_lines = image_w * image_h / 400;
			for(size_t k = 0; k < num_lines; ++k) {
				Vector2D halton = HaltonSequence((uint32_t) k + 1);
				real_t scalex = (image_view.x2 - image_view.x1) / (real_t) image_w;
				real_t scaley = (image_view.y2 - image_view.y1) / (real_t) image_h;
				for(size_t d = 0; d < 2; ++d) {
					real_t dir = (d == 0)? 1.0 : -1.0;
					real_t px = halton.x * (real_t) image_w, py = halton.y * (real_t) image_h;
					for(size_t l = 0; l < 15; ++l) {
						real_t x = image_view.x1 + px * scalex;
						real_t y = image_view.y1 + py * scaley;
						complex_t fx, fy;
						mesh->GetPointField(x, y, image_type, mode, fx, fy);
						real_t d = dir * sqrt(square(fx.real()) + square(fy.real()));
						if(d == 0.0)
							break;
						real_t px2 = px + fx.real() / d;
						real_t py2 = py - fy.real() / d; // Y-axis is upside down
						real_t pxx = (px + px2) * 0.5;
						real_t pyy = (py + py2) * 0.5;
						ptrdiff_t pxi = rintp(pxx);
						ptrdiff_t pyi = rintp(pyy);
						if(pxi < 0 || pxi > (ptrdiff_t) image_w || pyi < 0 || pyi > (ptrdiff_t) image_h)
							break;
						real_t pxf = pxx + 0.5 - (real_t) pxi;
						real_t pyf = pyy + 0.5 - (real_t) pyi;
						image_field_value[(size_t) (pyi + 0) * (image_w + 2) + (size_t) (pxi + 0)] += (1.0 - pyf) * (1.0 - pxf);
						image_field_value[(size_t) (pyi + 0) * (image_w + 2) + (size_t) (pxi + 1)] += (1.0 - pyf) * pxf;
						image_field_value[(size_t) (pyi + 1) * (image_w + 2) + (size_t) (pxi + 0)] += pyf * (1.0 - pxf);
						image_field_value[(size_t) (pyi + 1) * (image_w + 2) + (size_t) (pxi + 1)] += pyf * pxf;
						px = px2;
						py = py2;
					}
				}
			}
		}

		// color + contour plot
		real_t contours = 20.0;
		const ColorMap &cmap = COLORMAP_MAGMA;
		for(size_t j = 0; j < image_h; ++j) {
			uint32_t *row = (uint32_t*) image.scanLine((int) j);
			real_t *row_value = image_value.data() + j * image_w;
			Vector2D *row_gradient = image_gradient.data() + j * image_w;
			real_t *row_mesh_value = image_mesh_value.data() + j * image_w;
			real_t *row_field_value = image_field_value.data() + (j + 1) * (image_w + 2) + 1;
			for(size_t i = 0; i < image_w; ++i) {
				real_t value;
				if(logarithmic) {
					value = log(1e-8 + row_value[i]) * log_scale + 1.0;
				} else {
					value = fabs(row_value[i]);
				}
				Color plot_color = cmap(value);
				if(mesh_overlay) {
					plot_color = ColorMix(plot_color, cmap(row_mesh_value[i]), 0.2f);
				}
				if(contour_lines) {
					float alpha;
					if(field_lines) {
						alpha = 0.5f * std::min(1.0f, (float) row_field_value[i]);
					} else {
						real_t contour_range = hypot(row_gradient[i].x, row_gradient[i].y) * contours;
						if(logarithmic) {
							contour_range *= log_scale / (1e-8 + row_value[i]);
						}
						real_t temp = value * contours + 0.5;
						real_t temp2 = (temp - nearbyint(temp)) / contour_range;
						alpha = 0.5f * std::max(0.0f, 1.0f - (float) fabs(temp2));
					}
					Color line_color = {1.0f, 1.0f, 1.0f, alpha};
					plot_color = ColorBlend(plot_color, line_color);
				}
				row[i] = plot_color.ToUint32();
			}
		}

	}

}

void MeshRenderer::ProcessRequest() {

	// get new request
	bool new_request = false;
	size_t request_id;
	MeshRequestInfo request_info;
	{
		std::lock_guard<std::mutex> lock(m_meshviewer->m_request_mutex);
		Q_UNUSED(lock);
		if(m_meshviewer->m_request_id != m_meshviewer->m_response_id) {
			request_id = m_meshviewer->m_request_id;
			request_info = m_meshviewer->m_request_info;
			new_request = true;
		}
	}
	if(!new_request)
		return;

	// generate image
	QImage image;
	GenerateImage(request_info, image);

	// complete request
	bool completed_request = false;
	{
		std::lock_guard<std::mutex> lock(m_meshviewer->m_request_mutex);
		Q_UNUSED(lock);
		m_meshviewer->m_response_id = request_id;
		m_meshviewer->m_response_info = request_info;
		m_meshviewer->m_response_image = std::move(image);
		completed_request = true;
	}
	if(completed_request)
		emit CompletedRequest();

}

MeshWorker::MeshWorker(MeshViewer *meshviewer) : QThread(meshviewer) {
	m_meshviewer = meshviewer;
}

void MeshWorker::run() {
	MeshRenderer *renderer = new MeshRenderer(m_meshviewer);
	exec();
	renderer->deleteLater();
}

MeshViewer::MeshViewer(QWidget* parent)
	: QWidget(parent) {

	m_zoom = 0.0;
	m_image_type = MESHIMAGETYPE_MESH;
	m_mesh_overlay = false;
	m_contour_lines = false;
	m_mode = 0;

	m_request_id = 0;
	m_response_id = 0;
	m_request_info.m_mesh = nullptr;
	m_request_info.m_image_type = MESHIMAGETYPE_MESH;
	m_request_info.m_mesh_overlay = false;
	m_request_info.m_contour_lines = false;
	m_request_info.m_mode = 0;
	m_request_info.m_image_w = 0;
	m_request_info.m_image_h = 0;
	m_request_info.m_image_view = {0.0, 0.0, 0.0, 0.0};

	setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

	m_worker = new MeshWorker(this);
	m_worker->start();

}

MeshViewer::~MeshViewer() {
	m_worker->quit();
	m_worker->wait();
}

void MeshViewer::SetMesh(std::shared_ptr<GenericMesh> mesh) {
	m_mesh = mesh;
	update();
}

void MeshViewer::SetZoom(real_t zoom) {
	m_zoom = zoom;
	update();
}

void MeshViewer::SetImageType(MeshImageType image_type) {
	m_image_type = image_type;
	update();
}

void MeshViewer::SetMode(size_t mode) {
	m_mode = mode;
	update();
}

void MeshViewer::SetMeshOverlay(bool mesh_overlay) {
	m_mesh_overlay = mesh_overlay;
	update();
}

void MeshViewer::SetContourLines(bool contour_lines) {
	m_contour_lines = contour_lines;
	update();
}

QSize MeshViewer::minimumSizeHint() const {
	return QSize(100, 100);
}

QSize MeshViewer::sizeHint() const {
	return QSize(500, 500);
}

void MeshViewer::paintEvent(QPaintEvent *event) {
	Q_UNUSED(event);
	QPainter painter(this);

	// background
	painter.fillRect(0, 0, width(), height(), QColor(64, 64, 64));

	// get view parameters
	if(m_mesh == nullptr)
		return;
	Box2D world_box = m_mesh->GetWorldBox2D();
	Box2D world_focus = m_mesh->GetWorldFocus2D();
	real_t scale = exp10(-m_zoom);

	// calculate view
	Box2D view_world = AdjustAspectRatio(world_box, (real_t) width(), (real_t) height());
	Box2D view_focus = AdjustAspectRatio(world_focus, (real_t) width(), (real_t) height());
	Box2D view_full = {
		view_focus.x1 + (view_world.x1 - view_focus.x1) * scale,
		view_focus.x2 + (view_world.x2 - view_focus.x2) * scale,
		view_focus.y1 + (view_world.y1 - view_focus.y1) * scale,
		view_focus.y2 + (view_world.y2 - view_focus.y2) * scale,
	};
	if(view_full.Width() <= 0.0 || view_full.Height() <= 0.0)
		return;

	// calculate valid part of widget
	Box2D widget_full = {0.0, (real_t) width() * devicePixelRatioF(), 0.0, (real_t) height() * devicePixelRatioF()};
	if(widget_full.Width() < 2.0 || widget_full.Height() < 2.0)
		return;
	Box2D widget_valid = world_box.Map(view_full, widget_full).Clip(widget_full).NearbyInt();
	if(widget_valid.Width() < 2.0 || widget_valid.Height() < 2.0)
		return;

	// create image request
	MeshRequestInfo request_info;
	request_info.m_mesh = m_mesh;
	request_info.m_image_type = m_image_type;
	request_info.m_mesh_overlay = m_mesh_overlay;
	request_info.m_contour_lines = m_contour_lines;
	request_info.m_mode = m_mode;
	request_info.m_image_w = (size_t) widget_valid.Width();
	request_info.m_image_h = (size_t) widget_valid.Height();
	request_info.m_image_view = widget_valid.Map(widget_full, view_full);

	bool submit_request = true, emit_request = false;
	QImage image;
	Box2D image_view;
	{
		std::lock_guard<std::mutex> lock(m_request_mutex);
		Q_UNUSED(lock);

		// check whether we already have a usable response
		if(m_response_info.m_mesh == request_info.m_mesh &&
		   m_response_info.m_image_type == request_info.m_image_type &&
		   m_response_info.m_mesh_overlay == request_info.m_mesh_overlay &&
		   m_response_info.m_contour_lines == request_info.m_contour_lines &&
		   m_response_info.m_mode == request_info.m_mode) {
			image = m_response_image;
			image_view = m_response_info.m_image_view;
			submit_request = !(
				m_response_info.m_image_w == request_info.m_image_w &&
				m_response_info.m_image_h == request_info.m_image_h &&
				m_response_info.m_image_view == request_info.m_image_view
			);
		}

		// submit new request if required
		if(submit_request && !(
		   m_request_info.m_mesh == request_info.m_mesh &&
		   m_request_info.m_image_type == request_info.m_image_type &&
		   m_request_info.m_mesh_overlay == request_info.m_mesh_overlay &&
		   m_request_info.m_contour_lines == request_info.m_contour_lines &&
		   m_request_info.m_mode == request_info.m_mode &&
		   m_request_info.m_image_w == request_info.m_image_w &&
		   m_request_info.m_image_h == request_info.m_image_h &&
		   m_request_info.m_image_view == request_info.m_image_view)) {
			++m_request_id;
			m_request_info = request_info;
			emit_request = true;
		}

	}
	if(emit_request) {
		emit NewRequest();
	}

	// draw image
	if(!image.isNull()) {

		// draw image
		image.setDevicePixelRatio(devicePixelRatioF());
		Box2D target = image_view.Map(view_full, Box2D(0.0, (real_t) width(), (real_t) height(), 0.0));
		painter.setRenderHint(QPainter::SmoothPixmapTransform);
		painter.drawImage(QRectF(target.x1, target.y2, target.Width(), -target.Height()), image);

	}

	// draw loading test
	if(submit_request) {
		painter.setFont(QFont("Sans", 32, QFont::Bold));
		painter.setPen(QColor(255, 255, 255, 128));
		painter.drawText(0, 0, width(), height(), Qt::AlignHCenter | Qt::AlignVCenter, "Loading ...");
	}

}
