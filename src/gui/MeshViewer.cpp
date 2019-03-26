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

#include "ColorMap.h"
#include "GenericMesh.h"

inline Box2D AdjustAspectRatio(const Box2D &view, real_t w, real_t h) {
	if((view.x2 - view.x1) * h > (view.y2 - view.y1) * w) {
		real_t center = (view.y1 + view.y2) * 0.5;
		real_t offset = (view.x2 - view.x1) * 0.5 * h / w;
		return Box2D(view.x1, view.x2, center - offset, center + offset);
	} else {
		real_t center = (view.x1 + view.x2) * 0.5;
		real_t offset = (view.y2 - view.y1) * 0.5 * w / h;
		return Box2D(center - offset, center + offset, view.y1, view.y2);
	}
}

MeshViewer::MeshViewer(QWidget* parent)
	: QWidget(parent) {

	m_zoom = 0.0;
	m_image_type = MESHIMAGETYPE_MESH;
	m_mesh_overlay = false;
	m_mode = 0;

	setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

}

MeshViewer::~MeshViewer() {
	// nothing
}

void MeshViewer::SetMesh(std::unique_ptr<GenericMesh> mesh) {
	m_mesh = std::move(mesh);
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

QSize MeshViewer::minimumSizeHint() const {
	return QSize(100, 100);
}

QSize MeshViewer::sizeHint() const {
	return QSize(500, 500);
}

void MeshViewer::paintEvent(QPaintEvent* event) {
	Q_UNUSED(event);
	QPainter painter(this);

	// get size
	qreal fw = (qreal) width() * devicePixelRatioF();
	qreal fh = (qreal) height() * devicePixelRatioF();
	int iw = rinti(ceil(fw)), ih = rinti(ceil(fh));
	if(iw <= 0 || ih <= 0)
		return;

	// background
	painter.fillRect(0, 0, width(), height(), QColor(64, 64, 64));

	if(m_mesh == NULL)
		return;

	// check whether we can get an image
	bool has_image = (m_image_type == MESHIMAGETYPE_MESH)? m_mesh->IsInitialized() : m_mesh->IsSolved();
	if(!has_image)
		return;

	// get view parameters
	Box2D world_box = m_mesh->GetWorldBox2D();
	Box2D world_focus = m_mesh->GetWorldFocus2D();
	real_t scale = exp10(-m_zoom);

	// calculate view
	Box2D view_box = AdjustAspectRatio(world_box, (real_t) fw, (real_t) fh);
	Box2D view_focus = AdjustAspectRatio(world_focus, (real_t) fw, (real_t) fh);
	Box2D view = {
		view_focus.x1 + (view_box.x1 - view_focus.x1) * scale,
		view_focus.x2 + (view_box.x2 - view_focus.x2) * scale,
		view_focus.y1 + (view_box.y1 - view_focus.y1) * scale,
		view_focus.y2 + (view_box.y2 - view_focus.y2) * scale,
	};
	std::swap(view.y1, view.y2); // Y-axis is upside down

	// get image data
	std::vector<real_t> image_value;
	std::vector<Vector2D> image_gradient;
	m_mesh->GetImage2D(image_value, image_gradient, (size_t) iw, (size_t) ih, view, m_image_type, m_mode);

	// create image
	QImage image(iw, ih, QImage::Format_RGB32);
	image.setDevicePixelRatio(devicePixelRatioF());

	// generate image
	if(m_image_type == MESHIMAGETYPE_MESH) {

		// color plot
		const ColorMap &cmap = COLORMAP_GRAYSCALE;
		for(size_t j = 0; j < (size_t) ih; ++j) {
			uint32_t *row = (uint32_t*) image.scanLine((int) j);
			real_t *row_value = image_value.data() + j * (size_t) iw;
			for(size_t i = 0; i < (size_t) iw; ++i) {
				row[i] = cmap(row_value[i]).ToUint32();
			}
		}

	} else if(m_image_type == MESHIMAGETYPE_EPOT || m_image_type == MESHIMAGETYPE_MPOT) {

		// get mesh data
		std::vector<real_t> image_mesh_value;
		std::vector<Vector2D> image_mesh_gradient;
		if(m_mesh_overlay) {
			m_mesh->GetImage2D(image_mesh_value, image_mesh_gradient, (size_t) iw, (size_t) ih, view, MESHIMAGETYPE_MESH, m_mode);
		}

		// color + contour plot
		real_t contours = 20.0, contour_scale = contours * (view.x2 - view.x1) / (real_t) fw;
		const ColorMap &cmap = COLORMAP_MAGMA;
		for(size_t j = 0; j < (size_t) ih; ++j) {
			uint32_t *row = (uint32_t*) image.scanLine((int) j);
			real_t *row_value = image_value.data() + j * (size_t) iw;
			Vector2D *row_gradient = image_gradient.data() + j * (size_t) iw;
			real_t *row_mesh_value = image_mesh_value.data() + j * (size_t) iw;
			for(size_t i = 0; i < (size_t) iw; ++i) {
				real_t value = row_value[i];
				real_t contour_range = hypot(row_gradient[i].x, row_gradient[i].y) * contour_scale;
				real_t temp = value * contours + 0.5;
				real_t temp2 = (temp - nearbyint(temp)) / contour_range;
				Color plot_color = cmap(fabs(value));
				if(m_mesh_overlay) {
					plot_color = ColorMix(plot_color, cmap(row_mesh_value[i]), 0.2f);
				}
				Color contour_color = {1.0f, 1.0f, 1.0f, 0.5f * fmaxf(0.0f, 1.0f - (float) fabs(temp2))};
				row[i] = ColorBlend(plot_color, contour_color).ToUint32();
			}
		}

	} else if(m_image_type == MESHIMAGETYPE_ENERGY) {

		// get mesh data
		std::vector<real_t> image_mesh_value;
		std::vector<Vector2D> image_mesh_gradient;
		if(m_mesh_overlay) {
			m_mesh->GetImage2D(image_mesh_value, image_mesh_gradient, (size_t) iw, (size_t) ih, view, MESHIMAGETYPE_MESH, m_mode);
		}

		// color + contour plot
		real_t log_scale = 1.0 / log(1e4);
		real_t contours = 20.0, contour_scale = contours * (view.x2 - view.x1) / (real_t) fw;
		const ColorMap &cmap = COLORMAP_MAGMA;
		for(size_t j = 0; j < (size_t) ih; ++j) {
			uint32_t *row = (uint32_t*) image.scanLine((int) j);
			real_t *row_value = image_value.data() + j * (size_t) iw;
			Vector2D *row_gradient = image_gradient.data() + j * (size_t) iw;
			real_t *row_mesh_value = image_mesh_value.data() + j * (size_t) iw;
			for(size_t i = 0; i < (size_t) iw; ++i) {
				real_t value = log(fabs(row_value[i])) * log_scale + 1.0;
				real_t contour_range = hypot(row_gradient[i].x, row_gradient[i].y) / row_value[i] * log_scale * contour_scale;
				real_t temp = value * contours + 0.5;
				real_t temp2 = (temp - nearbyint(temp)) / contour_range;
				Color plot_color = cmap(value);
				if(m_mesh_overlay) {
					plot_color = ColorMix(plot_color, cmap(row_mesh_value[i]), 0.2f);
				}
				Color contour_color = {1.0f, 1.0f, 1.0f, 0.5f * fmaxf(0.0f, 1.0f - (float) fabs(temp2))};
				row[i] = ColorBlend(plot_color, contour_color).ToUint32();
			}
		}

	} else {

		// get mesh data
		std::vector<real_t> image_mesh_value;
		std::vector<Vector2D> image_mesh_gradient;
		if(m_mesh_overlay) {
			m_mesh->GetImage2D(image_mesh_value, image_mesh_gradient, (size_t) iw, (size_t) ih, view, MESHIMAGETYPE_MESH, m_mode);
		}

		// color plot
		real_t log_scale = 1.0 / log(1e4);
		const ColorMap &cmap = COLORMAP_MAGMA;
		for(size_t j = 0; j < (size_t) ih; ++j) {
			uint32_t *row = (uint32_t*) image.scanLine((int) j);
			real_t *row_value = image_value.data() + j * (size_t) iw;
			real_t *row_mesh_value = image_mesh_value.data() + j * (size_t) iw;
			for(size_t i = 0; i < (size_t) iw; ++i) {
				real_t value = row_value[i];
				Color plot_color = cmap(log(fabs(value)) * log_scale + 1.0);
				if(m_mesh_overlay) {
					plot_color = ColorMix(plot_color, cmap(row_mesh_value[i]), 0.2f);
				}
				row[i] = plot_color.ToUint32();
			}
		}

	}

	// draw image
	qreal cut_x1 = fw * clamp<qreal>((world_box.x1 - view.x1) / (view.x2 - view.x1), 0.0, 1.0);
	qreal cut_x2 = fw * clamp<qreal>((world_box.x2 - view.x1) / (view.x2 - view.x1), 0.0, 1.0);
	qreal cut_y1 = fh * clamp<qreal>((world_box.y2 - view.y1) / (view.y2 - view.y1), 0.0, 1.0);
	qreal cut_y2 = fh * clamp<qreal>((world_box.y1 - view.y1) / (view.y2 - view.y1), 0.0, 1.0);
	painter.drawImage(QPointF(cut_x1, cut_y1), image, QRectF(cut_x1, cut_y1, cut_x2 - cut_x1, cut_y2 - cut_y1));

}
