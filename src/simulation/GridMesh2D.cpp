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

#include "GridMesh2D.h"

#include "FemMatrix.h"
#include "MatrixMarket.h"
#include "MiscMath.h"
#include "StringHelper.h"

#include <cfloat>

#include <algorithm>
#include <chrono>
#include <deque>
#include <iostream>

#define SIMULATION_VERBOSE 1
#define SIMULATION_SAVE_MATRIXMARKET 0

complex_t Interp_Ax(complex_t e1x0, complex_t e1x1, complex_t e2x0, complex_t e2x1, real_t fx, real_t fy) {
	return lerp(e1x0, e1x1, fy) + lerp(e2x0, e2x1, fy) * (fx - 0.5);
}
complex_t Interp_Ax_Dy(complex_t e1x0, complex_t e1x1, complex_t e2x0, complex_t e2x1, real_t fx, real_t fy) {
	UNUSED(fy); return (e1x1 - e1x0) + (e2x1 - e2x0) * (fx - 0.5);
}

complex_t Interp_Ay(complex_t e1y0, complex_t e1y1, complex_t e2y0, complex_t e2y1, real_t fx, real_t fy) {
	return lerp(e1y0, e1y1, fx) + lerp(e2y0, e2y1, fx) * (fy - 0.5);
}
complex_t Interp_Ay_Dx(complex_t e1y0, complex_t e1y1, complex_t e2y0, complex_t e2y1, real_t fx, real_t fy) {
	UNUSED(fx); return (e1y1 - e1y0) + (e2y1 - e2y0) * (fy - 0.5);
}

complex_t Interp_Az(complex_t n00, complex_t n01, complex_t n10, complex_t n11, complex_t ex0, complex_t ex1, complex_t ey0, complex_t ey1, real_t fx, real_t fy) {
	return lerp(lerp(n00, n01, fx), lerp(n10, n11, fx), fy) + lerp(ex0, ex1, fy) * fx * (1.0 - fx) + lerp(ey0, ey1, fx) * fy * (1.0 - fy);
}
complex_t Interp_Az_Dx(complex_t n00, complex_t n01, complex_t n10, complex_t n11, complex_t ex0, complex_t ex1, complex_t ey0, complex_t ey1, real_t fx, real_t fy) {
	return lerp(n01 - n00, n11 - n10, fy) + lerp(ex0, ex1, fy) * (1.0 - 2.0 * fx) + (ey1 - ey0) * fy * (1.0 - fy);
}
complex_t Interp_Az_Dy(complex_t n00, complex_t n01, complex_t n10, complex_t n11, complex_t ex0, complex_t ex1, complex_t ey0, complex_t ey1, real_t fx, real_t fy) {
	return lerp(n10 - n00, n11 - n01, fx) + (ex1 - ex0) * fx * (1.0 - fx) + lerp(ey0, ey1, fx) * (1.0 - 2.0 * fy);
}

// Returns the maximum length of a vector of phasors.
real_t PhasorVectorMax(complex_t x, complex_t y, complex_t z) {
	real_t a = x.real() * x.real() + y.real() * y.real() + z.real() * z.real();
	real_t b = x.real() * x.imag() + y.real() * y.imag() + z.real() * z.imag();
	real_t c = x.imag() * x.imag() + y.imag() * y.imag() + z.imag() * z.imag();
	return sqrt(hypot(0.5 * (a - c), b) + 0.5 * (a + c));
}

template<class EigenSparseMatrix>
void EigenSparseFree(EigenSparseMatrix &matrix) {
	matrix.resize(0, 0);
	matrix.data().squeeze();
}

GridMesh2D::GridMesh2D(SolverType solver_type, ElementType element_type, const Box2D &world_box, const Box2D &world_focus, real_t max_frequency, real_t lambda_factor, real_t grid_inc, real_t grid_epsilon) {
	if(!FinitePositive(grid_inc))
		throw std::runtime_error("GridMesh2D error: grid_inc must be positive.");
	if(!FinitePositive(grid_epsilon))
		throw std::runtime_error("GridMesh2D error: grid_epsilon must be positive.");
	m_solver_type = solver_type;
	m_element_type = element_type;
	m_world_box = world_box.Normalize();
	m_world_focus = world_focus.Normalize();
	m_max_frequency = max_frequency;
	m_lambda_factor = lambda_factor;
	m_grid_inc = grid_inc;
	m_grid_epsilon = grid_epsilon;
	m_pml_box = m_world_box;
	m_pml_step = Box2D(REAL_MAX, REAL_MAX, REAL_MAX, REAL_MAX);
	m_pml_attenuation = 0.0;
	m_vars_static_e_free = 0;
	m_vars_static_e_fixed = 0;
	m_vars_static_m_free = 0;
	m_vars_static_m_fixed = 0;
	m_vars_full_em = 0;
}

GridMesh2D::~GridMesh2D() {
	// nothing
}

void GridMesh2D::SetPML(const Box2D &box, real_t step, real_t attenuation) {
	SetPML(box, Box2D(step, step, step, step), attenuation);
}

void GridMesh2D::SetPML(const Box2D &box, const Box2D &step, real_t attenuation) {
	m_pml_box = box.Clip(m_world_box).Normalize();
	m_pml_step = step;
	m_pml_attenuation = attenuation;
}

size_t GridMesh2D::AddPort(const Vector2D &anchor, bool infinite_area) {
	if(IsInitialized())
		throw std::runtime_error("GridMesh2D error: Can't add port after initialization.");
	m_ports.emplace_back(anchor, infinite_area);
	return m_ports.size() - 1;
}

void GridMesh2D::AddConductor(const Box2D &box, real_t step, const MaterialConductor *material, size_t port) {
	AddConductor(box, Box2D(step, step, step, step), material, port);
}

void GridMesh2D::AddConductor(const Box2D &box, const Box2D &step, const MaterialConductor *material, size_t port) {
	if(IsInitialized())
		throw std::runtime_error("GridMesh2D error: Can't add conductor after initialization.");
	if(!std::isfinite(box.x1) || !std::isfinite(box.x2) || !std::isfinite(box.y1) || !std::isfinite(box.y2))
		throw std::runtime_error("GridMesh2D error: Conductor box must be finite.");
	if(!FinitePositive(step.x1) || !FinitePositive(step.x2) || !FinitePositive(step.y1) || !FinitePositive(step.y2))
		throw std::runtime_error("GridMesh2D error: Conductor step must be positive.");
	if(material == NULL)
		throw std::runtime_error("GridMesh2D error: Material can't be NULL.");
	if(port >= m_ports.size())
		throw std::runtime_error("GridMesh2D error: Invalid port index.");
	m_conductors.emplace_back(box.Clip(m_world_box).Normalize(), step, material, port);
}

void GridMesh2D::AddDielectric(const Box2D &box, real_t step, const MaterialDielectric *material) {
	AddDielectric(box, Box2D(step, step, step, step), material);
}

void GridMesh2D::AddDielectric(const Box2D &box, const Box2D &step, const MaterialDielectric *material) {
	if(IsInitialized())
		throw std::runtime_error("GridMesh2D error: Can't add dielectric after initialization.");
	if(!std::isfinite(box.x1) || !std::isfinite(box.x2) || !std::isfinite(box.y1) || !std::isfinite(box.y2))
		throw std::runtime_error("GridMesh2D error: Dielectric box must be finite.");
	if(!FinitePositive(step.x1) || !FinitePositive(step.x2) || !FinitePositive(step.y1) || !FinitePositive(step.y2))
		throw std::runtime_error("GridMesh2D error: Dielectric step must be positive.");
	if(material == NULL)
		throw std::runtime_error("GridMesh2D error: Material can't be NULL.");
	m_dielectrics.emplace_back(box.Clip(m_world_box).Normalize(), step, material);
}

void GridMesh2D::AddIntegrationLine(const Box2D &box) {
	if(IsInitialized())
		throw std::runtime_error("GridMesh2D error: Can't add integration line after initialization.");
	if(!std::isfinite(box.x1) || !std::isfinite(box.x2) || !std::isfinite(box.y1) || !std::isfinite(box.y2))
		throw std::runtime_error("GridMesh2D error: Integration line must be finite.");
	if(fabs(box.x1 - box.x2) >= m_grid_epsilon && fabs(box.y1 - box.y2) >= m_grid_epsilon)
		throw std::runtime_error("GridMesh2D error: Integration line must be horizontal or vertical.");
	m_integration_lines.emplace_back(box.Clip(m_world_box).Normalize());
}

Box2D GridMesh2D::GetWorldBox2D() {
	return m_world_box;
}

Box2D GridMesh2D::GetWorldFocus2D() {
	return m_world_focus;
}

void GridMesh2D::GetImage2D(std::vector<real_t> &image_value, std::vector<Vector2D> &image_gradient, size_t width, size_t height, const Box2D &view, MeshImageType type, size_t mode) {
	if(!IsInitialized())
		throw std::runtime_error("GridMesh2D error: The mesh must be initialized first.");
	if(type != MESHIMAGETYPE_MESH) {
		if(!IsSolved())
			throw std::runtime_error("GridMesh2D error: The mesh must be solved first.");
		if(mode >= GetModeCount())
			throw std::runtime_error("GridMesh2D error: Invalid mode index.");
	}

	// clear image data
	image_value.clear();
	image_value.resize(width * height);
	if(type != MESHIMAGETYPE_MESH) {
		image_gradient.clear();
		image_gradient.resize(width * height);
	}

	// prepare grid
	std::vector<size_t> index_x, index_y;
	std::vector<real_t> frac_x, frac_y;
	if(type == MESHIMAGETYPE_MESH) {
		PrepareCellImage(index_x, frac_x, m_grid_x, m_midpoints_x, view.x1, view.x2, width);
		PrepareCellImage(index_y, frac_y, m_grid_y, m_midpoints_y, view.y1, view.y2, height);
	} else {
		PrepareNodeImage(index_x, frac_x, m_grid_x, view.x1, view.x2, width);
		PrepareNodeImage(index_y, frac_y, m_grid_y, view.y1, view.y2, height);
	}

	// generate data
	switch(type) {
		case MESHIMAGETYPE_MESH: {
			std::vector<real_t> cell_values;
			GetCellValues(cell_values, mode, type);
			for(size_t j = 0; j < height; ++j) {
				real_t *row_value = image_value.data() + j * width;
				for(size_t i = 0; i < width; ++i) {
					size_t ix = index_x[i], iy = index_y[j];
					real_t v00 = cell_values[GetCellIndex(ix    , iy    )];
					real_t v01 = cell_values[GetCellIndex(ix + 1, iy    )];
					real_t v10 = cell_values[GetCellIndex(ix    , iy + 1)];
					real_t v11 = cell_values[GetCellIndex(ix + 1, iy + 1)];
					real_t ex0 = (GetEdgeX(ix    , iy + 1).m_var_full_mt1 == INDEX_NONE)? std::max(fabs(frac_y[j]) - 0.5, 0.0) : 1.0;
					real_t ex1 = (GetEdgeX(ix + 1, iy + 1).m_var_full_mt1 == INDEX_NONE)? std::max(fabs(frac_y[j]) - 0.5, 0.0) : 1.0;
					real_t ey0 = (GetEdgeY(ix + 1, iy    ).m_var_full_mt1 == INDEX_NONE)? std::max(fabs(frac_x[i]) - 0.5, 0.0) : 1.0;
					real_t ey1 = (GetEdgeY(ix + 1, iy + 1).m_var_full_mt1 == INDEX_NONE)? std::max(fabs(frac_x[i]) - 0.5, 0.0) : 1.0;
					real_t fx = clamp(frac_x[i] + 0.5, 0.0, 1.0);
					real_t fy = clamp(frac_y[j] + 0.5, 0.0, 1.0);
					real_t v0 = v00 + (v01 - v00) * fx;
					real_t v1 = v10 + (v11 - v10) * fx;
					real_t ex = ex0 + (ex1 - ex0) * fx;
					real_t ey = ey0 + (ey1 - ey0) * fy;
					row_value[i] = (v0 + (v1 - v0) * fy) * ex * ey;
				}
			}
			break;
		}
		case MESHIMAGETYPE_EPOT: {
			SolutionField &field = m_solution_fields[mode];
			real_t scale = 0.0;
			for(size_t i = 0; i < m_nodes.size(); ++i) {
				scale = std::max(scale, std::abs(field.m_field_nodes[i].e1.real()));
			}
			std::cerr << "EPOT scale = " << scale << std::endl;
			scale = 1.0 / scale;
			for(size_t j = 0; j < height; ++j) {
				real_t *row_value = image_value.data() + j * width;
				Vector2D *row_gradient = image_gradient.data() + j * width;
				for(size_t i = 0; i < width; ++i) {
					size_t ix = index_x[i], iy = index_y[j];
					real_t fx = frac_x[i], fy = frac_y[j];
					complex_t en00 = field.m_field_nodes[GetNodeIndex(ix    , iy    )].e1;
					complex_t en01 = field.m_field_nodes[GetNodeIndex(ix + 1, iy    )].e1;
					complex_t en10 = field.m_field_nodes[GetNodeIndex(ix    , iy + 1)].e1;
					complex_t en11 = field.m_field_nodes[GetNodeIndex(ix + 1, iy + 1)].e1;
					complex_t eex0 = field.m_field_edges_x[GetEdgeXIndex(ix    , iy    )].e2;
					complex_t eex1 = field.m_field_edges_x[GetEdgeXIndex(ix    , iy + 1)].e2;
					complex_t eey0 = field.m_field_edges_y[GetEdgeYIndex(ix    , iy    )].e2;
					complex_t eey1 = field.m_field_edges_y[GetEdgeYIndex(ix + 1, iy    )].e2;
					row_value[i] = Interp_Az(en00, en01, en10, en11, eex0, eex1, eey0, eey1, fx, fy).real() * scale;
					row_gradient[i].x = Interp_Az_Dx(en00, en01, en10, en11, eex0, eex1, eey0, eey1, fx, fy).real() / (m_grid_x[ix + 1] - m_grid_x[ix]) * scale;
					row_gradient[i].y = Interp_Az_Dy(en00, en01, en10, en11, eex0, eex1, eey0, eey1, fx, fy).real() / (m_grid_y[iy + 1] - m_grid_y[iy]) * scale;
				}
			}
			break;
		}
		case MESHIMAGETYPE_MPOT: {
			SolutionField &field = m_solution_fields[mode];
			real_t scale = 0.0;
			for(size_t i = 0; i < m_nodes.size(); ++i) {
				scale = std::max(scale, std::abs(field.m_field_nodes[i].m1.real()));
			}
			std::cerr << "MPOT scale = " << scale << std::endl;
			scale = 1.0 / scale;
			for(size_t j = 0; j < height; ++j) {
				real_t *row_value = image_value.data() + j * width;
				Vector2D *row_gradient = image_gradient.data() + j * width;
				for(size_t i = 0; i < width; ++i) {
					size_t ix = index_x[i], iy = index_y[j];
					real_t fx = frac_x[i], fy = frac_y[j];
					complex_t mn00 = field.m_field_nodes[GetNodeIndex(ix    , iy    )].m1;
					complex_t mn01 = field.m_field_nodes[GetNodeIndex(ix + 1, iy    )].m1;
					complex_t mn10 = field.m_field_nodes[GetNodeIndex(ix    , iy + 1)].m1;
					complex_t mn11 = field.m_field_nodes[GetNodeIndex(ix + 1, iy + 1)].m1;
					complex_t mex0 = field.m_field_edges_x[GetEdgeXIndex(ix    , iy    )].m2;
					complex_t mex1 = field.m_field_edges_x[GetEdgeXIndex(ix    , iy + 1)].m2;
					complex_t mey0 = field.m_field_edges_y[GetEdgeYIndex(ix    , iy    )].m2;
					complex_t mey1 = field.m_field_edges_y[GetEdgeYIndex(ix + 1, iy    )].m2;
					row_value[i] = Interp_Az(mn00, mn01, mn10, mn11, mex0, mex1, mey0, mey1, fx, fy).real() * scale;
					row_gradient[i].x = Interp_Az_Dx(mn00, mn01, mn10, mn11, mex0, mex1, mey0, mey1, fx, fy).real() / (m_grid_x[ix + 1] - m_grid_x[ix]) * scale;
					row_gradient[i].y = Interp_Az_Dy(mn00, mn01, mn10, mn11, mex0, mex1, mey0, mey1, fx, fy).real() / (m_grid_y[iy + 1] - m_grid_y[iy]) * scale;
				}
			}
			break;
		}
		case MESHIMAGETYPE_EFIELD:
		case MESHIMAGETYPE_EFIELD_X:
		case MESHIMAGETYPE_EFIELD_Y:
		case MESHIMAGETYPE_EFIELD_Z: {
			real_t scale = 0.0;
			size_t best_ix = 0, best_iy = 0;
			for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
				for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
					real_t scale_prev = scale;
					complex_t fex, fey, fez, fmx, fmy, fmz;
					GetCellField(mode, ix, iy, 0.0, 0.0, fex, fey, fez, fmx, fmy, fmz);
					scale = std::max(scale, PhasorVectorMax(fex, fey, fez));
					GetCellField(mode, ix, iy, 0.0, 1.0, fex, fey, fez, fmx, fmy, fmz);
					scale = std::max(scale, PhasorVectorMax(fex, fey, fez));
					GetCellField(mode, ix, iy, 1.0, 0.0, fex, fey, fez, fmx, fmy, fmz);
					scale = std::max(scale, PhasorVectorMax(fex, fey, fez));
					GetCellField(mode, ix, iy, 1.0, 1.0, fex, fey, fez, fmx, fmy, fmz);
					scale = std::max(scale, PhasorVectorMax(fex, fey, fez));
					if(scale > scale_prev) {
						best_ix = ix;
						best_iy = iy;
					}
				}
			}
			std::cerr << "EFIELD scale = " << scale << " (" << best_ix << "," << best_iy << ")" << std::endl;
			scale = 1.0 / scale;
			switch(type) {
				case MESHIMAGETYPE_EFIELD: {
					for(size_t j = 0; j < height; ++j) {
						real_t *row_value = image_value.data() + j * width;
						for(size_t i = 0; i < width; ++i) {
							complex_t fex, fey, fez, fmx, fmy, fmz;
							GetCellField(mode, index_x[i], index_y[j], frac_x[i], frac_y[j], fex, fey, fez, fmx, fmy, fmz);
							row_value[i] = sqrt(square(fex.real()) + square(fey.real()) + square(fez.real())) * scale;
						}
					}
					break;
				}
				case MESHIMAGETYPE_EFIELD_X: {
					for(size_t j = 0; j < height; ++j) {
						real_t *row_value = image_value.data() + j * width;
						for(size_t i = 0; i < width; ++i) {
							complex_t fex, fey, fez, fmx, fmy, fmz;
							GetCellField(mode, index_x[i], index_y[j], frac_x[i], frac_y[j], fex, fey, fez, fmx, fmy, fmz);
							row_value[i] = fex.real() * scale;
						}
					}
					break;
				}
				case MESHIMAGETYPE_EFIELD_Y: {
					for(size_t j = 0; j < height; ++j) {
						real_t *row_value = image_value.data() + j * width;
						for(size_t i = 0; i < width; ++i) {
							complex_t fex, fey, fez, fmx, fmy, fmz;
							GetCellField(mode, index_x[i], index_y[j], frac_x[i], frac_y[j], fex, fey, fez, fmx, fmy, fmz);
							row_value[i] = fey.real() * scale;
						}
					}
					break;
				}
				case MESHIMAGETYPE_EFIELD_Z: {
					for(size_t j = 0; j < height; ++j) {
						real_t *row_value = image_value.data() + j * width;
						for(size_t i = 0; i < width; ++i) {
							complex_t fex, fey, fez, fmx, fmy, fmz;
							GetCellField(mode, index_x[i], index_y[j], frac_x[i], frac_y[j], fex, fey, fez, fmx, fmy, fmz);
							row_value[i] = fez.real() * scale;
						}
					}
					break;
				}
				default: assert(false); break;
			}
			break;
		}
		case MESHIMAGETYPE_MFIELD:
		case MESHIMAGETYPE_MFIELD_X:
		case MESHIMAGETYPE_MFIELD_Y:
		case MESHIMAGETYPE_MFIELD_Z: {
			real_t scale = 0.0;
			for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
				for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
					complex_t fex, fey, fez, fmx, fmy, fmz;
					GetCellField(mode, ix, iy, 0.0, 0.0, fex, fey, fez, fmx, fmy, fmz);
					scale = std::max(scale, PhasorVectorMax(fmx, fmy, fmz));
					GetCellField(mode, ix, iy, 0.0, 1.0, fex, fey, fez, fmx, fmy, fmz);
					scale = std::max(scale, PhasorVectorMax(fmx, fmy, fmz));
					GetCellField(mode, ix, iy, 1.0, 0.0, fex, fey, fez, fmx, fmy, fmz);
					scale = std::max(scale, PhasorVectorMax(fmx, fmy, fmz));
					GetCellField(mode, ix, iy, 1.0, 1.0, fex, fey, fez, fmx, fmy, fmz);
					scale = std::max(scale, PhasorVectorMax(fmx, fmy, fmz));
				}
			}
			std::cerr << "MFIELD scale = " << scale << std::endl;
			scale = 1.0 / scale;
			switch(type) {
				case MESHIMAGETYPE_MFIELD: {
					for(size_t j = 0; j < height; ++j) {
						real_t *row_value = image_value.data() + j * width;
						for(size_t i = 0; i < width; ++i) {
							complex_t fex, fey, fez, fmx, fmy, fmz;
							GetCellField(mode, index_x[i], index_y[j], frac_x[i], frac_y[j], fex, fey, fez, fmx, fmy, fmz);
							row_value[i] = sqrt(square(fmx.real()) + square(fmy.real()) + square(fmz.real())) * scale;
						}
					}
					break;
				}
				case MESHIMAGETYPE_MFIELD_X: {
					for(size_t j = 0; j < height; ++j) {
						real_t *row_value = image_value.data() + j * width;
						for(size_t i = 0; i < width; ++i) {
							complex_t fex, fey, fez, fmx, fmy, fmz;
							GetCellField(mode, index_x[i], index_y[j], frac_x[i], frac_y[j], fex, fey, fez, fmx, fmy, fmz);
							row_value[i] = fmx.real() * scale;
						}
					}
					break;
				}
				case MESHIMAGETYPE_MFIELD_Y: {
					for(size_t j = 0; j < height; ++j) {
						real_t *row_value = image_value.data() + j * width;
						for(size_t i = 0; i < width; ++i) {
							complex_t fex, fey, fez, fmx, fmy, fmz;
							GetCellField(mode, index_x[i], index_y[j], frac_x[i], frac_y[j], fex, fey, fez, fmx, fmy, fmz);
							row_value[i] = fmy.real() * scale;
						}
					}
					break;
				}
				case MESHIMAGETYPE_MFIELD_Z: {
					for(size_t j = 0; j < height; ++j) {
						real_t *row_value = image_value.data() + j * width;
						for(size_t i = 0; i < width; ++i) {
							complex_t fex, fey, fez, fmx, fmy, fmz;
							GetCellField(mode, index_x[i], index_y[j], frac_x[i], frac_y[j], fex, fey, fez, fmx, fmy, fmz);
							row_value[i] = fmz.real() * scale;
						}
					}
					break;
				}
				default: assert(false); break;
			}
			break;
		}
		case MESHIMAGETYPE_ENERGY:
		case MESHIMAGETYPE_CURRENT: {
			std::vector<std::array<real_t, 4>> cellnode_values;
			GetCellNodeValues(cellnode_values, mode, type);
			for(size_t j = 0; j < height; ++j) {
				real_t *row_value = image_value.data() + j * width;
				Vector2D *row_gradient = image_gradient.data() + j * width;
				for(size_t i = 0; i < width; ++i) {
					size_t ix = index_x[i], iy = index_y[j];
					size_t cell_index = GetCellIndex(ix, iy);
					real_t v00 = cellnode_values[cell_index][0];
					real_t v01 = cellnode_values[cell_index][1];
					real_t v10 = cellnode_values[cell_index][2];
					real_t v11 = cellnode_values[cell_index][3];
					real_t v0 = v00 + (v01 - v00) * frac_x[i];
					real_t v1 = v10 + (v11 - v10) * frac_x[i];
					row_value[i] = v0 + (v1 - v0) * frac_y[j];
					real_t gx0 = v01 - v00, gx1 = v11 - v10;
					real_t gy0 = v10 - v00, gy1 = v11 - v01;
					row_gradient[i].x = (gx0 + (gx1 - gx0) * frac_y[j]) / (m_grid_x[ix + 1] - m_grid_x[ix]);
					row_gradient[i].y = (gy0 + (gy1 - gy0) * frac_x[i]) / (m_grid_y[iy + 1] - m_grid_y[iy]);
				}
			}
			break;
		}
	}

}

void GridMesh2D::DoInitialize() {
	assert(!IsInitialized());

#if SIMULATION_VERBOSE
	auto t1 = std::chrono::high_resolution_clock::now();
#endif
	InitGrid();
#if SIMULATION_VERBOSE
	auto t2 = std::chrono::high_resolution_clock::now();
#endif
	InitCells();
#if SIMULATION_VERBOSE
	auto t3 = std::chrono::high_resolution_clock::now();
#endif
	InitVariables();
#if SIMULATION_VERBOSE
	auto t4 = std::chrono::high_resolution_clock::now();
#endif

#if SIMULATION_VERBOSE
	std::cerr << "Mesh stats:"
			  << " grid=" << m_grid_x.size() << "x" << m_grid_y.size()
			  << " nodes=" << m_nodes.size()
			  << " cells=" << m_cells.size()
			  << std::endl;
	std::cerr << "Static E vars:"
			  << " free=" << m_vars_static_e_free
			  << " fixed=" << m_vars_static_e_fixed
			  << std::endl;
	std::cerr << "Static M vars:"
			  << " free=" << m_vars_static_m_free
			  << " fixed=" << m_vars_static_m_fixed
			  << std::endl;
	std::cerr << "Full EM vars:"
			  << " free=" << m_vars_full_em
			  << std::endl;
	std::cerr << "Init time:"
			  << " grid=" << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us"
			  << " cells=" << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << "us"
			  << " vars=" << std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() << "us"
			  << std::endl;
#endif

}

void GridMesh2D::DoSolve() {
	assert(IsInitialized());

#if SIMULATION_VERBOSE
	auto t1 = std::chrono::high_resolution_clock::now();
#endif
	BuildMatrices();
#if SIMULATION_VERBOSE
	auto t2 = std::chrono::high_resolution_clock::now();
#endif
	SolveStaticModes();
#if SIMULATION_VERBOSE
	auto t3 = std::chrono::high_resolution_clock::now();
#endif
	SolveStaticEigenModes();
#if SIMULATION_VERBOSE
	auto t4 = std::chrono::high_resolution_clock::now();
#endif
	SolveFullEigenModes();
#if SIMULATION_VERBOSE
	auto t5 = std::chrono::high_resolution_clock::now();
#endif

#if SIMULATION_VERBOSE
	std::cerr << "GridMesh2D solve time:"
			  << " build=" << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us"
			  << " static=" << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << "us"
			  << " static-eigen=" << std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() << "us"
			  << " full-eigen=" << std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count() << "us"
			  << std::endl;
#endif

	size_t mode = 0;
	size_t cx = m_grid_x.size() - 1, cy = m_grid_y.size() - 1;
	std::vector<complex_t> fex(4 * cx * cy), fey(4 * cx * cy), fez(4 * cx * cy), fmx(4 * cx * cy), fmy(4 * cx * cy), fmz(4 * cx * cy);
	for(size_t iy = 0; iy < cy; ++iy) {
		for(size_t ix = 0; ix < cx; ++ix) {
			for(size_t fy = 0; fy < 2; ++fy) {
				for(size_t fx = 0; fx < 2; ++fx) {
					size_t ii = cx * (2 * iy + fy) + (2 * ix + fx);
					GetCellField(mode, ix, iy, (real_t) fx, (real_t) fy, fex[ii], fey[ii], fez[ii], fmx[ii], fmy[ii], fmz[ii]);
				}
			}
		}
	}
	auto savefield = [](const std::string &filename, const std::vector<complex_t> &field) {
		std::filebuf stream;
		if(stream.open(filename, std::ios_base::out | std::ios_base::binary) == NULL)
			throw std::runtime_error("Could not open file '" + filename + "' for writing.");
		stream.sputn((const char*) field.data(), field.size() * sizeof(complex_t));
	};
	savefield("field_fex.dat", fex);
	savefield("field_fey.dat", fey);
	savefield("field_fez.dat", fez);
	savefield("field_fmx.dat", fmx);
	savefield("field_fmy.dat", fmy);
	savefield("field_fmz.dat", fmz);

}

void GridMesh2D::DoCleanup() {
	EigenSparseFree(m_matrix_static_e[0]);
	EigenSparseFree(m_matrix_static_e[1]);
	EigenSparseFree(m_matrix_static_e[2]);
	EigenSparseFree(m_matrix_static_e[3]);
	EigenSparseFree(m_matrix_static_m[0]);
	EigenSparseFree(m_matrix_static_m[1]);
	EigenSparseFree(m_matrix_static_m[2]);
	EigenSparseFree(m_matrix_static_m[3]);
	EigenSparseFree(m_matrix_full_em[0]);
	EigenSparseFree(m_matrix_full_em[1]);
}

size_t GridMesh2D::GetPortCount() {
	return m_ports.size();
}

void GridMesh2D::InitGrid() {
	assert(!IsInitialized());

	// add grid lines
	std::vector<GridLine> grid_x, grid_y;
	GridAddBox(grid_x, grid_y, m_world_box, Box2D(REAL_MAX, REAL_MAX, REAL_MAX, REAL_MAX));
	GridAddBox(grid_x, grid_y, m_pml_box, m_pml_step);
	for(Conductor &conductor : m_conductors) {
		GridAddBox(grid_x, grid_y, conductor.m_box, conductor.m_step);
	}
	for(Dielectric &dielectric : m_dielectrics) {
		GridAddBox(grid_x, grid_y, dielectric.m_box, dielectric.m_step);
	}
	for(Box2D &integration_line : m_integration_lines) {
		GridAddBox(grid_x, grid_y, integration_line, Box2D(REAL_MAX, REAL_MAX, REAL_MAX, REAL_MAX));
	}

	// sort grid lines
	std::sort(grid_x.begin(), grid_x.end());
	std::sort(grid_y.begin(), grid_y.end());

	// calculate midpoints
	std::vector<real_t> midpoints_x, midpoints_y;
	GridMidpoints(midpoints_x, grid_x);
	GridMidpoints(midpoints_y, grid_y);

	// set maximum step size
	real_t max_step_vacuum = m_lambda_factor / (m_max_frequency * sqrt(VACUUM_PERMITTIVITY * VACUUM_PERMEABILITY));
	std::vector<real_t> max_step_x(grid_x.size() - 1, max_step_vacuum), max_step_y(grid_y.size() - 1, max_step_vacuum);
	for(Dielectric &dielectric : m_dielectrics) {
		MaterialDielectric::Properties properties;
		dielectric.m_material->GetProperties(m_max_frequency, properties);
		real_t step_x = m_lambda_factor / (m_max_frequency * sqrt(properties.permittivity_x.real() * VACUUM_PERMEABILITY));
		real_t step_y = m_lambda_factor / (m_max_frequency * sqrt(properties.permittivity_y.real() * VACUUM_PERMEABILITY));
		size_t ix1 = (size_t) (std::upper_bound(midpoints_x.begin(), midpoints_x.end(), dielectric.m_box.x1) - midpoints_x.begin());
		size_t ix2 = (size_t) (std::upper_bound(midpoints_x.begin(), midpoints_x.end(), dielectric.m_box.x2) - midpoints_x.begin());
		size_t iy1 = (size_t) (std::upper_bound(midpoints_y.begin(), midpoints_y.end(), dielectric.m_box.y1) - midpoints_y.begin());
		size_t iy2 = (size_t) (std::upper_bound(midpoints_y.begin(), midpoints_y.end(), dielectric.m_box.y2) - midpoints_y.begin());
		for(size_t ix = ix1; ix < ix2; ++ix) {
			max_step_x[ix] = std::min(max_step_x[ix], step_x);
		}
		for(size_t iy = iy1; iy < iy2; ++iy) {
			max_step_y[iy] = std::min(max_step_y[iy], step_y);
		}
	}

	// refine grid
	GridRefine(m_grid_x, grid_x, max_step_x, m_grid_inc, m_grid_epsilon);
	GridRefine(m_grid_y, grid_y, max_step_y, m_grid_inc, m_grid_epsilon);

	if(m_grid_x.size() < 2 || m_grid_y.size() < 2)
		throw std::runtime_error("GridMesh2D error: The mesh must have at least 2 grid lines.");

	// calculate midpoints
	GridMidpoints(m_midpoints_x, m_grid_x);
	GridMidpoints(m_midpoints_y, m_grid_y);

}

void GridMesh2D::InitCells() {
	assert(!IsInitialized());

	// allocate nodes and cells
	m_nodes.resize(m_grid_x.size() * m_grid_y.size());
	m_edges_x.resize((m_grid_x.size() - 1) * m_grid_y.size());
	m_edges_y.resize(m_grid_x.size() * (m_grid_y.size() - 1));
	m_cells.resize((m_grid_x.size() - 1) * (m_grid_y.size() - 1));

	// apply conductors
	for(size_t conductor_index = 0; conductor_index < m_conductors.size(); ++conductor_index) {
		Conductor &conductor = m_conductors[conductor_index];
		size_t ix1 = (size_t) (std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), conductor.m_box.x1) - m_midpoints_x.begin());
		size_t ix2 = (size_t) (std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), conductor.m_box.x2) - m_midpoints_x.begin());
		size_t iy1 = (size_t) (std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), conductor.m_box.y1) - m_midpoints_y.begin());
		size_t iy2 = (size_t) (std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), conductor.m_box.y2) - m_midpoints_y.begin());
		for(size_t iy = iy1; iy <= iy2; ++iy) {
			for(size_t ix = ix1; ix <= ix2; ++ix) {
				Node &node = GetNode(ix, iy);
				if(node.m_conductor != INDEX_NONE && m_conductors[node.m_conductor].m_port != conductor.m_port) {
					throw std::runtime_error(MakeString("Port ", conductor.m_port, " makes contact with port ", m_conductors[node.m_conductor].m_port, "."));
				}
				node.m_conductor = conductor_index;
			}
		}
		for(size_t iy = iy1; iy <= iy2; ++iy) {
			for(size_t ix = ix1; ix < ix2; ++ix) {
				Edge &edge = GetEdgeX(ix, iy);
				edge.m_conductor = conductor_index;
			}
		}
		for(size_t iy = iy1; iy < iy2; ++iy) {
			for(size_t ix = ix1; ix <= ix2; ++ix) {
				Edge &edge = GetEdgeY(ix, iy);
				edge.m_conductor = conductor_index;
			}
		}
		for(size_t iy = iy1; iy < iy2; ++iy) {
			for(size_t ix = ix1; ix < ix2; ++ix) {
				Cell &cell = GetCell(ix, iy);
				cell.m_conductor = conductor_index;
			}
		}
	}

	// apply dielectrics
	for(size_t dielectric_index = 0; dielectric_index < m_dielectrics.size(); ++dielectric_index) {
		Dielectric &dielectric = m_dielectrics[dielectric_index];
		size_t ix1 = (size_t) (std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), dielectric.m_box.x1) - m_midpoints_x.begin());
		size_t ix2 = (size_t) (std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), dielectric.m_box.x2) - m_midpoints_x.begin());
		size_t iy1 = (size_t) (std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), dielectric.m_box.y1) - m_midpoints_y.begin());
		size_t iy2 = (size_t) (std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), dielectric.m_box.y2) - m_midpoints_y.begin());
		for(size_t iy = iy1; iy < iy2; ++iy) {
			for(size_t ix = ix1; ix < ix2; ++ix) {
				Cell &cell = GetCell(ix, iy);
				cell.m_dielectric = dielectric_index;
			}
		}
	}

	// mark surfaces
	for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
		for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
			Cell &cell = GetCell(ix, iy);
			if(cell.m_conductor != INDEX_NONE)
				continue;
			Node &node00 = GetNode(ix    , iy    );
			Node &node01 = GetNode(ix + 1, iy    );
			Node &node10 = GetNode(ix    , iy + 1);
			Node &node11 = GetNode(ix + 1, iy + 1);
			Edge &edgex0 = GetEdgeX(ix    , iy    );
			Edge &edgex1 = GetEdgeX(ix    , iy + 1);
			Edge &edgey0 = GetEdgeY(ix    , iy    );
			Edge &edgey1 = GetEdgeY(ix + 1, iy    );
			if(node00.m_conductor != INDEX_NONE)
				node00.m_surface = true;
			if(node01.m_conductor != INDEX_NONE)
				node01.m_surface = true;
			if(node10.m_conductor != INDEX_NONE)
				node10.m_surface = true;
			if(node11.m_conductor != INDEX_NONE)
				node11.m_surface = true;
			if(edgex0.m_conductor != INDEX_NONE)
				edgex0.m_surface = true;
			if(edgex1.m_conductor != INDEX_NONE)
				edgex1.m_surface = true;
			if(edgey0.m_conductor != INDEX_NONE)
				edgey0.m_surface = true;
			if(edgey1.m_conductor != INDEX_NONE)
				edgey1.m_surface = true;
		}
	}

}

void GridMesh2D::InitVariables() {
	assert(!IsInitialized());

	if(m_ports.size() == 0)
		throw std::runtime_error("GridMesh2D error: The mesh has no ports.");

	// assign variables to ports
	for(size_t i = 0; i < m_ports.size(); ++i) {
		Port &port = m_ports[i];
		port.m_var_static_e = INDEX_OFFSET + m_vars_static_e_fixed++;
		port.m_var_static_m = INDEX_OFFSET + m_vars_static_m_fixed++;
	}

	// assign variables to nodes
	for(size_t iy = 0; iy < m_grid_y.size(); ++iy) {
		for(size_t ix = 0; ix < m_grid_x.size(); ++ix) {
			Node &node = GetNode(ix, iy);
			if(node.m_conductor == INDEX_NONE) {
				node.m_var_static_e1 = m_vars_static_e_free++;
				node.m_var_static_m1 = m_vars_static_m_free++;
			} else if(node.m_surface) {
				node.m_var_static_e1 = m_ports[m_conductors[node.m_conductor].m_port].m_var_static_e;
				node.m_var_static_m1 = m_vars_static_m_free++;
			} else {
				node.m_var_static_e1 = m_ports[m_conductors[node.m_conductor].m_port].m_var_static_e;
				node.m_var_static_m1 = m_ports[m_conductors[node.m_conductor].m_port].m_var_static_m;
			}
		}
	}

	// assign variables to edges
	if(m_element_type == ELEMENTTYPE_QUADRATIC) {
		for(size_t iy = 0; iy < m_grid_y.size(); ++iy) {
			for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
				Edge &edge = GetEdgeX(ix, iy);
				if(edge.m_conductor == INDEX_NONE) {
					edge.m_var_static_e2 = m_vars_static_e_free++;
					edge.m_var_static_m2 = m_vars_static_m_free++;
				} else if(edge.m_surface) {
					edge.m_var_static_m2 = m_vars_static_m_free++;
				}
			}
		}
		for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
			for(size_t ix = 0; ix < m_grid_x.size(); ++ix) {
				Edge &edge = GetEdgeY(ix, iy);
				if(edge.m_conductor == INDEX_NONE) {
					edge.m_var_static_e2 = m_vars_static_e_free++;
					edge.m_var_static_m2 = m_vars_static_m_free++;
				} else if(edge.m_surface) {
					edge.m_var_static_m2 = m_vars_static_m_free++;
				}
			}
		}
	}

	// assign EM variables
	if(m_solver_type == SOLVERTYPE_FULLWAVE) {

		// allocate flood fill queue
		size_t placeholder_count = 0;
		std::deque<std::pair<size_t, size_t>> tree_node_queue;

		// mark nodes outside conductors with placeholders
		for(size_t iy = 0; iy < m_grid_y.size(); ++iy) {
			for(size_t ix = 0; ix < m_grid_x.size(); ++ix) {
				Node &node = GetNode(ix, iy);
				if(node.m_conductor == INDEX_NONE) {
					node.m_var_full_m1 = INDEX_OFFSET; // placeholder
					++placeholder_count;
				} else {
					tree_node_queue.emplace_back(ix, iy);
				}
			}
		}

		// mark edges on or outside conductors with placeholders
		for(size_t iy = 0; iy < m_grid_y.size(); ++iy) {
			for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
				Edge &edge = GetEdgeX(ix, iy);
				if(edge.m_conductor == INDEX_NONE || edge.m_surface) {
					edge.m_var_full_mt1 = INDEX_OFFSET; // placeholder
				}
			}
		}
		for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
			for(size_t ix = 0; ix < m_grid_x.size(); ++ix) {
				Edge &edge = GetEdgeY(ix, iy);
				if(edge.m_conductor == INDEX_NONE || edge.m_surface) {
					edge.m_var_full_mt1 = INDEX_OFFSET; // placeholder
				}
			}
		}

		// mark integration lines
		for(Box2D &integration_line : m_integration_lines) {
			size_t ix1 = (size_t) (std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), integration_line.x1) - m_midpoints_x.begin());
			size_t ix2 = (size_t) (std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), integration_line.x2) - m_midpoints_x.begin());
			size_t iy1 = (size_t) (std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), integration_line.y1) - m_midpoints_y.begin());
			size_t iy2 = (size_t) (std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), integration_line.y2) - m_midpoints_y.begin());
			for(size_t iy = iy1; iy <= iy2; ++iy) {
				for(size_t ix = ix1; ix <= ix2; ++ix) {
					Node &node = GetNode(ix, iy);
					if(node.m_var_full_m1 == INDEX_OFFSET) {
						node.m_var_full_m1 = INDEX_NONE;
						--placeholder_count;
						tree_node_queue.emplace_back(ix, iy);
					}
				}
			}
			for(size_t iy = iy1; iy <= iy2; ++iy) {
				for(size_t ix = ix1; ix < ix2; ++ix) {
					Edge &edge = GetEdgeX(ix, iy);
					edge.m_var_full_mt1 = INDEX_NONE;
				}
			}
			for(size_t iy = iy1; iy < iy2; ++iy) {
				for(size_t ix = ix1; ix <= ix2; ++ix) {
					Edge &edge = GetEdgeY(ix, iy);
					edge.m_var_full_mt1 = INDEX_NONE;
				}
			}
		}

		// breadth-first flood fill
		while(!tree_node_queue.empty()) {
			size_t ix, iy;
			std::tie(ix, iy) = tree_node_queue.front();
			tree_node_queue.pop_front();
			if(ix > 0) {
				Node &node = GetNode(ix - 1, iy);
				if(node.m_var_full_m1 == INDEX_OFFSET) {
					node.m_var_full_m1 = INDEX_NONE;
					--placeholder_count;
					tree_node_queue.emplace_back(ix - 1, iy);
					Edge &edge = GetEdgeX(ix - 1, iy);
					edge.m_var_full_mt1 = INDEX_NONE;
				}
			}
			if(ix < m_grid_x.size() - 1) {
				Node &node = GetNode(ix + 1, iy);
				if(node.m_var_full_m1 == INDEX_OFFSET) {
					node.m_var_full_m1 = INDEX_NONE;
					--placeholder_count;
					tree_node_queue.emplace_back(ix + 1, iy);
					Edge &edge = GetEdgeX(ix, iy);
					edge.m_var_full_mt1 = INDEX_NONE;
				}
			}
			if(iy > 0) {
				Node &node = GetNode(ix, iy - 1);
				if(node.m_var_full_m1 == INDEX_OFFSET) {
					node.m_var_full_m1 = INDEX_NONE;
					--placeholder_count;
					tree_node_queue.emplace_back(ix, iy - 1);
					Edge &edge = GetEdgeY(ix, iy - 1);
					edge.m_var_full_mt1 = INDEX_NONE;
				}
			}
			if(iy < m_grid_y.size() - 1) {
				Node &node = GetNode(ix, iy + 1);
				if(node.m_var_full_m1 == INDEX_OFFSET) {
					node.m_var_full_m1 = INDEX_NONE;
					--placeholder_count;
					tree_node_queue.emplace_back(ix, iy + 1);
					Edge &edge = GetEdgeY(ix, iy);
					edge.m_var_full_mt1 = INDEX_NONE;
				}
			}
		}

		// sanity check
		if(placeholder_count != 0)
			throw std::runtime_error("GridMesh2D error: The mesh contains unreachable nodes.");

		// get absolute reference
		//size_t ref_ix = (size_t) (std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), m_ports[0].m_anchor.x) - m_midpoints_x.begin());
		//size_t ref_iy = (size_t) (std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), m_ports[0].m_anchor.y) - m_midpoints_y.begin());

		// assign variables to ports (except the first one)
		for(size_t i = 1; i < m_ports.size(); ++i) {
			Port &port = m_ports[i];
			port.m_var_full_e = m_vars_full_em++;
		}

		// assign variables to nodes
		for(size_t iy = 0; iy < m_grid_y.size(); ++iy) {
			for(size_t ix = 0; ix < m_grid_x.size(); ++ix) {
				Node &node = GetNode(ix, iy);
				if(node.m_conductor == INDEX_NONE) {
					node.m_var_full_e1 = m_vars_full_em++;
					node.m_var_full_m1 = m_vars_full_em++;
				} else {
					size_t port = m_conductors[node.m_conductor].m_port;
					node.m_var_full_e1 = m_ports[port].m_var_full_e;
					if(node.m_surface) {
						/*size_t ref_ix = (size_t) (std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), m_ports[port].m_anchor.x) - m_midpoints_x.begin());
						size_t ref_iy = (size_t) (std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), m_ports[port].m_anchor.y) - m_midpoints_y.begin());
						if(ix != ref_ix || iy != ref_iy) {
							node.m_var_full_e = m_vars_full++;
						}*/
						//node.m_var_full_e = m_vars_full++;
						node.m_var_full_m1 = m_vars_full_em++;
					} else {
						node.m_var_full_m1 = m_ports[port].m_var_full_e;
					}
				}
			}
		}

		// assign variables to edges
		for(size_t iy = 0; iy < m_grid_y.size(); ++iy) {
			for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
				Edge &edge = GetEdgeX(ix, iy);
				if(edge.m_var_full_mt1 == INDEX_OFFSET) {
					edge.m_var_full_mt1 = m_vars_full_em++;
					if(m_element_type == ELEMENTTYPE_QUADRATIC)
						edge.m_var_full_mt2 = m_vars_full_em++;
				}
			}
		}
		for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
			for(size_t ix = 0; ix < m_grid_x.size(); ++ix) {
				Edge &edge = GetEdgeY(ix, iy);
				if(edge.m_var_full_mt1 == INDEX_OFFSET) {
					edge.m_var_full_mt1 = m_vars_full_em++;
					if(m_element_type == ELEMENTTYPE_QUADRATIC)
						edge.m_var_full_mt2 = m_vars_full_em++;
				}
			}
		}
		if(m_element_type == ELEMENTTYPE_QUADRATIC) {
			for(size_t iy = 0; iy < m_grid_y.size(); ++iy) {
				for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
					Edge &edge = GetEdgeX(ix, iy);
					if(edge.m_conductor == INDEX_NONE) {
						edge.m_var_full_e2 = m_vars_full_em++;
						edge.m_var_full_m2 = m_vars_full_em++;
					} else if(edge.m_surface) {
						edge.m_var_full_m2 = m_vars_full_em++;
					}
				}
			}
			for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
				for(size_t ix = 0; ix < m_grid_x.size(); ++ix) {
					Edge &edge = GetEdgeY(ix, iy);
					if(edge.m_conductor == INDEX_NONE) {
						edge.m_var_full_e2 = m_vars_full_em++;
						edge.m_var_full_m2 = m_vars_full_em++;
					} else if(edge.m_surface) {
						edge.m_var_full_m2 = m_vars_full_em++;
					}
				}
			}
		}

	}

	// avoid problems later
	if(m_vars_static_e_free == 0)
		throw std::runtime_error("GridMesh2D error: There are no static-E free variables.");
	if(m_vars_static_e_fixed == 0)
		throw std::runtime_error("GridMesh2D error: There are no static-E fixed variables.");
	if(m_vars_static_m_free == 0)
		throw std::runtime_error("GridMesh2D error: There are no static-M free variables.");
	if(m_vars_static_m_fixed == 0)
		throw std::runtime_error("GridMesh2D error: There are no static-M fixed variables.");
	if(m_solver_type == SOLVERTYPE_FULLWAVE && m_vars_full_em == 0)
		throw std::runtime_error("GridMesh2D error: There are no full-EM variables.");

}

void GridMesh2D::BuildMatrices() {
	assert(IsInitialized() && !IsSolved());

	real_t omega = 2.0 * M_PI * GetFrequency();

	// load conductor properties
	m_conductor_properties.clear();
	m_conductor_properties.resize(m_conductors.size());
	for(size_t i = 0; i < m_conductors.size(); ++i) {
		m_conductors[i].m_material->GetProperties(GetFrequency(), m_conductor_properties[i]);
	}

	// load dielectric properties
	m_dielectric_properties.clear();
	m_dielectric_properties.resize(m_dielectrics.size());
	for(size_t i = 0; i < m_dielectrics.size(); ++i) {
		m_dielectrics[i].m_material->GetProperties(GetFrequency(), m_dielectric_properties[i]);
	}

	// prepare PML
	/*complex_t pml_mult_x1 = m_pml_attenuation / ((m_pml_box.x1 - m_world_box.x1) * complex_t(1.0, -2.0 * GetFrequency() * (m_pml_box.x1 - m_world_box.x1) / SPEED_OF_LIGHT));
	complex_t pml_mult_x2 = m_pml_attenuation / ((m_world_box.x2 - m_pml_box.x2) * complex_t(1.0, -2.0 * GetFrequency() * (m_world_box.x2 - m_pml_box.x2) / SPEED_OF_LIGHT));
	complex_t pml_mult_y1 = m_pml_attenuation / ((m_pml_box.y1 - m_world_box.y1) * complex_t(1.0, -2.0 * GetFrequency() * (m_pml_box.y1 - m_world_box.y1) / SPEED_OF_LIGHT));
	complex_t pml_mult_y2 = m_pml_attenuation / ((m_world_box.y2 - m_pml_box.y2) * complex_t(1.0, -2.0 * GetFrequency() * (m_world_box.y2 - m_pml_box.y2) / SPEED_OF_LIGHT));*/

	// allocate matrices
	std::vector<real_t> port_dc_conductances(m_ports.size(), 0.0);
	SparseBlockMatrix<complex_t> matrix_static_e, matrix_static_m;
	SparseMatrix<complex_t> matrix_full_em[2];

	matrix_static_e.Reset(m_vars_static_e_free, m_vars_static_e_fixed, m_vars_static_e_free, m_vars_static_e_fixed);
	matrix_static_m.Reset(m_vars_static_m_free, m_vars_static_m_fixed, m_vars_static_m_free, m_vars_static_m_fixed);
	matrix_full_em[0].Reset(m_vars_full_em, m_vars_full_em);
	matrix_full_em[1].Reset(m_vars_full_em, m_vars_full_em);

	// build matrices
	for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
		for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {

			// get cell
			Cell &cell = GetCell(ix, iy);

			// get cell size
			real_t delta_x = m_grid_x[ix + 1] - m_grid_x[ix];
			real_t delta_y = m_grid_y[iy + 1] - m_grid_y[iy];

			// skip cells that are inside conductors
			if(cell.m_conductor != INDEX_NONE) {
				size_t port = m_conductors[cell.m_conductor].m_port;
				real_t conductivity = m_conductor_properties[cell.m_conductor].conductivity;
				port_dc_conductances[port] += conductivity * delta_x * delta_y;
				continue;
			}

			// get neighboring nodes and edges
			Node &node00 = GetNode(ix    , iy    );
			Node &node01 = GetNode(ix + 1, iy    );
			Node &node10 = GetNode(ix    , iy + 1);
			Node &node11 = GetNode(ix + 1, iy + 1);
			Edge &edgex0 = GetEdgeX(ix    , iy    );
			Edge &edgex1 = GetEdgeX(ix    , iy + 1);
			Edge &edgey0 = GetEdgeY(ix    , iy    );
			Edge &edgey1 = GetEdgeY(ix + 1, iy    );

			// get dielectric properties
			complex_t permittivity_x(VACUUM_PERMITTIVITY, 0.0), permittivity_y(VACUUM_PERMITTIVITY, 0.0), permittivity_z(VACUUM_PERMITTIVITY, 0.0);
			if(cell.m_dielectric != INDEX_NONE) {
				permittivity_x = m_dielectric_properties[cell.m_dielectric].permittivity_x;
				permittivity_y = m_dielectric_properties[cell.m_dielectric].permittivity_y;
				permittivity_z = m_dielectric_properties[cell.m_dielectric].permittivity_z;
			}
			complex_t permeability_x(VACUUM_PERMEABILITY, 0.0), permeability_y(VACUUM_PERMEABILITY, 0.0), permeability_z(VACUUM_PERMEABILITY, 0.0);

			// get PML properties
			/*real_t center_x = 0.5 * (m_grid_x[ix] + m_grid_x[ix + 1]);
			real_t center_y = 0.5 * (m_grid_y[iy] + m_grid_y[iy + 1]);
			complex_t pml_sx = 1.0, pml_sy = 1.0;
			if(center_x < m_pml_box.x1) pml_sx += (m_pml_box.x1 - center_x) * pml_mult_x1;
			if(center_x > m_pml_box.x2) pml_sx += (center_x - m_pml_box.x2) * pml_mult_x2;
			if(center_y < m_pml_box.y1) pml_sy += (m_pml_box.y1 - center_y) * pml_mult_y1;
			if(center_y > m_pml_box.y2) pml_sy += (center_y - m_pml_box.y2) * pml_mult_y2;
			complex_t delta_pml_x = delta_x * pml_sx;
			complex_t delta_pml_y = delta_y * pml_sy;*/

			// add to static E matrix
			size_t vars_static_epot_rect[8] = {
				node00.m_var_static_e1, node01.m_var_static_e1, node10.m_var_static_e1, node11.m_var_static_e1,
				edgex0.m_var_static_e2, edgex1.m_var_static_e2, edgey0.m_var_static_e2, edgey1.m_var_static_e2,
			};
			FemMatrix_StaticEPot_Rect(matrix_static_e, vars_static_epot_rect, delta_x, delta_y, permittivity_x, permittivity_y);

			// add to static M matrix
			size_t vars_static_mpot_rect[8] = {
				node00.m_var_static_m1, node01.m_var_static_m1, node10.m_var_static_m1, node11.m_var_static_m1,
				edgex0.m_var_static_m2, edgex1.m_var_static_m2, edgey0.m_var_static_m2, edgey1.m_var_static_m2,
			};
			FemMatrix_StaticMPot_Rect(matrix_static_m, vars_static_mpot_rect, delta_x, delta_y, permeability_x, permeability_y);

			// add to full EM matrix
			if(m_solver_type == SOLVERTYPE_FULLWAVE) {
				size_t vars_full_empot_rect[24] = {
					node00.m_var_full_e1, node01.m_var_full_e1, node10.m_var_full_e1, node11.m_var_full_e1,
					edgex0.m_var_full_e2, edgex1.m_var_full_e2, edgey0.m_var_full_e2, edgey1.m_var_full_e2,
					edgex0.m_var_full_mt1, edgex1.m_var_full_mt1, edgex0.m_var_full_mt2, edgex1.m_var_full_mt2,
					edgey0.m_var_full_mt1, edgey1.m_var_full_mt1, edgey0.m_var_full_mt2, edgey1.m_var_full_mt2,
					node00.m_var_full_m1, node01.m_var_full_m1, node10.m_var_full_m1, node11.m_var_full_m1,
					edgex0.m_var_full_m2, edgex1.m_var_full_m2, edgey0.m_var_full_m2, edgey1.m_var_full_m2,
				};
				FemMatrix_FullEMPot_Rect(matrix_full_em, vars_full_empot_rect, delta_x, delta_y, omega,
						permittivity_x, permittivity_y, permittivity_z,
						permeability_x, permeability_y, permeability_z);
			}

			// process conductor surfaces
			if(edgex0.m_conductor != INDEX_NONE) {
				complex_t impedance = m_conductor_properties[edgex0.m_conductor].impedance;
				size_t vars_static_mpot_xline[4] = {
					node00.m_var_static_m1, node01.m_var_static_m1, edgex0.m_var_static_m2,
					m_ports[m_conductors[edgex0.m_conductor].m_port].m_var_static_m,
				};
				FemMatrix_StaticMPot_XLine(matrix_static_m, vars_static_mpot_xline, delta_x, omega, impedance);
				if(m_solver_type == SOLVERTYPE_FULLWAVE) {
					size_t vars_full_empot_xline[8] = {
						node00.m_var_full_e1, node01.m_var_full_e1, edgex0.m_var_full_e2,
						edgex0.m_var_full_mt1, edgex0.m_var_full_mt2,
						node00.m_var_full_m1, node01.m_var_full_m1, edgex0.m_var_full_m2,
					};
					FemMatrix_FullEMPot_XLine(matrix_full_em, vars_full_empot_xline, delta_x, omega, impedance);
				}
			}
			if(edgex1.m_conductor != INDEX_NONE) {
				complex_t impedance = m_conductor_properties[edgex1.m_conductor].impedance;
				size_t vars_static_mpot_xline[4] = {
					node10.m_var_static_m1, node11.m_var_static_m1, edgex1.m_var_static_m2,
					m_ports[m_conductors[edgex1.m_conductor].m_port].m_var_static_m,
				};
				FemMatrix_StaticMPot_XLine(matrix_static_m, vars_static_mpot_xline, delta_x, omega, impedance);
				if(m_solver_type == SOLVERTYPE_FULLWAVE) {
					size_t vars_full_empot_xline[8] = {
						node10.m_var_full_e1, node11.m_var_full_e1, edgex1.m_var_full_e2,
						edgex1.m_var_full_mt1, edgex1.m_var_full_mt2,
						node10.m_var_full_m1, node11.m_var_full_m1, edgex1.m_var_full_m2,
					};
					FemMatrix_FullEMPot_XLine(matrix_full_em, vars_full_empot_xline, delta_x, omega, impedance);
				}
			}
			if(edgey0.m_conductor != INDEX_NONE) {
				complex_t impedance = m_conductor_properties[edgey0.m_conductor].impedance;
				size_t vars_static_mpot_yline[4] = {
					node00.m_var_static_m1, node10.m_var_static_m1, edgey0.m_var_static_m2,
					m_ports[m_conductors[edgey0.m_conductor].m_port].m_var_static_m,
				};
				FemMatrix_StaticMPot_YLine(matrix_static_m, vars_static_mpot_yline, delta_y, omega, impedance);
				if(m_solver_type == SOLVERTYPE_FULLWAVE) {
					size_t vars_full_empot_yline[8] = {
						node00.m_var_full_e1, node10.m_var_full_e1, edgey0.m_var_full_e2,
						edgey0.m_var_full_mt1, edgey0.m_var_full_mt2,
						node00.m_var_full_m1, node10.m_var_full_m1, edgey0.m_var_full_m2,
					};
					FemMatrix_FullEMPot_YLine(matrix_full_em, vars_full_empot_yline, delta_y, omega, impedance);
				}
			}
			if(edgey1.m_conductor != INDEX_NONE) {
				complex_t impedance = m_conductor_properties[edgey1.m_conductor].impedance;
				size_t vars_static_mpot_yline[4] = {
					node01.m_var_static_m1, node11.m_var_static_m1, edgey1.m_var_static_m2,
					m_ports[m_conductors[edgey1.m_conductor].m_port].m_var_static_m,
				};
				FemMatrix_StaticMPot_YLine(matrix_static_m, vars_static_mpot_yline, delta_y, omega, impedance);
				if(m_solver_type == SOLVERTYPE_FULLWAVE) {
					size_t vars_full_empot_yline[8] = {
						node01.m_var_full_e1, node11.m_var_full_e1, edgey1.m_var_full_e2,
						edgey1.m_var_full_mt1, edgey1.m_var_full_mt2,
						node01.m_var_full_m1, node11.m_var_full_m1, edgey1.m_var_full_m2,
					};
					FemMatrix_FullEMPot_YLine(matrix_full_em, vars_full_empot_yline, delta_y, omega, impedance);
				}
			}

		}
	}

	// convert port DC conductance to resistance
	m_vector_dc_resistances.resize((Eigen::Index) m_ports.size());
	for(size_t i = 0; i < m_ports.size(); ++i) {
		m_vector_dc_resistances[(Eigen::Index) i] = (m_ports[i].m_infinite_area)? 0.0 : 1.0 / port_dc_conductances[i];
	}

	// convert to Eigen sparse matrices
	matrix_static_e.GetMatrixA().ToEigen(m_matrix_static_e[0]);
	matrix_static_e.GetMatrixB().ToEigen(m_matrix_static_e[1]);
	matrix_static_e.GetMatrixC().ToEigen(m_matrix_static_e[2]);
	matrix_static_e.GetMatrixD().ToEigen(m_matrix_static_e[3]);
	matrix_static_m.GetMatrixA().ToEigen(m_matrix_static_m[0]);
	matrix_static_m.GetMatrixB().ToEigen(m_matrix_static_m[1]);
	matrix_static_m.GetMatrixC().ToEigen(m_matrix_static_m[2]);
	matrix_static_m.GetMatrixD().ToEigen(m_matrix_static_m[3]);
	if(m_solver_type == SOLVERTYPE_FULLWAVE) {
		matrix_full_em[0].ToEigen(m_matrix_full_em[0]);
		matrix_full_em[1].ToEigen(m_matrix_full_em[1]);
	}

#if SIMULATION_SAVE_MATRIXMARKET
	MatrixMarket::Save("matrix_epot.mtx", m_matrix_static_e[0], MatrixMarket::TYPE_SYMMETRIC);
	MatrixMarket::Save("matrix_mpot.mtx", m_matrix_static_m[0], MatrixMarket::TYPE_SYMMETRIC);
	if(m_solver_type == SOLVERTYPE_FULLWAVE) {
		MatrixMarket::Save("matrix_empot0.mtx", m_matrix_full_em[0], MatrixMarket::TYPE_GENERAL);
		MatrixMarket::Save("matrix_empot1.mtx", m_matrix_full_em[1], MatrixMarket::TYPE_GENERAL);
	}
#endif

}

void GridMesh2D::SolveStaticModes() {
	assert(IsInitialized() && !IsSolved());

	real_t omega = 2.0 * M_PI * GetFrequency();

	// factorize electric potential matrix
	if(m_lu_static_e.colsPermutation().size() == 0) {
		m_lu_static_e.isSymmetric(true);
		m_lu_static_e.setPivotThreshold(0.0);
		m_lu_static_e.analyzePattern(m_matrix_static_e[0]);
	}
	m_lu_static_e.factorize(m_matrix_static_e[0]);
	if(m_lu_static_e.info() != Eigen::Success)
		throw std::runtime_error("Static E matrix factorization failed!");

	// solve electric potential matrix
	Eigen::MatrixXc rhs_static_e = -m_matrix_static_e[1] * GetModes();
	m_solution_static_e = m_lu_static_e.solve(rhs_static_e);

	// factorize magnetic potential matrix
	if(m_lu_static_m.colsPermutation().size() == 0) {
		m_lu_static_m.isSymmetric(true);
		m_lu_static_m.setPivotThreshold(0.0);
		m_lu_static_m.analyzePattern(m_matrix_static_m[0]);
	}
	m_lu_static_m.factorize(m_matrix_static_m[0]);
	if(m_lu_static_m.info() != Eigen::Success)
		throw std::runtime_error("Static M matrix factorization failed!");

	// solve magnetic potential matrix
	Eigen::MatrixXc rhs_static_m = -m_matrix_static_m[1] * GetModes();
	m_solution_static_m = m_lu_static_m.solve(rhs_static_m);

	// calculate residuals
	Eigen::MatrixXc residual_epot = m_matrix_static_e[2] * m_solution_static_e + m_matrix_static_e[3] * GetModes();
	Eigen::MatrixXc residual_mpot = m_matrix_static_m[2] * m_solution_static_m + m_matrix_static_m[3] * GetModes();
	Eigen::MatrixXc charge_matrix = GetModes().transpose() * residual_epot;
	Eigen::MatrixXc current_matrix = GetModes().transpose() * residual_mpot;
	Eigen::MatrixXc current_matrix_inv = current_matrix.inverse();

	// calculate DC losses and combine with surface losses
	/*Eigen::MatrixXr dc_loss_matrix = residual_mpot.transpose() * m_vector_dc_resistances.asDiagonal() * residual_mpot;
	Eigen::MatrixXr combined_loss_matrix = surface_loss_matrix.cwiseMax(dc_loss_matrix);*/

	std::cerr << "charges =\n" << charge_matrix << std::endl;
	std::cerr << "currents =\n" << current_matrix << std::endl;
	std::cerr << std::endl;

	// calculate L, C, R and G
	m_inductance_matrix = current_matrix_inv.real();
	m_capacitance_matrix = charge_matrix.real();
	m_resistance_matrix = -omega * current_matrix_inv.imag();
	m_conductance_matrix = -omega * charge_matrix.imag();

}

void GridMesh2D::SolveStaticEigenModes() {

	std::cerr << "inductance =\n" << m_inductance_matrix << std::endl;
	std::cerr << "capacitance =\n" << m_capacitance_matrix << std::endl;
	std::cerr << "resistance =\n" << m_resistance_matrix << std::endl;
	std::cerr << "conductance =\n" << m_conductance_matrix << std::endl;
	std::cerr << std::endl;

	real_t omega = 2.0 * M_PI * GetFrequency();

	// convert to impedance and admittance matrices
	Eigen::MatrixXc impedance(GetModeCount(), GetModeCount()), admittance(GetModeCount(), GetModeCount());
	impedance.real() = m_resistance_matrix;
	impedance.imag() = omega * m_inductance_matrix;
	admittance.real() = m_conductance_matrix;
	admittance.imag() = omega * m_capacitance_matrix;

	// calculate eigenmodes
	Eigen::MatrixXc matrix_zy = impedance * admittance;
	Eigen::ComplexEigenSolver<Eigen::MatrixXc> eigensolver(matrix_zy);
	if(eigensolver.info() != Eigen::Success)
		throw std::runtime_error("Eigenmode decomposition failed!");
	auto &eigval = eigensolver.eigenvalues();
	auto &eigvec = eigensolver.eigenvectors();

	std::cerr << "eigenvalues =\n" << eigval << std::endl;
	std::cerr << "eigenvectors =\n" << eigvec << std::endl;
	std::cerr << std::endl;

	// make sure we get the correct complex square root (positive imaginary part)
	Eigen::VectorXc eigval_sqrt = (-eigval).cwiseSqrt() * complex_t(0.0, 1.0);

	// calculate impedances, admittance and (approximated) characteristic impedance
	Eigen::MatrixXc matrix_zy_invsqrt = eigvec * eigval_sqrt.cwiseInverse().asDiagonal() * eigvec.inverse();
	m_characteristic_impedance_matrix = matrix_zy_invsqrt * impedance;
	Eigen::VectorXc diag_impedance = m_characteristic_impedance_matrix.diagonal();
	Eigen::VectorXc diag_admittance = m_characteristic_impedance_matrix.inverse().diagonal();
	m_characteristic_impedances = (diag_impedance.array() / diag_admittance.array()).sqrt().matrix();

	// calculate propagation constants
	m_propagation_constants = (-matrix_zy.diagonal()).cwiseSqrt() * complex_t(0.0, 1.0);

	std::cerr << "m_characteristic_impedance_matrix =\n" << m_characteristic_impedance_matrix << std::endl;
	std::cerr << "m_characteristic_impedances =\n" << m_characteristic_impedances << std::endl;
	std::cerr << "m_propagation_constants =\n" << m_propagation_constants << std::endl;
	std::cerr << std::endl;

	// reorder eigenmodes to match user-provided modes and calculate eigenmode propagation constants
	m_eigenmodes.resize((Eigen::Index) GetModeCount(), (Eigen::Index) GetModeCount());
	m_eigenmode_propagation_constants.resize((Eigen::Index) GetModeCount());
	std::vector<size_t> modemap(GetModeCount());
	for(size_t i = 0; i < GetModeCount(); ++i) {
		modemap[i] = i;
	}
	for(size_t i = 0; i < GetModeCount(); ++i) {
		size_t best_index = i;
		complex_t best_value = complex_t(0.0, 0.0);
		for(size_t j = i; j < GetModeCount(); ++j) {
			complex_t value = eigvec((Eigen::Index) i, (Eigen::Index) modemap[j]);
			if(std::norm(value) > std::norm(best_value)) {
				best_index = j;
				best_value = value;
			}
		}
		std::swap(modemap[i], modemap[best_index]);
		m_eigenmodes.col((Eigen::Index) i) = eigvec.col((Eigen::Index) modemap[i]) / best_value;
		m_eigenmode_propagation_constants[(Eigen::Index) i] = eigval_sqrt[(Eigen::Index) modemap[i]];
	}

	std::cerr << "m_eigenmodes =\n" << m_eigenmodes << std::endl;
	std::cerr << "m_eigenmode_propagation_constants =\n" << m_eigenmode_propagation_constants << std::endl;
	std::cerr << std::endl;

	// save fields
	Eigen::VectorXc eigenmode = m_eigenmodes.col(0);
	Eigen::VectorXc fixed_values = GetModes() * eigenmode;
	Eigen::VectorXc solution_e = m_solution_static_e * eigenmode;
	Eigen::VectorXc solution_m = m_solution_static_m * eigenmode;
	m_solution_fields.resize(GetModeCount());
	for(size_t mode = 0; mode < GetModeCount(); ++mode) {
		SolutionField &fields = m_solution_fields[mode];
		fields.m_propagation_constant = m_eigenmode_propagation_constants[(Eigen::Index) mode];
		fields.m_effective_index = fields.m_propagation_constant * SPEED_OF_LIGHT / complex_t(0.0, omega);
		fields.m_field_nodes.resize(m_nodes.size());
		fields.m_field_edges_x.resize(m_edges_x.size());
		fields.m_field_edges_y.resize(m_edges_y.size());
		complex_t scale_factor_e = 1.0;
		complex_t scale_factor_m = fields.m_effective_index / SPEED_OF_LIGHT;
		for(size_t i = 0; i < m_nodes.size(); ++i) {
			Node &node = m_nodes[i];
			fields.m_field_nodes[i].e1 = ((node.m_var_static_e1 < INDEX_OFFSET)? solution_e[(Eigen::Index) node.m_var_static_e1] : fixed_values[(Eigen::Index) (node.m_var_static_e1 - INDEX_OFFSET)]) * scale_factor_e;
			fields.m_field_nodes[i].m1 = ((node.m_var_static_m1 < INDEX_OFFSET)? solution_m[(Eigen::Index) node.m_var_static_m1] : fixed_values[(Eigen::Index) (node.m_var_static_m1 - INDEX_OFFSET)]) * scale_factor_m;
		}
		for(size_t i = 0; i < m_edges_x.size(); ++i) {
			Edge &edge = m_edges_x[i];
			fields.m_field_edges_x[i].e2 = (edge.m_var_static_e2 == INDEX_NONE)? 0.0 : solution_e[(Eigen::Index) edge.m_var_static_e2] * scale_factor_e;
			fields.m_field_edges_x[i].m2 = (edge.m_var_static_m2 == INDEX_NONE)? 0.0 : solution_m[(Eigen::Index) edge.m_var_static_m2] * scale_factor_m;
			fields.m_field_edges_x[i].mt1 = 0.0;
			fields.m_field_edges_x[i].mt2 = 0.0;
		}
		for(size_t i = 0; i < m_edges_y.size(); ++i) {
			Edge &edge = m_edges_y[i];
			fields.m_field_edges_y[i].e2 = (edge.m_var_static_e2 == INDEX_NONE)? 0.0 : solution_e[(Eigen::Index) edge.m_var_static_e2] * scale_factor_e;
			fields.m_field_edges_y[i].m2 = (edge.m_var_static_m2 == INDEX_NONE)? 0.0 : solution_m[(Eigen::Index) edge.m_var_static_m2] * scale_factor_m;
			fields.m_field_edges_y[i].mt1 = 0.0;
			fields.m_field_edges_y[i].mt2 = 0.0;
		}
	}

}

void GridMesh2D::SolveFullEigenModes() {

	// skip if the full solver is disabled
	if(m_solver_type != SOLVERTYPE_FULLWAVE)
		return;

	real_t omega = 2.0 * M_PI * GetFrequency();

	// refine eigenmodes from static solver
	for(size_t mode = 0; mode < GetModeCount(); ++mode) {
		SolutionField &fields = m_solution_fields[mode];

		// get static eigenmode and propagation constant
		Eigen::VectorXc eigenmode = m_eigenmodes.col((Eigen::Index) mode);
		complex_t eigenvalue = fields.m_effective_index;

		// calculate static potentials for this mode
		Eigen::VectorXc fixed_values = GetModes() * eigenmode;
		Eigen::VectorXc static_epot = m_solution_static_e * eigenmode;
		Eigen::VectorXc static_mpot = m_solution_static_m * eigenmode;

		complex_t scale_factor_e = 1.0;
		complex_t scale_factor_m = fields.m_effective_index / SPEED_OF_LIGHT;
		complex_t scale_factor_mt = 1.0 / SPEED_OF_LIGHT;

		// calculate initial guess
		// TODO: deal with duplicates
		Eigen::VectorXc eigenvector(m_vars_full_em);
		Eigen::VectorXc orthofactors(m_vars_full_em);
		orthofactors.fill(1.0);
		for(size_t i = 0; i < m_nodes.size(); ++i) {
			Node &node = m_nodes[i];
			if(node.m_conductor == INDEX_NONE || node.m_surface) {
				if(node.m_var_full_e1 != INDEX_NONE) {
					eigenvector[(Eigen::Index) node.m_var_full_e1] = fields.m_field_nodes[i].e1 / scale_factor_e;
				}
				if(node.m_var_full_m1 != INDEX_NONE) {
					eigenvector[(Eigen::Index) node.m_var_full_m1] = fields.m_field_nodes[i].m1 / scale_factor_m;
					orthofactors[(Eigen::Index) node.m_var_full_m1] = 0.0;
				}
			}
		}
		for(size_t i = 0; i < m_edges_x.size(); ++i) {
			Edge &edge = m_edges_x[i];
			if(edge.m_conductor == INDEX_NONE || edge.m_surface) {
				if(edge.m_var_full_e2 != INDEX_NONE) {
					eigenvector[(Eigen::Index) edge.m_var_full_e2] = fields.m_field_edges_x[i].e2 / scale_factor_e;
				}
				if(edge.m_var_full_m2 != INDEX_NONE) {
					eigenvector[(Eigen::Index) edge.m_var_full_m2] = fields.m_field_edges_x[i].m2 / scale_factor_m;
					orthofactors[(Eigen::Index) edge.m_var_full_m2] = 0.0;
				}
				if(edge.m_var_full_mt1 != INDEX_NONE) {
					eigenvector[(Eigen::Index) edge.m_var_full_mt1] = fields.m_field_edges_x[i].mt1 / scale_factor_mt;
				}
				if(edge.m_var_full_mt2 != INDEX_NONE) {
					eigenvector[(Eigen::Index) edge.m_var_full_mt2] = fields.m_field_edges_x[i].mt2 / scale_factor_mt;
				}
			}
		}
		for(size_t i = 0; i < m_edges_y.size(); ++i) {
			Edge &edge = m_edges_y[i];
			if(edge.m_conductor == INDEX_NONE || edge.m_surface) {
				if(edge.m_var_full_e2 != INDEX_NONE) {
					eigenvector[(Eigen::Index) edge.m_var_full_e2] = fields.m_field_edges_y[i].e2 / scale_factor_e;
				}
				if(edge.m_var_full_m2 != INDEX_NONE) {
					eigenvector[(Eigen::Index) edge.m_var_full_m2] = fields.m_field_edges_y[i].m2 / scale_factor_m;
					orthofactors[(Eigen::Index) edge.m_var_full_m2] = 0.0;
				}
				if(edge.m_var_full_mt1 != INDEX_NONE) {
					eigenvector[(Eigen::Index) edge.m_var_full_mt1] = fields.m_field_edges_y[i].mt1 / scale_factor_mt;
				}
				if(edge.m_var_full_mt2 != INDEX_NONE) {
					eigenvector[(Eigen::Index) edge.m_var_full_mt2] = fields.m_field_edges_y[i].mt2 / scale_factor_mt;
				}
			}
		}
		eigenvector *= 1.0 / eigenvector.norm();

		// calculate residual
		Eigen::VectorXc resid = m_matrix_full_em[0] * eigenvector + square(eigenvalue) * (m_matrix_full_em[1] * eigenvector);
		std::cerr << "eigenvalue = " << eigenvalue << ", norm(resid) = " << resid.norm() << std::endl;

		// M = [1] = B
		// Q = [0] = A

		//propagation_constant = 0.0;
		//propagation_constant.real(0.0);
		//propagation_constant.imag(propagation_constant.imag() * 0.9);

		Eigen::VectorXc eigenvector_init = eigenvector;

		size_t iters = 3; // 4
		for(size_t it = 0; it < iters; ++it) {

			// improve eigenvalue
			if(it != 0) {
				complex_t b = eigenvector.transpose() * m_matrix_full_em[1] * eigenvector;
				complex_t a = eigenvector.transpose() * m_matrix_full_em[0] * eigenvector;
				eigenvalue = std::sqrt(-a / b);
			}

			// calculate residual
			resid = m_matrix_full_em[0] * eigenvector + square(eigenvalue) * (m_matrix_full_em[1] * eigenvector);
			std::cerr << "eigenvalue = " << eigenvalue << ", norm(resid) = " << resid.norm() << ", dot(eig, init) = " << ((eigenvector_init.adjoint() * eigenvector)[0]) << std::endl;

			if(it == iters - 1)
				break;

			// generate updated matrix
			Eigen::SparseMatrix<complex_t> mat = m_matrix_full_em[0] + square(eigenvalue) * m_matrix_full_em[1];

			// factorize matrix
			if(m_lu_full_em.colsPermutation().size() == 0) {
				m_lu_full_em.isSymmetric(true);
				m_lu_full_em.setPivotThreshold(0.0);
				m_lu_full_em.analyzePattern(mat);
			}
			m_lu_full_em.factorize(mat);
			if(m_lu_full_em.info() != Eigen::Success)
				throw std::runtime_error("Sparse matrix factorization failed!");

			// improve eigenvector
			for(size_t k = 0; k < 10; ++k) {
				Eigen::VectorXc temp1 = m_matrix_full_em[1] * eigenvector;
				Eigen::VectorXc temp2 = temp1.cwiseProduct(orthofactors);
				//std::cerr << "ortho: " << temp1.norm() << " " << temp2.norm() << " " << (temp1 - temp2).norm() << std::endl;
				eigenvector = m_lu_full_em.solve(temp2);
				if(m_lu_full_em.info() != Eigen::Success)
					throw std::runtime_error("Sparse matrix solving failed!");
				eigenvector *= 1.0 / eigenvector.norm();
			}

		}

		// normalize phase
		// TODO: change this?
		eigenvector *= 1.0 / std::sqrt((eigenvector.transpose() * eigenvector)[0]);

		// norm of parts
		complex_t norms_e = 0.0, norms_m = 0.0, norms_mt = 0.0;
		complex_t resid_e = 0.0, resid_m = 0.0, resid_mt = 0.0;
		for(Node &node : m_nodes) {
			if(node.m_var_full_e1 != INDEX_NONE) {
				norms_e += square(eigenvector[(Eigen::Index) node.m_var_full_e1]);
				resid_e += square(resid((Eigen::Index) node.m_var_full_e1, 0));
			}
			if(node.m_var_full_m1 != INDEX_NONE) {
				norms_m += square(eigenvector[(Eigen::Index) node.m_var_full_m1]);
				resid_m += square(resid((Eigen::Index) node.m_var_full_m1, 0));
			}
		}
		for(Edge &edge : m_edges_x) {
			if(edge.m_var_full_mt1 != INDEX_NONE) {
				norms_mt += square(eigenvector[(Eigen::Index) edge.m_var_full_mt1]);
				resid_mt += square(resid((Eigen::Index) edge.m_var_full_mt1, 0));
			}
		}
		for(Edge &edge : m_edges_y) {
			if(edge.m_var_full_mt1 != INDEX_NONE) {
				norms_mt += square(eigenvector[(Eigen::Index) edge.m_var_full_mt1]);
				resid_mt += square(resid((Eigen::Index) edge.m_var_full_mt1, 0));
			}
		}
		std::cerr << "norms: e = " << std::sqrt(norms_e) << ", m = " << std::sqrt(norms_m) << ", mt = " << std::sqrt(norms_mt) << std::endl;
		std::cerr << "resid: e = " << std::sqrt(resid_e) << ", m = " << std::sqrt(resid_m) << ", mt = " << std::sqrt(resid_mt) << std::endl;

		// update scale factors
		/*Eigen::VectorXc scalefactors2(m_vars_full);
		scalefactors2.fill(1.0);
		for(Node &node : m_nodes) {
			if(node.m_var_full_e != INDEX_NONE) {
				scalefactors2[(Eigen::Index) node.m_var_full_e] = complex_t(0.0, 1.0);
			}
			if(node.m_var_full_m != INDEX_NONE) {
				scalefactors2[(Eigen::Index) node.m_var_full_m] = -propagation_constant / omega;
			}
		}
		scalefactors = scalefactors.cwiseProduct(scalefactors2);*/

		scale_factor_e = 1.0;
		scale_factor_m = eigenvalue / SPEED_OF_LIGHT;
		scale_factor_mt = 1.0 / SPEED_OF_LIGHT;

		//std::cerr << "eigenvector_init:" << eigenvector_init << std::endl;
		//std::cerr << "eigenvector:" << eigenvector << std::endl;

		// TODO: remove
		m_propagation_constants[(Eigen::Index) mode] = eigenvalue / SPEED_OF_LIGHT * complex_t(0.0, omega);

		// save fields
		fields.m_propagation_constant = eigenvalue / SPEED_OF_LIGHT * complex_t(0.0, omega);
		fields.m_effective_index = fields.m_propagation_constant * SPEED_OF_LIGHT / complex_t(0.0, omega);
		fields.m_field_nodes.resize(m_nodes.size());
		fields.m_field_edges_x.resize(m_edges_x.size());
		fields.m_field_edges_y.resize(m_edges_y.size());
		for(size_t i = 0; i < m_nodes.size(); ++i) {
			Node &node = m_nodes[i];
			fields.m_field_nodes[i].e1 = (node.m_var_full_e1 == INDEX_NONE)? 0.0 : eigenvector[(Eigen::Index) node.m_var_full_e1] * scale_factor_e;
			fields.m_field_nodes[i].m1 = (node.m_var_full_m1 == INDEX_NONE)? 0.0 : eigenvector[(Eigen::Index) node.m_var_full_m1] * scale_factor_m;
		}
		for(size_t i = 0; i < m_edges_x.size(); ++i) {
			Edge &edge = m_edges_x[i];
			fields.m_field_edges_x[i].e2 = (edge.m_var_full_e2 == INDEX_NONE)? 0.0 : eigenvector[(Eigen::Index) edge.m_var_full_e2] * scale_factor_e;
			fields.m_field_edges_x[i].m2 = (edge.m_var_full_m2 == INDEX_NONE)? 0.0 : eigenvector[(Eigen::Index) edge.m_var_full_m2] * scale_factor_m;
			fields.m_field_edges_x[i].mt1 = (edge.m_var_full_mt1 == INDEX_NONE)? 0.0 : eigenvector[(Eigen::Index) edge.m_var_full_mt1] * scale_factor_mt;
			fields.m_field_edges_x[i].mt2 = (edge.m_var_full_mt2 == INDEX_NONE)? 0.0 : eigenvector[(Eigen::Index) edge.m_var_full_mt2] * scale_factor_mt;
		}
		for(size_t i = 0; i < m_edges_y.size(); ++i) {
			Edge &edge = m_edges_y[i];
			fields.m_field_edges_y[i].e2 = (edge.m_var_full_e2 == INDEX_NONE)? 0.0 : eigenvector[(Eigen::Index) edge.m_var_full_e2] * scale_factor_e;
			fields.m_field_edges_y[i].m2 = (edge.m_var_full_m2 == INDEX_NONE)? 0.0 : eigenvector[(Eigen::Index) edge.m_var_full_m2] * scale_factor_m;
			fields.m_field_edges_y[i].mt1 = (edge.m_var_full_mt1 == INDEX_NONE)? 0.0 : eigenvector[(Eigen::Index) edge.m_var_full_mt1] * scale_factor_mt;
			fields.m_field_edges_y[i].mt2 = (edge.m_var_full_mt2 == INDEX_NONE)? 0.0 : eigenvector[(Eigen::Index) edge.m_var_full_mt2] * scale_factor_mt;
		}

	}

}

void GridMesh2D::GetCellField(size_t mode, size_t ix, size_t iy, real_t fx, real_t fy, complex_t &fex, complex_t &fey, complex_t &fez, complex_t &fmx, complex_t &fmy, complex_t &fmz) {
	SolutionField &fields = m_solution_fields[mode];
	real_t delta_x = m_grid_x[ix + 1] - m_grid_x[ix];
	real_t delta_y = m_grid_y[iy + 1] - m_grid_y[iy];
	complex_t s(0.0, 2.0 * M_PI * GetFrequency());
	complex_t en00 = fields.m_field_nodes[GetNodeIndex(ix    , iy    )].e1;
	complex_t en01 = fields.m_field_nodes[GetNodeIndex(ix + 1, iy    )].e1;
	complex_t en10 = fields.m_field_nodes[GetNodeIndex(ix    , iy + 1)].e1;
	complex_t en11 = fields.m_field_nodes[GetNodeIndex(ix + 1, iy + 1)].e1;
	complex_t eex0 = fields.m_field_edges_x[GetEdgeXIndex(ix    , iy    )].e2;
	complex_t eex1 = fields.m_field_edges_x[GetEdgeXIndex(ix    , iy + 1)].e2;
	complex_t eey0 = fields.m_field_edges_y[GetEdgeYIndex(ix    , iy    )].e2;
	complex_t eey1 = fields.m_field_edges_y[GetEdgeYIndex(ix + 1, iy    )].e2;
	complex_t mn00 = fields.m_field_nodes[GetNodeIndex(ix    , iy    )].m1;
	complex_t mn01 = fields.m_field_nodes[GetNodeIndex(ix + 1, iy    )].m1;
	complex_t mn10 = fields.m_field_nodes[GetNodeIndex(ix    , iy + 1)].m1;
	complex_t mn11 = fields.m_field_nodes[GetNodeIndex(ix + 1, iy + 1)].m1;
	complex_t mex0 = fields.m_field_edges_x[GetEdgeXIndex(ix    , iy    )].m2;
	complex_t mex1 = fields.m_field_edges_x[GetEdgeXIndex(ix    , iy + 1)].m2;
	complex_t mey0 = fields.m_field_edges_y[GetEdgeYIndex(ix    , iy    )].m2;
	complex_t mey1 = fields.m_field_edges_y[GetEdgeYIndex(ix + 1, iy    )].m2;
	complex_t mt1x0 = fields.m_field_edges_x[GetEdgeXIndex(ix    , iy    )].mt1;
	complex_t mt1x1 = fields.m_field_edges_x[GetEdgeXIndex(ix    , iy + 1)].mt1;
	complex_t mt1y0 = fields.m_field_edges_y[GetEdgeYIndex(ix    , iy    )].mt1;
	complex_t mt1y1 = fields.m_field_edges_y[GetEdgeYIndex(ix + 1, iy    )].mt1;
	complex_t mt2x0 = fields.m_field_edges_x[GetEdgeXIndex(ix    , iy    )].mt2;
	complex_t mt2x1 = fields.m_field_edges_x[GetEdgeXIndex(ix    , iy + 1)].mt2;
	complex_t mt2y0 = fields.m_field_edges_y[GetEdgeYIndex(ix    , iy    )].mt2;
	complex_t mt2y1 = fields.m_field_edges_y[GetEdgeYIndex(ix + 1, iy    )].mt2;
	fex = -Interp_Az_Dx(en00, en01, en10, en11, eex0, eex1, eey0, eey1, fx, fy) / delta_x - s * Interp_Ax(mt1x0, mt1x1, mt2x0, mt2x1, fx, fy);
	fey = -Interp_Az_Dy(en00, en01, en10, en11, eex0, eex1, eey0, eey1, fx, fy) / delta_y - s * Interp_Ay(mt1y0, mt1y1, mt2y0, mt2y1, fx, fy);
	fez = fields.m_propagation_constant * Interp_Az(en00, en01, en10, en11, eex0, eex1, eey0, eey1, fx, fy) - s * Interp_Az(mn00, mn01, mn10, mn11, mex0, mex1, mey0, mey1, fx, fy);
	fmx =  Interp_Az_Dy(mn00, mn01, mn10, mn11, mex0, mex1, mey0, mey1, fx, fy) / delta_y + fields.m_propagation_constant * Interp_Ay(mt1x0, mt1x1, mt2x0, mt2x1, fx, fy);
	fmy = -Interp_Az_Dx(mn00, mn01, mn10, mn11, mex0, mex1, mey0, mey1, fx, fy) / delta_x - fields.m_propagation_constant * Interp_Ax(mt1y0, mt1y1, mt2y0, mt2y1, fx, fy);
	fmz = Interp_Ay_Dx(mt1y0, mt1y1, mt2y0, mt2y1, fx, fy) / delta_x - Interp_Ax_Dy(mt1x0, mt1x1, mt2x0, mt2x1, fx, fy) / delta_y;
}

void GridMesh2D::GetCellValues(std::vector<real_t> &cell_values, size_t mode, MeshImageType type) {
	assert(IsInitialized());
	assert(mode < GetModeCount());
	assert(type == MESHIMAGETYPE_MESH);
	UNUSED(mode);
	UNUSED(type);
	size_t num_materials = 0;
	const MaterialDielectric *last_material = NULL;
	std::vector<real_t> dielectric_values(m_dielectrics.size());
	for(size_t i = 0; i < m_dielectrics.size(); ++i) {
		if(m_dielectrics[i].m_material != last_material) {
			++num_materials;
			last_material = m_dielectrics[i].m_material;
		}
		dielectric_values[i] = (real_t) num_materials;
	}
	for(size_t i = 0; i < m_dielectrics.size(); ++i) {
		dielectric_values[i] *= 0.9 / (real_t) (num_materials + 1);
	}
	cell_values.resize(m_cells.size());
	for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
		for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
			size_t cell_index = GetCellIndex(ix, iy);
			Cell &cell = m_cells[cell_index];
			real_t val = (cell.m_conductor != INDEX_NONE)? 0.0 : (cell.m_dielectric != INDEX_NONE)? dielectric_values[cell.m_dielectric] : 0.9;
			if((ix & 1) == (iy & 1))
				val += 0.1;
			cell_values[cell_index] = val;
		}
	}
}

void GridMesh2D::GetCellNodeValues(std::vector<std::array<real_t, 4>> &cellnode_values, size_t mode, MeshImageType type) {
	assert(IsInitialized() && IsSolved());
	assert(mode < GetModeCount());
	cellnode_values.clear();
	cellnode_values.resize(m_cells.size());
	switch(type) {
		case MESHIMAGETYPE_EFIELD: {
			break;
		}
		case MESHIMAGETYPE_MFIELD: {
			break;
		}
		case MESHIMAGETYPE_ENERGY: {
			complex_t *solution_epot = m_solution_static_e.data() + m_solution_static_e.outerStride() * (ptrdiff_t) mode;
			complex_t *solution_mpot = m_solution_static_m.data() + m_solution_static_m.outerStride() * (ptrdiff_t) mode;
			const real_t *fixed_values = GetModes().data() + GetModes().outerStride() * (ptrdiff_t) mode;
			std::vector<std::vector<real_t>> dielectric_values(m_dielectrics.size() + 1);
			std::vector<std::vector<int32_t>> dielectric_count(m_dielectrics.size() + 1);
			for(size_t i = 0; i < dielectric_values.size(); ++i) {
				dielectric_values[i].resize(m_nodes.size(), 0.0);
				dielectric_count[i].resize(m_nodes.size(), 0);
			}
			for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
				for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
					Cell &cell = GetCell(ix, iy);
					if(cell.m_conductor != INDEX_NONE)
						continue;
					size_t node00_index = GetNodeIndex(ix    , iy    );
					size_t node01_index = GetNodeIndex(ix + 1, iy    );
					size_t node10_index = GetNodeIndex(ix    , iy + 1);
					size_t node11_index = GetNodeIndex(ix + 1, iy + 1);
					Node &node00 = m_nodes[node00_index];
					Node &node01 = m_nodes[node01_index];
					Node &node10 = m_nodes[node10_index];
					Node &node11 = m_nodes[node11_index];
					//Node &node00 = GetNode(ix    , iy    );
					//Node &node01 = GetNode(ix + 1, iy    );
					//Node &node10 = GetNode(ix    , iy + 1);
					//Node &node11 = GetNode(ix + 1, iy + 1);
					complex_t node00_value_epot = (node00.m_var_static_e1 < INDEX_OFFSET)? solution_epot[node00.m_var_static_e1] : fixed_values[node00.m_var_static_e1 - INDEX_OFFSET];
					complex_t node01_value_epot = (node01.m_var_static_e1 < INDEX_OFFSET)? solution_epot[node01.m_var_static_e1] : fixed_values[node01.m_var_static_e1 - INDEX_OFFSET];
					complex_t node10_value_epot = (node10.m_var_static_e1 < INDEX_OFFSET)? solution_epot[node10.m_var_static_e1] : fixed_values[node10.m_var_static_e1 - INDEX_OFFSET];
					complex_t node11_value_epot = (node11.m_var_static_e1 < INDEX_OFFSET)? solution_epot[node11.m_var_static_e1] : fixed_values[node11.m_var_static_e1 - INDEX_OFFSET];
					complex_t node00_value_mpot = (node00.m_var_static_m1 < INDEX_OFFSET)? solution_mpot[node00.m_var_static_m1] : fixed_values[node00.m_var_static_m1 - INDEX_OFFSET];
					complex_t node01_value_mpot = (node01.m_var_static_m1 < INDEX_OFFSET)? solution_mpot[node01.m_var_static_m1] : fixed_values[node01.m_var_static_m1 - INDEX_OFFSET];
					complex_t node10_value_mpot = (node10.m_var_static_m1 < INDEX_OFFSET)? solution_mpot[node10.m_var_static_m1] : fixed_values[node10.m_var_static_m1 - INDEX_OFFSET];
					complex_t node11_value_mpot = (node11.m_var_static_m1 < INDEX_OFFSET)? solution_mpot[node11.m_var_static_m1] : fixed_values[node11.m_var_static_m1 - INDEX_OFFSET];
					complex_t ex0 = node01_value_epot - node00_value_epot, ex1 = node11_value_epot - node10_value_epot;
					complex_t ey0 = node10_value_epot - node00_value_epot, ey1 = node11_value_epot - node01_value_epot;
					complex_t hx0 = node01_value_mpot - node00_value_mpot, hx1 = node11_value_mpot - node10_value_mpot;
					complex_t hy0 = node10_value_mpot - node00_value_mpot, hy1 = node11_value_mpot - node01_value_mpot;
					real_t scale_x = 1.0 / square(m_grid_x[ix + 1] - m_grid_x[ix]);
					real_t scale_y = 1.0 / square(m_grid_y[iy + 1] - m_grid_y[iy]);
					size_t dielectric_index = (cell.m_dielectric == INDEX_NONE)? 0 : cell.m_dielectric + 1;
					dielectric_values[dielectric_index][node00_index] += (ex0 * std::conj(hx0)).real() * scale_x + (ey0 * std::conj(hy0)).real() * scale_y;
					dielectric_values[dielectric_index][node01_index] += (ex0 * std::conj(hx0)).real() * scale_x + (ey1 * std::conj(hy1)).real() * scale_y;
					dielectric_values[dielectric_index][node10_index] += (ex1 * std::conj(hx1)).real() * scale_x + (ey0 * std::conj(hy0)).real() * scale_y;
					dielectric_values[dielectric_index][node11_index] += (ex1 * std::conj(hx1)).real() * scale_x + (ey1 * std::conj(hy1)).real() * scale_y;
					++dielectric_count[dielectric_index][node00_index];
					++dielectric_count[dielectric_index][node01_index];
					++dielectric_count[dielectric_index][node10_index];
					++dielectric_count[dielectric_index][node11_index];
					//cellnode_values[cell_index][0] = ex0 * hx0 * scale_x + ey0 * hy0 * scale_y;
					//cellnode_values[cell_index][1] = ex0 * hx0 * scale_x + ey1 * hy1 * scale_y;
					//cellnode_values[cell_index][2] = ex1 * hx1 * scale_x + ey0 * hy0 * scale_y;
					//cellnode_values[cell_index][3] = ex1 * hx1 * scale_x + ey1 * hy1 * scale_y;
				}
			}
			for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
				for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
					size_t cell_index = GetCellIndex(ix, iy);
					Cell &cell = m_cells[cell_index];
					if(cell.m_conductor != INDEX_NONE)
						continue;
					size_t node00_index = GetNodeIndex(ix    , iy    );
					size_t node01_index = GetNodeIndex(ix + 1, iy    );
					size_t node10_index = GetNodeIndex(ix    , iy + 1);
					size_t node11_index = GetNodeIndex(ix + 1, iy + 1);
					size_t dielectric_index = (cell.m_dielectric == INDEX_NONE)? 0 : cell.m_dielectric + 1;
					cellnode_values[cell_index][0] = dielectric_values[dielectric_index][node00_index] / (real_t) dielectric_count[dielectric_index][node00_index];
					cellnode_values[cell_index][1] = dielectric_values[dielectric_index][node01_index] / (real_t) dielectric_count[dielectric_index][node01_index];
					cellnode_values[cell_index][2] = dielectric_values[dielectric_index][node10_index] / (real_t) dielectric_count[dielectric_index][node10_index];
					cellnode_values[cell_index][3] = dielectric_values[dielectric_index][node11_index] / (real_t) dielectric_count[dielectric_index][node11_index];
				}
			}
			real_t max_value = 0.0;
			for(size_t i = 0; i < cellnode_values.size(); ++i) {
				for(size_t j = 0; j < 4; ++j) {
					max_value = std::max(max_value, fabs(cellnode_values[i][j]));
				}
			}
			real_t scale = 1.0 / max_value;
			for(size_t i = 0; i < cellnode_values.size(); ++i) {
				for(size_t j = 0; j < 4; ++j) {
					cellnode_values[i][j] *= scale;
				}
			}
			break;
		}
		case MESHIMAGETYPE_CURRENT: {
			/*real_t *solution_values = m_solution_full_em.data() + m_solution_full_em.outerStride() * (ptrdiff_t) mode;
			real_t max_value = 0.0;
			for(size_t i = 0; i < m_vars_surf; ++i) {
				max_value = std::max(max_value, fabs(solution_values[i]));
			}
			real_t scale = 1.0 / max_value;
			for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
				for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
					size_t cell_index = GetCellIndex(ix, iy);
					if(m_cells[cell_index].m_conductor == INDEX_NONE)
						continue;
					Node &node00 = GetNode(ix    , iy    );
					Node &node01 = GetNode(ix + 1, iy    );
					Node &node10 = GetNode(ix    , iy + 1);
					Node &node11 = GetNode(ix + 1, iy + 1);
					cellnode_values[cell_index][0] = (node00.m_var_surf == INDEX_NONE)? 0.0 : solution_values[node00.m_var_surf] * scale;
					cellnode_values[cell_index][1] = (node01.m_var_surf == INDEX_NONE)? 0.0 : solution_values[node01.m_var_surf] * scale;
					cellnode_values[cell_index][2] = (node10.m_var_surf == INDEX_NONE)? 0.0 : solution_values[node10.m_var_surf] * scale;
					cellnode_values[cell_index][3] = (node11.m_var_surf == INDEX_NONE)? 0.0 : solution_values[node11.m_var_surf] * scale;
				}
			}*/
			// TODO: reimplement
			break;
		}
		default: assert(false);
	}
}

void GridMesh2D::GridAddBox(std::vector<GridLine> &grid_x, std::vector<GridLine> &grid_y, const Box2D &box, const Box2D &step) {
	grid_x.emplace_back(box.x1, step.x1);
	grid_x.emplace_back(box.x2, step.x2);
	grid_y.emplace_back(box.y1, step.y1);
	grid_y.emplace_back(box.y2, step.y2);
}

void GridMesh2D::GridRefine(std::vector<real_t> &result, std::vector<GridLine> &grid, const std::vector<real_t> &max_step, real_t inc, real_t epsilon) {
	if(grid.size() < 2)
		throw std::runtime_error("GridMesh2D error: The mesh must have at least 2 grid lines.");

	// propagate minimum step size
	real_t min_step = grid[0].m_step;
	for(size_t i = 1; i < grid.size(); ++i) {
		if(min_step != REAL_MAX)
			min_step += (grid[i].m_value - grid[i - 1].m_value) * inc;
		min_step = std::min(min_step, grid[i].m_step);
		grid[i].m_step = min_step;
	}
	for(size_t i = grid.size() - 1; i > 0; ) {
		--i;
		if(min_step != REAL_MAX)
			min_step += (grid[i + 1].m_value - grid[i].m_value) * inc;
		min_step = std::min(min_step, grid[i].m_step);
		grid[i].m_step = min_step;
	}

	// combine lines that are closer than 'epsilon'
	result.clear();
	real_t result_back_step = 0.0;
	real_t avg_sum = grid[0].m_value;
	real_t avg_prev = grid[0].m_value;
	real_t avg_step = grid[0].m_step;
	size_t avg_count = 1;
	for(size_t i = 1; i < grid.size(); ++i) {
		const GridLine &gridline = grid[i];
		if(gridline.m_value - avg_prev <= epsilon) {
			avg_sum += gridline.m_value;
			avg_prev = gridline.m_value;
			avg_step = std::min(avg_step, gridline.m_step);
			++avg_count;
		} else if(result.empty()) {
			result.push_back(avg_sum / (real_t) avg_count);
			result_back_step = avg_step;
			avg_sum = gridline.m_value;
			avg_prev = gridline.m_value;
			avg_step = gridline.m_step;
			avg_count = 1;
		} else {
			GridRefine2(result, result.back(), avg_sum / (real_t) avg_count, result_back_step, avg_step, max_step[i - 1], inc);
			result_back_step = avg_step;
			avg_sum = gridline.m_value;
			avg_prev = gridline.m_value;
			avg_step = gridline.m_step;
			avg_count = 1;
		}
	}

	// complete the grid
	if(result.empty())
		throw std::runtime_error("GridMesh2D error: The mesh must have at least 2 grid lines.");
	GridRefine2(result, result.back(), avg_sum / (real_t) avg_count, result_back_step, avg_step, max_step.back(), inc);
	if(result.size() < 2)
		throw std::runtime_error("GridMesh2D error: The mesh must have at least 2 grid lines.");

}

void GridMesh2D::GridRefine2(std::vector<real_t> &result, real_t x1, real_t x2, real_t step1, real_t step2, real_t max_step, real_t inc) {
	assert(FiniteLess(x1, x2));
	assert(FinitePositive(step1));
	assert(FinitePositive(step2));
	assert(FinitePositive(inc));

	// the goal is to ensure that the grid size near important grid lines is never larger than 'step'
	// (frac**(n-1) - frac**n) / (1 - frac**n) = st
	// frac**(n-1) - frac**n = st * (1 - frac**n)
	// frac**(n-1) - frac**n = st - st * frac**n
	// frac**(n-1) + (st - 1) * frac**n = st
	// frac**(n-1) = st / (1 + (st - 1) * frac)
	// n = log2(st / (1 + (st - 1) * frac)) / log2(frac) + 1

	if(x2 - x1 < std::min(step1, step2)) {
		result.push_back(x2);
		return;
	}

	// split into parts
	real_t frac = 1.0 / (1.0 + inc);
	real_t log2frac = log2(frac);
	real_t delta = x2 - x1;
	real_t split1 = clamp((max_step - step1) / inc, 0.0, delta);
	real_t split2 = clamp(delta - (max_step - step2) / inc, 0.0, delta);
	if(split2 - split1 < max_step * 0.5) {
		real_t split = ((step2 - step1) / inc + delta) * 0.5;
		if(split < step1 * 0.5) {
			split = 0.0;
		} else if(split > delta - step2 * 0.5) {
			split = delta;
		}
		split1 = split;
		split2 = split;
	} else {
		if(split1 < step1 * 0.5) {
			split1 = 0.0;
		}
		if(split2 > delta - step2 * 0.5) {
			split2 = delta;
		}
	}
	real_t delta1 = split1;
	real_t delta2 = split2 - split1;
	real_t delta3 = delta - split2;

	// part 1
	if(delta1 != 0.0) {
		real_t st = step1 / delta1;
		size_t n = (size_t) std::max<ptrdiff_t>(1, rintp(log2(st / (1.0 + (st - 1.0) * frac)) / log2frac + 1.5));
		real_t vmin = exp2(log2frac * (real_t) n);
		for(size_t i = n - 1; i > 0; --i) {
			result.push_back(x1 + delta1 * (exp2(log2frac * (real_t) i) - vmin) / (1.0 - vmin));
		}
		result.push_back(x1 + delta1);
	}

	// part 2
	if(delta2 != 0.0) {
		size_t n = (size_t) std::max<ptrdiff_t>(1, rintp(delta2 / max_step + 0.5));
		for(size_t i = 1; i <= n; ++i) {
			result.push_back(x1 + delta1 + delta2 * ((real_t) i) / ((real_t) n));
		}
	}

	// part 3
	if(delta3 != 0.0) {
		real_t st = step2 / delta3;
		size_t n = (size_t) std::max<ptrdiff_t>(1, rintp(log2(st / (1.0 + (st - 1.0) * frac)) / log2frac + 1.5));
		real_t vmin = exp2(log2frac * (real_t) n);
		for(size_t i = 1; i < n; ++i) {
			result.push_back(x2 - delta3 * (exp2(log2frac * (real_t) i) - vmin) / (1.0 - vmin));
		}
		result.push_back(x2);
	}

}

void GridMesh2D::GridMidpoints(std::vector<real_t> &result, const std::vector<real_t> &grid) {
	assert(grid.size() >= 2);
	result.clear();
	result.resize(grid.size() - 1);
	for(size_t i = 0; i < result.size(); ++i) {
		result[i] = (grid[i] + grid[i + 1]) * 0.5;
	}
}

void GridMesh2D::GridMidpoints(std::vector<real_t> &result, const std::vector<GridLine> &grid) {
	assert(grid.size() >= 2);
	result.clear();
	result.resize(grid.size() - 1);
	for(size_t i = 0; i < result.size(); ++i) {
		result[i] = (grid[i].m_value + grid[i + 1].m_value) * 0.5;
	}
}

void GridMesh2D::PrepareNodeImage(std::vector<size_t> &index, std::vector<real_t> &frac, const std::vector<real_t> &grid, real_t value1, real_t value2, size_t size) {
	index.resize(size);
	frac.resize(size);
	auto range_begin = grid.begin() + 1, range_end = grid.end() - 1;
	for(size_t i = 0; i < size; ++i) {
		real_t value = value1 + (value2 - value1) * ((real_t) i + 0.5) / (real_t) size;
		size_t j = (size_t) (std::upper_bound(range_begin, range_end, value) - range_begin);
		index[i] = j;
		frac[i] = (value - grid[j]) / (grid[j + 1] - grid[j]);
	}
}

void GridMesh2D::PrepareCellImage(std::vector<size_t> &index, std::vector<real_t> &frac, const std::vector<real_t> &grid, const std::vector<real_t> &midpoints, real_t value1, real_t value2, size_t size) {
	index.resize(size);
	frac.resize(size);
	real_t scale = (real_t) size / fabs(value2 - value1);
	auto range_begin = midpoints.begin() + 1, range_end = midpoints.end() - 1;
	for(size_t i = 0; i < size; ++i) {
		real_t value = value1 + (value2 - value1) * ((real_t) i + 0.5) / (real_t) size;
		size_t j = (size_t) (std::upper_bound(range_begin, range_end, value) - range_begin);
		index[i] = j;
		frac[i] = clamp((value - grid[j + 1]) * scale, -1.5, 1.5);
	}
}
