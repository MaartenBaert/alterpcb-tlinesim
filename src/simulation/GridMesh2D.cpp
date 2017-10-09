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

#include "MiscMath.h"
#include "StringHelper.h"

#include <cfloat>

#include <algorithm>

#define SIMULATION_VERBOSE 0

#if SIMULATION_VERBOSE
#include <chrono>
#include <iostream>
#endif

GridMesh2D::GridMesh2D(const Box2D &world_box, const Box2D &world_focus, real_t grid_inc, real_t grid_epsilon) {
	if(!FinitePositive(grid_inc))
		throw std::runtime_error("GridMesh2D error: grid_inc must be positive.");
	if(!FinitePositive(grid_epsilon))
		throw std::runtime_error("GridMesh2D error: grid_epsilon must be positive.");
	m_world_box = world_box.Normalized();
	m_world_focus = world_focus.Normalized();
	m_grid_inc = grid_inc;
	m_grid_epsilon = grid_epsilon;
	m_initialized = false;
	m_vars_real = 0;
	m_vars_fixed = 0;
	m_solved = false;
	m_solution_frequency = 0.0;
	m_mode_count = 0;
}

GridMesh2D::~GridMesh2D() {
	// nothing
}

size_t GridMesh2D::AddPort(GridMesh2D::PortType type) {
	if(m_initialized)
		throw std::runtime_error("GridMesh2D error: Can't add port after initialization.");
	m_ports.emplace_back(type);
	return m_ports.size() - 1;
}

void GridMesh2D::AddConductor(const Box2D &box, real_t step, const MaterialConductor *material, size_t port) {
	AddConductor(box, step, step, step, step, material, port);
}

void GridMesh2D::AddConductor(const Box2D &box, real_t step_left, real_t step_right, real_t step_top, real_t step_bottom, const MaterialConductor *material, size_t port) {
	if(m_initialized)
		throw std::runtime_error("GridMesh2D error: Can't add conductor after initialization.");
	if(!std::isfinite(box.x1) || !std::isfinite(box.y1) || !std::isfinite(box.x2) || !std::isfinite(box.y2))
		throw std::runtime_error("GridMesh2D error: Conductor box must be finite.");
	if(!FinitePositive(step_left) || !FinitePositive(step_right) || !FinitePositive(step_top) || !FinitePositive(step_bottom))
		throw std::runtime_error("GridMesh2D error: Conductor step must be positive.");
	m_conductors.emplace_back(box.Normalized(), step_left, step_right, step_top, step_bottom, material, port);
}

void GridMesh2D::AddDielectric(const Box2D &box, real_t step, const MaterialDielectric *material) {
	AddDielectric(box, step, step, step, step, material);
}

void GridMesh2D::AddDielectric(const Box2D &box, real_t step_left, real_t step_right, real_t step_top, real_t step_bottom, const MaterialDielectric *material) {
	if(m_initialized)
		throw std::runtime_error("GridMesh2D error: Can't add dielectric after initialization.");
	if(!std::isfinite(box.x1) || !std::isfinite(box.y1) || !std::isfinite(box.x2) || !std::isfinite(box.y2))
		throw std::runtime_error("GridMesh2D error: Dielectric box must be finite.");
	if(!FinitePositive(step_left) || !FinitePositive(step_right) || !FinitePositive(step_top) || !FinitePositive(step_bottom))
		throw std::runtime_error("GridMesh2D error: Dielectric step must be positive.");
	m_dielectrics.emplace_back(box.Normalized(), step_left, step_right, step_top, step_bottom, material);
}

void GridMesh2D::Initialize() {
	if(m_initialized)
		throw std::runtime_error("GridMesh2D error: The mesh has already been initialized.");

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
	std::cerr << "GridMesh2D stats:"
			  << " grid=" << m_grid_x.size() << "x" << m_grid_y.size()
			  << " nodes=" << m_nodes.size()
			  << " cells=" << m_cells.size()
			  << " vars=" << m_vars_real << "+" << m_vars_fixed
			  << std::endl;
	std::cerr << "GridMesh2D init time:"
			  << " grid=" << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us"
			  << " cells=" << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << "us"
			  << " vars=" << std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3).count() << "us"
			  << std::endl;
#endif

	m_initialized = true;

}

void GridMesh2D::Solve(std::vector<real_t> &charges, std::vector<real_t> &currents, const std::vector<real_t> &modes, size_t mode_count, real_t frequency) {
	if(!m_initialized)
		throw std::runtime_error("GridMesh2D error: The mesh must be initialized first.");
	if(modes.size() != m_vars_fixed * mode_count)
		throw std::runtime_error(MakeString("GridMesh2D error: Expected modes array with ", m_vars_fixed * mode_count, " elements, got ", modes.size(), " instead."));

	m_solved = false;
	m_solution_frequency = frequency;
	m_mode_count = mode_count;
	m_modes = modes;

#if SIMULATION_VERBOSE
	auto t1 = std::chrono::high_resolution_clock::now();
#endif
	GenerateMatrices();
#if SIMULATION_VERBOSE
	auto t2 = std::chrono::high_resolution_clock::now();
#endif
	SolveModes(charges, currents);
#if SIMULATION_VERBOSE
	auto t3 = std::chrono::high_resolution_clock::now();
#endif

#if SIMULATION_VERBOSE
	std::cerr << "GridMesh2D solve time:"
			  << " build=" << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() << "ms"
			  << " solve=" << std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() << "ms"
			  << std::endl;
#endif

	m_solved = true;

}

void GridMesh2D::Cleanup() {
	m_cholmod_sparse.Reset();
	m_cholmod_factor.Reset();
	m_cholmod_rhs.Reset();
	m_cholmod_ws1.Reset();
	m_cholmod_ws2.Reset();
}

Box2D GridMesh2D::GetWorldBox2D() {
	return m_world_box;
}

bool GridMesh2D::IsInitialized() {
	return m_initialized;
}

bool GridMesh2D::IsSolved() {
	return m_solved;
}

Box2D GridMesh2D::GetWorldFocus2D() {
	return m_world_focus;
}

size_t GridMesh2D::GetModeCount() {
	if(!m_initialized || !m_solved)
		return 0;
	return m_mode_count;
}

bool GridMesh2D::GetImage2D(std::vector<real_t> &image_value, std::vector<Vector2D> &image_gradient, size_t width, size_t height, const Box2D &view, MeshImageType type, size_t mode) {
	if(!m_initialized)
		return false;
	if(type == MESHIMAGETYPE_FIELD_E || type == MESHIMAGETYPE_FIELD_H) {
		if(!m_solved || mode >= m_mode_count)
			return false;
	}

	// clear image data
	image_value.clear();
	image_value.resize(width * height);
	if(type == MESHIMAGETYPE_FIELD_E || type == MESHIMAGETYPE_FIELD_H) {
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
					real_t v0 = v00 + (v01 - v00) * frac_x[i];
					real_t v1 = v10 + (v11 - v10) * frac_x[i];
					row_value[i] = v0 + (v1 - v0) * frac_y[j];
				}
			}
			break;
		}
		case MESHIMAGETYPE_FIELD_E:
		case MESHIMAGETYPE_FIELD_H: {
			std::vector<real_t> node_values;
			GetNodeValues(node_values, mode, type);
			for(size_t j = 0; j < height; ++j) {
				real_t *row_value = image_value.data() + j * width;
				Vector2D *row_gradient = image_gradient.data() + j * width;
				for(size_t i = 0; i < width; ++i) {
					size_t ix = index_x[i], iy = index_y[j];
					real_t v00 = node_values[GetNodeIndex(ix    , iy    )];
					real_t v01 = node_values[GetNodeIndex(ix + 1, iy    )];
					real_t v10 = node_values[GetNodeIndex(ix    , iy + 1)];
					real_t v11 = node_values[GetNodeIndex(ix + 1, iy + 1)];
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

	return true;
}

void GridMesh2D::InitGrid() {
	assert(!m_initialized);

	// add grid lines
	std::vector<GridLine> grid_x, grid_y;
	GridAddBox(grid_x, grid_y, m_world_box, REAL_MAX, REAL_MAX, REAL_MAX, REAL_MAX);
	for(Conductor &conductor : m_conductors) {
		GridAddBox(grid_x, grid_y, conductor.m_box, conductor.m_step_left, conductor.m_step_right, conductor.m_step_top, conductor.m_step_bottom);
	}
	for(Dielectric &dielectric : m_dielectrics) {
		GridAddBox(grid_x, grid_y, dielectric.m_box, dielectric.m_step_left, dielectric.m_step_right, dielectric.m_step_top, dielectric.m_step_bottom);
	}

	// refine grid
	GridRefine(m_grid_x, grid_x, m_grid_inc, m_grid_epsilon);
	GridRefine(m_grid_y, grid_y, m_grid_inc, m_grid_epsilon);

	if(m_grid_x.size() < 2 || m_grid_y.size() < 2)
		throw std::runtime_error("GridMesh2D error: The mesh must have at least 2 grid lines.");

	/*std::cerr << "m_grid_x =";
	for(real_t val : m_grid_x) {
		std::cerr << " " << val;
	}
	std::cerr << std::endl;
	std::cerr << "m_grid_y =";
	for(real_t val : m_grid_y) {
		std::cerr << " " << val;
	}
	std::cerr << std::endl;*/

	// calculate midpoints
	GridMidpoints(m_midpoints_x, m_grid_x);
	GridMidpoints(m_midpoints_y, m_grid_y);

}

void GridMesh2D::InitCells() {
	assert(!m_initialized);

	// generate nodes and cells
	m_nodes.clear();
	m_nodes.resize(m_grid_x.size() * m_grid_y.size());
	m_cells.clear();
	m_cells.resize((m_grid_x.size() - 1) * (m_grid_y.size() - 1));

	// apply conductors
	for(size_t conductor_index = 0; conductor_index < m_conductors.size(); ++conductor_index) {
		Conductor &conductor = m_conductors[conductor_index];
		size_t ix1 = std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), conductor.m_box.x1) - m_midpoints_x.begin();
		size_t ix2 = std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), conductor.m_box.x2) - m_midpoints_x.begin();
		size_t iy1 = std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), conductor.m_box.y1) - m_midpoints_y.begin();
		size_t iy2 = std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), conductor.m_box.y2) - m_midpoints_y.begin();
		for(size_t iy = iy1; iy <= iy2; ++iy) {
			for(size_t ix = ix1; ix <= ix2; ++ix) {
				Node &node = GetNode(ix, iy);
				if(node.m_port != INDEX_NONE && node.m_port != conductor.m_port) {
					throw std::runtime_error(MakeString("Port ", conductor.m_port, " makes contact with port ", node.m_port, "."));
				}
				node.m_port = conductor.m_port;
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
		size_t ix1 = std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), dielectric.m_box.x1) - m_midpoints_x.begin();
		size_t ix2 = std::upper_bound(m_midpoints_x.begin(), m_midpoints_x.end(), dielectric.m_box.x2) - m_midpoints_x.begin();
		size_t iy1 = std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), dielectric.m_box.y1) - m_midpoints_y.begin();
		size_t iy2 = std::upper_bound(m_midpoints_y.begin(), m_midpoints_y.end(), dielectric.m_box.y2) - m_midpoints_y.begin();
		for(size_t iy = iy1; iy < iy2; ++iy) {
			for(size_t ix = ix1; ix < ix2; ++ix) {
				Cell &cell = GetCell(ix, iy);
				cell.m_dielectric = dielectric_index;
			}
		}
	}

}

void GridMesh2D::InitVariables() {
	assert(!m_initialized);

	// assign variables to ports
	for(size_t i = 0; i < m_ports.size(); ++i) {
		Port &port = m_ports[i];
		switch(port.m_type) {
			case PORTTYPE_FIXED: {
				port.m_var = INDEX_OFFSET + m_vars_fixed++;
				break;
			}
			case PORTTYPE_FLOATING: {
				port.m_var = m_vars_real++;
				break;
			}
		}
	}

	// assign variables to nodes
	for(size_t iy = 0; iy < m_grid_y.size(); ++iy) {
		for(size_t ix = 0; ix < m_grid_x.size(); ++ix) {
			Node &node = GetNode(ix, iy);
			if(node.m_port == INDEX_NONE) {
				node.m_var = m_vars_real++;
			} else {
				node.m_var = m_ports[node.m_port].m_var;
			}
		}
	}

	if(m_vars_real == 0)
		throw std::runtime_error("GridMesh2D error: The mesh has no free variables.");
	if(m_vars_fixed == 0)
		throw std::runtime_error("GridMesh2D error: The mesh has no fixed variables.");

}

void GridMesh2D::GenerateMatrices() {
	assert(m_initialized && !m_solved);

	// load conductor properties
	std::vector<MaterialConductorProperties> conductor_properties(m_conductors.size());
	for(size_t i = 0; i < m_conductors.size(); ++i) {
		GetConductorProperties(conductor_properties[i], m_conductors[i].m_material, m_solution_frequency);
	}

	// load dielectric properties
	std::vector<MaterialDielectricProperties> dielectric_properties(m_dielectrics.size());
	for(size_t i = 0; i < m_dielectrics.size(); ++i) {
		GetDielectricProperties(dielectric_properties[i], m_dielectrics[i].m_material, m_solution_frequency);
	}

	// generate matrices
	m_matrix_e.Reset(m_vars_real, m_vars_fixed, 5);
	m_matrix_h.Reset(m_vars_real, m_vars_fixed, 5);
	for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
		for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {

			// get cell, skip cells that are inside conductors
			Cell &cell = GetCell(ix, iy);
			if(cell.m_conductor != INDEX_NONE)
				continue;

			// get neighboring nodes
			Node &node00 = GetNode(ix    , iy    );
			Node &node01 = GetNode(ix + 1, iy    );
			Node &node10 = GetNode(ix    , iy + 1);
			Node &node11 = GetNode(ix + 1, iy + 1);

			// get dielectric properties
			std::complex<real_t> permittivity_x(1.0, 0.0), permittivity_y(1.0, 0.0);
			if(cell.m_dielectric != INDEX_NONE) {
				permittivity_x = dielectric_properties[cell.m_dielectric].m_permittivity_x;
				permittivity_y = dielectric_properties[cell.m_dielectric].m_permittivity_y;
			}

			// calculate scale factors
			// TODO: use complex permittivity
			real_t delta_x = m_grid_x[ix + 1] - m_grid_x[ix];
			real_t delta_y = m_grid_y[iy + 1] - m_grid_y[iy];
			real_t scale_x_e = 1.0 / 6.0 * permittivity_x.real() * delta_y / delta_x;
			real_t scale_y_e = 1.0 / 6.0 * permittivity_y.real() * delta_x / delta_y;
			real_t scale_x_h = 1.0 / 6.0 * delta_y / delta_x;
			real_t scale_y_h = 1.0 / 6.0 * delta_x / delta_y;

			// calculate E coefficients
			real_t coef_x_e = -2.0 * scale_x_e + scale_y_e;
			real_t coef_y_e = -2.0 * scale_y_e + scale_x_e;
			real_t coef_d_e = -scale_x_e - scale_y_e;
			real_t coef_s_e = scale_x_e + scale_y_e;

			// calculate H coefficients
			real_t coef_x_h = -2.0 * scale_x_h + scale_y_h;
			real_t coef_y_h = -2.0 * scale_y_h + scale_x_h;
			real_t coef_d_h = -scale_x_h - scale_y_h;
			real_t coef_s_h = scale_x_h + scale_y_h;

			// add to E matrix
			m_matrix_e.Add(node00.m_var, node01.m_var, coef_x_e);
			m_matrix_e.Add(node10.m_var, node11.m_var, coef_x_e);
			m_matrix_e.Add(node00.m_var, node10.m_var, coef_y_e);
			m_matrix_e.Add(node01.m_var, node11.m_var, coef_y_e);
			m_matrix_e.Add(node00.m_var, node11.m_var, coef_d_e);
			m_matrix_e.Add(node01.m_var, node10.m_var, coef_d_e);
			m_matrix_e.Add(node00.m_var, node00.m_var, coef_s_e);
			m_matrix_e.Add(node01.m_var, node01.m_var, coef_s_e);
			m_matrix_e.Add(node10.m_var, node10.m_var, coef_s_e);
			m_matrix_e.Add(node11.m_var, node11.m_var, coef_s_e);

			// add to H matrix
			m_matrix_h.Add(node00.m_var, node01.m_var, coef_x_h);
			m_matrix_h.Add(node10.m_var, node11.m_var, coef_x_h);
			m_matrix_h.Add(node00.m_var, node10.m_var, coef_y_h);
			m_matrix_h.Add(node01.m_var, node11.m_var, coef_y_h);
			m_matrix_h.Add(node00.m_var, node11.m_var, coef_d_h);
			m_matrix_h.Add(node01.m_var, node10.m_var, coef_d_h);
			m_matrix_h.Add(node00.m_var, node00.m_var, coef_s_h);
			m_matrix_h.Add(node01.m_var, node01.m_var, coef_s_h);
			m_matrix_h.Add(node10.m_var, node10.m_var, coef_s_h);
			m_matrix_h.Add(node11.m_var, node11.m_var, coef_s_h);

		}
	}

}

void GridMesh2D::SolveModes(std::vector<real_t> &charges, std::vector<real_t> &currents) {
	assert(m_initialized && !m_solved);

	// preallocate matrix for right-hand-side
	m_cholmod_rhs.ResetReal(m_vars_real, m_mode_count);

	// factorize E matrix
	m_cholmod_sparse.Reset(m_matrix_e);
	m_cholmod_factor.Factorize(m_cholmod_sparse);

	// solve E matrix
	std::fill_n(m_cholmod_rhs.GetRealData(), m_cholmod_rhs.GetStride() * m_mode_count, 0.0);
	for(size_t i = 0; i < m_mode_count; ++i) {
		real_t *rhs = m_cholmod_rhs.GetRealData() + m_cholmod_rhs.GetStride() * i;
		const real_t *fixed_values = m_modes.data() + m_vars_fixed * i;
		m_matrix_e.SubtractFromRhs(rhs, fixed_values);
	}
	m_cholmod_factor.Solve(m_cholmod_solution_e, m_cholmod_rhs, m_cholmod_ws1, m_cholmod_ws2);

	// factorize H matrix
	m_cholmod_sparse.Reset(m_matrix_h);
	m_cholmod_factor.Factorize(m_cholmod_sparse);

	// solve H matrix
	std::fill_n(m_cholmod_rhs.GetRealData(), m_cholmod_rhs.GetStride() * m_mode_count, 0.0);
	for(size_t i = 0; i < m_mode_count; ++i) {
		m_matrix_h.SubtractFromRhs(m_cholmod_rhs.GetRealData() + m_cholmod_rhs.GetStride() * i, m_modes.data() + m_vars_fixed * i);
	}
	m_cholmod_factor.Solve(m_cholmod_solution_h, m_cholmod_rhs, m_cholmod_ws1, m_cholmod_ws2);

	// calculate the residual charges
	charges.resize(m_vars_fixed * m_mode_count);
	for(size_t i = 0; i < m_mode_count; ++i) {
		real_t *residual = charges.data() + m_vars_fixed * i;
		real_t *solution = m_cholmod_solution_e.GetRealData() + m_cholmod_solution_e.GetStride() * i;
		real_t *fixed_values = m_modes.data() + m_vars_fixed * i;
		m_matrix_e.CalculateResidual(residual, solution, fixed_values);
	}

	// calculate the residual currents
	currents.resize(m_vars_fixed * m_mode_count);
	for(size_t i = 0; i < m_mode_count; ++i) {
		real_t *residual = currents.data() + m_vars_fixed * i;
		real_t *solution = m_cholmod_solution_h.GetRealData() + m_cholmod_solution_h.GetStride() * i;
		real_t *fixed_values = m_modes.data() + m_vars_fixed * i;
		m_matrix_h.CalculateResidual(residual, solution, fixed_values);
	}

}

void GridMesh2D::GetCellValues(std::vector<real_t> &cell_values, size_t mode, MeshImageType type) {
	assert(m_initialized);
	assert(mode < m_mode_count);
	assert(type == MESHIMAGETYPE_MESH);
	UNUSED(mode);
	UNUSED(type);
	cell_values.resize(m_cells.size());
	for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
		for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
			size_t cell_index = GetCellIndex(ix, iy);
			Cell &cell = m_cells[cell_index];
			real_t val = (cell.m_conductor != INDEX_NONE)? 0.0 : (cell.m_dielectric != INDEX_NONE)? 0.45 : 0.9;
			if((ix & 1) == (iy & 1))
				val += 0.1;
			cell_values[cell_index] = val;
		}
	}
}

void GridMesh2D::GetNodeValues(std::vector<real_t> &node_values, size_t mode, MeshImageType type) {
	assert(m_initialized && m_solved);
	assert(mode < m_mode_count);
	assert(type == MESHIMAGETYPE_FIELD_E || type == MESHIMAGETYPE_FIELD_H);
	node_values.resize(m_nodes.size());
	CholmodDenseMatrix &solution = (type == MESHIMAGETYPE_FIELD_E)? m_cholmod_solution_e : m_cholmod_solution_h;
	real_t *solution_values = solution.GetRealData() + solution.GetStride() * mode;
	real_t *fixed_values = m_modes.data() + m_vars_fixed * mode;
	for(size_t i = 0; i < m_nodes.size(); ++i) {
		size_t var = m_nodes[i].m_var;
		node_values[i] = (var < INDEX_OFFSET)? solution_values[var] : fixed_values[var - INDEX_OFFSET];
	}
}

void GridMesh2D::GridAddBox(std::vector<GridLine> &grid_x, std::vector<GridLine> &grid_y, const Box2D &box, real_t step_left, real_t step_right, real_t step_top, real_t step_bottom) {
	grid_x.emplace_back(box.x1, step_left);
	grid_x.emplace_back(box.x2, step_right);
	grid_y.emplace_back(box.y1, step_top);
	grid_y.emplace_back(box.y2, step_bottom);
}

void GridMesh2D::GridRefine(std::vector<real_t> &result, std::vector<GridLine> &grid, real_t inc, real_t epsilon) {
	if(grid.size() < 2)
		throw std::runtime_error("GridMesh2D error: The mesh must have at least 2 grid lines.");

	// sort the list
	std::sort(grid.begin(), grid.end());

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
			avg_step = fmin(avg_step, gridline.m_step);
			++avg_count;
		} else if(result.empty()) {
			result.push_back(avg_sum / (real_t) avg_count);
			result_back_step = avg_step;
			avg_sum = gridline.m_value;
			avg_prev = gridline.m_value;
			avg_step = gridline.m_step;
			avg_count = 1;
		} else {
			GridRefine2(result, result.back(), avg_sum / (real_t) avg_count, result_back_step, avg_step, inc);
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
	GridRefine2(result, result.back(), avg_sum / (real_t) avg_count, result_back_step, avg_step, inc);
	if(result.size() < 2)
		throw std::runtime_error("GridMesh2D error: The mesh must have at least 2 grid lines.");

}

void GridMesh2D::GridRefine2(std::vector<real_t> &result, real_t x1, real_t x2, real_t step1, real_t step2, real_t inc) {
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

	if(x2 - x1 < fmin(step1, step2)) {
		result.push_back(x2);
		return;
	}

	// split into two halves
	real_t frac = 1.0 / (1.0 + inc);
	real_t log2frac = log2(frac);
	real_t split = (step2 - step1) / (inc * (x2 - x1)) * 0.5 + 0.5;
	real_t delta1 = (x2 - x1) * split;
	real_t delta2 = (x2 - x1) * (1.0 - split);

	// eliminate one half if empty
	if(delta1 < step1 * 0.5) {
		delta1 = 0.0;
		delta2 = x2 - x1;
	} else if(delta2 < step2 * 0.5) {
		delta1 = x2 - x1;
		delta2 = 0.0;
	}

	// first half
	if(delta1 != 0.0) {
		real_t st1 = step1 / delta1;
		size_t n1 = std::max<ptrdiff_t>(1, rints(log2(st1 / (1.0 + (st1 - 1.0) * frac)) / log2frac + 1.5));
		real_t vmin1 = exp2(log2frac * (real_t) n1);
		for(size_t i = n1 - 1; i > 0; --i) {
			result.push_back(x1 + delta1 * (exp2(log2frac * (real_t) i) - vmin1) / (1.0 - vmin1));
		}
		result.push_back(x1 + delta1);
	}

	// second half
	if(delta2 != 0.0) {
		real_t st2 = step2 / delta2;
		size_t n2 = std::max<ptrdiff_t>(1, rints(log2(st2 / (1.0 + (st2 - 1.0) * frac)) / log2frac + 1.5));
		real_t vmin2 = exp2(log2frac * (real_t) n2);
		for(size_t i = 1; i < n2; ++i) {
			result.push_back(x2 - delta2 * (exp2(log2frac * (real_t) i) - vmin2) / (1.0 - vmin2));
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

void GridMesh2D::PrepareNodeImage(std::vector<size_t> &index, std::vector<real_t> &frac, const std::vector<real_t> &grid, real_t value1, real_t value2, size_t size) {
	index.resize(size);
	frac.resize(size);
	auto range_begin = grid.begin() + 1, range_end = grid.end() - 1;
	for(size_t i = 0; i < size; ++i) {
		real_t value = value1 + (value2 - value1) * ((real_t) i + 0.5) / (real_t) size;
		size_t j = std::upper_bound(range_begin, range_end, value) - range_begin;
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
		size_t j = std::upper_bound(range_begin, range_end, value) - range_begin;
		index[i] = j;
		frac[i] = clamp((value - grid[j + 1]) * scale + 0.5, 0.0, 1.0);
	}
}
