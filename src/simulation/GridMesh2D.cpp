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
	m_vars_surf = 0;
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
	if(material == NULL)
		throw std::runtime_error("GridMesh2D error: Material can't be NULL.");
	if(port >= m_ports.size())
		throw std::runtime_error("GridMesh2D error: Invalid port index.");
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
	if(material == NULL)
		throw std::runtime_error("GridMesh2D error: Material can't be NULL.");
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
	InitSurfaceMatrix();
#if SIMULATION_VERBOSE
	auto t5 = std::chrono::high_resolution_clock::now();
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
			  << " surf=" << std::chrono::duration_cast<std::chrono::microseconds>(t5 - t4).count() << "us"
			  << std::endl;
#endif

	m_initialized = true;

}

void GridMesh2D::Solve(std::vector<real_t> &charges, std::vector<real_t> &currents, std::vector<real_t> &dielectric_losses,
					   std::vector<real_t> &resistive_losses, const std::vector<real_t> &modes, size_t mode_count, real_t frequency) {
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
	SolveModes(charges, currents, dielectric_losses, resistive_losses);
#if SIMULATION_VERBOSE
	auto t3 = std::chrono::high_resolution_clock::now();
#endif

#if SIMULATION_VERBOSE
	std::cerr << "GridMesh2D solve time:"
			  << " build=" << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us"
			  << " solve=" << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << "us"
			  << std::endl;
#endif

	m_solved = true;

}

void GridMesh2D::Cleanup() {
	m_cholmod_sparse.Reset();
	m_cholmod_sparse_surf.Reset();
	m_cholmod_factor.Reset();
	m_cholmod_factor_surf.Reset();
	m_cholmod_rhs.Reset();
	m_cholmod_rhs_surf.Reset();
	m_cholmod_ws1.Reset();
	m_cholmod_ws1_surf.Reset();
	m_cholmod_ws2.Reset();
	m_cholmod_ws2_surf.Reset();
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
	if(type == MESHIMAGETYPE_FIELD_E || type == MESHIMAGETYPE_FIELD_H || type == MESHIMAGETYPE_CURRENT) {
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
		case MESHIMAGETYPE_CURRENT: {
			std::vector<real_t> node_values;
			GetNodeValues(node_values, mode, type);
			for(size_t j = 0; j < height; ++j) {
				real_t *row_value = image_value.data() + j * width;
				for(size_t i = 0; i < width; ++i) {
					size_t ix = index_x[i], iy = index_y[j];
					real_t v00 = node_values[GetNodeIndex(ix    , iy    )];
					real_t v01 = node_values[GetNodeIndex(ix + 1, iy    )];
					real_t v10 = node_values[GetNodeIndex(ix    , iy + 1)];
					real_t v11 = node_values[GetNodeIndex(ix + 1, iy + 1)];
					real_t v0 = v00 + (v01 - v00) * frac_x[i];
					real_t v1 = v10 + (v11 - v10) * frac_x[i];
					row_value[i] = v0 + (v1 - v0) * frac_y[j];
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
	m_nodes.resize(m_grid_x.size() * m_grid_y.size());
	m_edges_h.resize((m_grid_x.size() - 1) * m_grid_y.size());
	m_edges_v.resize(m_grid_x.size() * (m_grid_y.size() - 1));
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
		for(size_t iy = iy1; iy <= iy2; ++iy) {
			for(size_t ix = ix1; ix < ix2; ++ix) {
				Edge &edge = GetEdgeH(ix, iy);
				edge.m_conductor = conductor_index;
			}
		}
		for(size_t iy = iy1; iy < iy2; ++iy) {
			for(size_t ix = ix1; ix <= ix2; ++ix) {
				Edge &edge = GetEdgeV(ix, iy);
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

	// assign variables to surface nodes
	// TODO: this is incorrect for planar conductors that are not located on the edge of the world box
	for(size_t iy = 0; iy < m_grid_y.size() - 1; ++iy) {
		for(size_t ix = 0; ix < m_grid_x.size() - 1; ++ix) {
			Cell &cell = GetCell(ix, iy);
			if(cell.m_conductor != INDEX_NONE)
				continue;
			Node &node00 = GetNode(ix    , iy    );
			Node &node01 = GetNode(ix + 1, iy    );
			Node &node10 = GetNode(ix    , iy + 1);
			Node &node11 = GetNode(ix + 1, iy + 1);
			if(node00.m_port != INDEX_NONE && node00.m_var_surf == INDEX_NONE) {
				node00.m_var_surf = m_vars_surf++;
			}
			if(node01.m_port != INDEX_NONE && node01.m_var_surf == INDEX_NONE) {
				node01.m_var_surf = m_vars_surf++;
			}
			if(node10.m_port != INDEX_NONE && node10.m_var_surf == INDEX_NONE) {
				node10.m_var_surf = m_vars_surf++;
			}
			if(node11.m_port != INDEX_NONE && node11.m_var_surf == INDEX_NONE) {
				node11.m_var_surf = m_vars_surf++;
			}
		}
	}

	// avoid problems later
	if(m_vars_real == 0)
		throw std::runtime_error("GridMesh2D error: The mesh has no free variables.");
	if(m_vars_fixed == 0)
		throw std::runtime_error("GridMesh2D error: The mesh has no fixed variables.");
	if(m_vars_surf == 0)
		throw std::runtime_error("GridMesh2D error: The mesh has no surface variables.");

}

void GridMesh2D::InitSurfaceMatrix() {

	// generate matrix
	m_matrix_surf.Reset(m_vars_surf, 0, 2);
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

			// get cell size
			real_t delta_x = m_grid_x[ix + 1] - m_grid_x[ix];
			real_t delta_y = m_grid_y[iy + 1] - m_grid_y[iy];

			// process conductor surfaces
			Edge &edgeh0 = GetEdgeH(ix, iy);
			Edge &edgeh1 = GetEdgeH(ix, iy + 1);
			Edge &edgev0 = GetEdgeV(ix, iy);
			Edge &edgev1 = GetEdgeV(ix + 1, iy);
			if(edgeh0.m_conductor != INDEX_NONE) {
				assert(node00.m_var_surf != INDEX_NONE);
				assert(node01.m_var_surf != INDEX_NONE);
				real_t coef = 1.0 / 6.0 * delta_x;
				m_matrix_surf.Add(node00.m_var_surf, node00.m_var_surf, coef);
				m_matrix_surf.Add(node01.m_var_surf, node01.m_var_surf, coef);
				m_matrix_surf.Add(node00.m_var_surf, node01.m_var_surf, coef);
			}
			if(edgeh1.m_conductor != INDEX_NONE) {
				assert(node10.m_var_surf != INDEX_NONE);
				assert(node11.m_var_surf != INDEX_NONE);
				real_t coef = 1.0 / 6.0 * delta_x;
				m_matrix_surf.Add(node10.m_var_surf, node10.m_var_surf, coef);
				m_matrix_surf.Add(node11.m_var_surf, node11.m_var_surf, coef);
				m_matrix_surf.Add(node10.m_var_surf, node11.m_var_surf, coef);
			}
			if(edgev0.m_conductor != INDEX_NONE) {
				assert(node00.m_var_surf != INDEX_NONE);
				assert(node10.m_var_surf != INDEX_NONE);
				real_t coef = 1.0 / 6.0 * delta_y;
				m_matrix_surf.Add(node00.m_var_surf, node00.m_var_surf, coef);
				m_matrix_surf.Add(node10.m_var_surf, node10.m_var_surf, coef);
				m_matrix_surf.Add(node00.m_var_surf, node10.m_var_surf, coef);
			}
			if(edgev1.m_conductor != INDEX_NONE) {
				assert(node01.m_var_surf != INDEX_NONE);
				assert(node11.m_var_surf != INDEX_NONE);
				real_t coef = 1.0 / 6.0 * delta_y;
				m_matrix_surf.Add(node01.m_var_surf, node01.m_var_surf, coef);
				m_matrix_surf.Add(node11.m_var_surf, node11.m_var_surf, coef);
				m_matrix_surf.Add(node01.m_var_surf, node11.m_var_surf, coef);
			}

		}
	}

	// factorize surface matrix
	m_cholmod_sparse_surf.Reset(m_matrix_surf);
	m_cholmod_factor_surf.Factorize(m_cholmod_sparse_surf);

}

void GridMesh2D::GenerateMatrices() {
	assert(m_initialized && !m_solved);

	// load conductor properties
	m_conductor_properties.clear();
	m_conductor_properties.resize(m_conductors.size());
	for(size_t i = 0; i < m_conductors.size(); ++i) {
		GetConductorProperties(m_conductor_properties[i], m_conductors[i].m_material, m_solution_frequency);
	}

	// load dielectric properties
	m_dielectric_properties.clear();
	m_dielectric_properties.resize(m_dielectrics.size());
	for(size_t i = 0; i < m_dielectrics.size(); ++i) {
		GetDielectricProperties(m_dielectric_properties[i], m_dielectrics[i].m_material, m_solution_frequency);
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

			// get cell size
			real_t delta_x = m_grid_x[ix + 1] - m_grid_x[ix];
			real_t delta_y = m_grid_y[iy + 1] - m_grid_y[iy];

			// get dielectric properties
			std::complex<real_t> permittivity_x(VACUUM_PERMITTIVITY, 0.0), permittivity_y(VACUUM_PERMITTIVITY, 0.0);
			if(cell.m_dielectric != INDEX_NONE) {
				permittivity_x = m_dielectric_properties[cell.m_dielectric].m_permittivity_x;
				permittivity_y = m_dielectric_properties[cell.m_dielectric].m_permittivity_y;
			}

			// calculate scale factors
			real_t scale_x_e = 1.0 / 6.0 * permittivity_x.real() * delta_y / delta_x;
			real_t scale_y_e = 1.0 / 6.0 * permittivity_y.real() * delta_x / delta_y;
			real_t scale_x_h = 1.0 / 6.0 / VACUUM_PERMEABILITY * delta_y / delta_x;
			real_t scale_y_h = 1.0 / 6.0 / VACUUM_PERMEABILITY * delta_x / delta_y;

			// calculate E coefficients
			real_t coef_s_e = scale_x_e + scale_y_e;
			real_t coef_x_e = -2.0 * scale_x_e + scale_y_e;
			real_t coef_y_e = -2.0 * scale_y_e + scale_x_e;
			real_t coef_d_e = -scale_x_e - scale_y_e;

			// calculate H coefficients
			real_t coef_s_h = scale_x_h + scale_y_h;
			real_t coef_x_h = -2.0 * scale_x_h + scale_y_h;
			real_t coef_y_h = -2.0 * scale_y_h + scale_x_h;
			real_t coef_d_h = -scale_x_h - scale_y_h;

			// add to E matrix
			m_matrix_e.Add(node00.m_var, node00.m_var, coef_s_e);
			m_matrix_e.Add(node01.m_var, node01.m_var, coef_s_e);
			m_matrix_e.Add(node10.m_var, node10.m_var, coef_s_e);
			m_matrix_e.Add(node11.m_var, node11.m_var, coef_s_e);
			m_matrix_e.Add(node00.m_var, node01.m_var, coef_x_e);
			m_matrix_e.Add(node10.m_var, node11.m_var, coef_x_e);
			m_matrix_e.Add(node00.m_var, node10.m_var, coef_y_e);
			m_matrix_e.Add(node01.m_var, node11.m_var, coef_y_e);
			m_matrix_e.Add(node00.m_var, node11.m_var, coef_d_e);
			m_matrix_e.Add(node01.m_var, node10.m_var, coef_d_e);

			// add to H matrix
			m_matrix_h.Add(node00.m_var, node00.m_var, coef_s_h);
			m_matrix_h.Add(node01.m_var, node01.m_var, coef_s_h);
			m_matrix_h.Add(node10.m_var, node10.m_var, coef_s_h);
			m_matrix_h.Add(node11.m_var, node11.m_var, coef_s_h);
			m_matrix_h.Add(node00.m_var, node01.m_var, coef_x_h);
			m_matrix_h.Add(node10.m_var, node11.m_var, coef_x_h);
			m_matrix_h.Add(node00.m_var, node10.m_var, coef_y_h);
			m_matrix_h.Add(node01.m_var, node11.m_var, coef_y_h);
			m_matrix_h.Add(node00.m_var, node11.m_var, coef_d_h);
			m_matrix_h.Add(node01.m_var, node10.m_var, coef_d_h);

		}
	}

}

void GridMesh2D::SolveModes(std::vector<real_t> &charges, std::vector<real_t> &currents, std::vector<real_t> &dielectric_losses, std::vector<real_t> &resistive_losses) {
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

	// calculate the surface charges
	charges.clear();
	charges.resize(m_vars_fixed * m_mode_count);
	for(size_t i = 0; i < m_mode_count; ++i) {
		real_t *residual = charges.data() + m_vars_fixed * i;
		real_t *solution = m_cholmod_solution_e.GetRealData() + m_cholmod_solution_e.GetStride() * i;
		real_t *fixed_values = m_modes.data() + m_vars_fixed * i;
		m_matrix_e.CalculateResidual(residual, solution, fixed_values);
	}

	// factorize H matrix
	m_cholmod_sparse.Reset(m_matrix_h);
	m_cholmod_factor.Factorize(m_cholmod_sparse);

	// solve H matrix
	std::fill_n(m_cholmod_rhs.GetRealData(), m_cholmod_rhs.GetStride() * m_mode_count, 0.0);
	for(size_t i = 0; i < m_mode_count; ++i) {
		m_matrix_h.SubtractFromRhs(m_cholmod_rhs.GetRealData() + m_cholmod_rhs.GetStride() * i, m_modes.data() + m_vars_fixed * i);
	}
	m_cholmod_factor.Solve(m_cholmod_solution_h, m_cholmod_rhs, m_cholmod_ws1, m_cholmod_ws2);

	// calculate the surface currents
	currents.clear();
	currents.resize(m_vars_fixed * m_mode_count);
	for(size_t i = 0; i < m_mode_count; ++i) {
		real_t *residual = currents.data() + m_vars_fixed * i;
		real_t *solution = m_cholmod_solution_h.GetRealData() + m_cholmod_solution_h.GetStride() * i;
		real_t *fixed_values = m_modes.data() + m_vars_fixed * i;
		m_matrix_h.CalculateResidual(residual, solution, fixed_values);
	}

	// calculate surface residuals and dielectric losses
	real_t solution_omega = 2.0 * M_PI * m_solution_frequency;
	dielectric_losses.clear();
	dielectric_losses.resize(m_mode_count);
	m_cholmod_rhs_surf.ResetReal(m_vars_surf, m_mode_count);
	std::fill_n(m_cholmod_rhs_surf.GetRealData(), m_cholmod_rhs_surf.GetStride() * m_mode_count, 0.0);
	for(size_t i = 0; i < m_mode_count; ++i) {
		real_t *solution_e = m_cholmod_solution_e.GetRealData() + m_cholmod_solution_e.GetStride() * i;
		real_t *solution_h = m_cholmod_solution_h.GetRealData() + m_cholmod_solution_h.GetStride() * i;
		real_t *fixed_values = m_modes.data() + m_vars_fixed * i;
		real_t *surf_resid_h = m_cholmod_rhs_surf.GetRealData() + m_cholmod_rhs_surf.GetStride() * i;
		real_t total_dielectric_loss = 0.0;
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

				// get node values
				real_t node00_value_e = (node00.m_var < INDEX_OFFSET)? solution_e[node00.m_var] : fixed_values[node00.m_var - INDEX_OFFSET];
				real_t node01_value_e = (node01.m_var < INDEX_OFFSET)? solution_e[node01.m_var] : fixed_values[node01.m_var - INDEX_OFFSET];
				real_t node10_value_e = (node10.m_var < INDEX_OFFSET)? solution_e[node10.m_var] : fixed_values[node10.m_var - INDEX_OFFSET];
				real_t node11_value_e = (node11.m_var < INDEX_OFFSET)? solution_e[node11.m_var] : fixed_values[node11.m_var - INDEX_OFFSET];
				real_t node00_value_h = (node00.m_var < INDEX_OFFSET)? solution_h[node00.m_var] : fixed_values[node00.m_var - INDEX_OFFSET];
				real_t node01_value_h = (node01.m_var < INDEX_OFFSET)? solution_h[node01.m_var] : fixed_values[node01.m_var - INDEX_OFFSET];
				real_t node10_value_h = (node10.m_var < INDEX_OFFSET)? solution_h[node10.m_var] : fixed_values[node10.m_var - INDEX_OFFSET];
				real_t node11_value_h = (node11.m_var < INDEX_OFFSET)? solution_h[node11.m_var] : fixed_values[node11.m_var - INDEX_OFFSET];

				// get cell size
				real_t delta_x = m_grid_x[ix + 1] - m_grid_x[ix];
				real_t delta_y = m_grid_y[iy + 1] - m_grid_y[iy];

				// get dielectric properties
				std::complex<real_t> permittivity_x(VACUUM_PERMITTIVITY, 0.0), permittivity_y(VACUUM_PERMITTIVITY, 0.0);
				if(cell.m_dielectric != INDEX_NONE) {
					permittivity_x = m_dielectric_properties[cell.m_dielectric].m_permittivity_x;
					permittivity_y = m_dielectric_properties[cell.m_dielectric].m_permittivity_y;
				}

				// calculate scale factors
				real_t scale_x_h = 1.0 / 6.0 / VACUUM_PERMEABILITY * delta_y / delta_x;
				real_t scale_y_h = 1.0 / 6.0 / VACUUM_PERMEABILITY * delta_x / delta_y;

				// calculate H coefficients
				real_t coef_s_h = scale_x_h + scale_y_h;
				real_t coef_x_h = -2.0 * scale_x_h + scale_y_h;
				real_t coef_y_h = -2.0 * scale_y_h + scale_x_h;
				real_t coef_d_h = -scale_x_h - scale_y_h;

				// calculate residual
				if(node00.m_var_surf != INDEX_NONE)
					surf_resid_h[node00.m_var_surf] += node00_value_h * coef_s_h * 2.0 + node01_value_h * coef_x_h + node10_value_h * coef_y_h + node11_value_h * coef_d_h;
				if(node01.m_var_surf != INDEX_NONE)
					surf_resid_h[node01.m_var_surf] += node01_value_h * coef_s_h * 2.0 + node00_value_h * coef_x_h + node11_value_h * coef_y_h + node10_value_h * coef_d_h;
				if(node10.m_var_surf != INDEX_NONE)
					surf_resid_h[node10.m_var_surf] += node10_value_h * coef_s_h * 2.0 + node11_value_h * coef_x_h + node00_value_h * coef_y_h + node01_value_h * coef_d_h;
				if(node11.m_var_surf != INDEX_NONE)
					surf_resid_h[node11.m_var_surf] += node11_value_h * coef_s_h * 2.0 + node10_value_h * coef_x_h + node01_value_h * coef_y_h + node00_value_h * coef_d_h;

				// calculate dielectric loss
				real_t loss_x = -1.0 / 3.0 * permittivity_x.imag() * solution_omega * delta_y / delta_x;
				real_t loss_y = -1.0 / 3.0 * permittivity_y.imag() * solution_omega * delta_x / delta_y;
				real_t dx0 = node01_value_e - node00_value_e, dx1 = node11_value_e - node10_value_e;
				real_t dy0 = node10_value_e - node00_value_e, dy1 = node11_value_e - node01_value_e;
				total_dielectric_loss += (sqr(dx0 + dx1) - dx0 * dx1) * loss_x + (sqr(dy0 + dy1) - dy0 * dy1) * loss_y;

			}
		}
		dielectric_losses[i] = total_dielectric_loss;
	}

	// solve surface matrix
	m_cholmod_factor_surf.Solve(m_cholmod_solution_surf, m_cholmod_rhs_surf, m_cholmod_ws1_surf, m_cholmod_ws2_surf);

	// calculate resistive losses
	resistive_losses.clear();
	resistive_losses.resize(m_mode_count);
	for(size_t i = 0; i < m_mode_count; ++i) {
		real_t *solution_surf = m_cholmod_solution_surf.GetRealData() + m_cholmod_solution_surf.GetStride() * i;
		real_t total_resistive_loss = 0.0;
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

				// get cell size
				real_t delta_x = m_grid_x[ix + 1] - m_grid_x[ix];
				real_t delta_y = m_grid_y[iy + 1] - m_grid_y[iy];

				// process conductor surfaces
				Edge &edgeh0 = GetEdgeH(ix, iy);
				Edge &edgeh1 = GetEdgeH(ix, iy + 1);
				Edge &edgev0 = GetEdgeV(ix, iy);
				Edge &edgev1 = GetEdgeV(ix + 1, iy);
				if(edgeh0.m_conductor != INDEX_NONE) {
					real_t loss_x = 1.0 / 3.0 * delta_x * m_conductor_properties[edgeh0.m_conductor].m_surface_resistivity;
					real_t curr0 = solution_surf[node00.m_var_surf], curr1 = solution_surf[node01.m_var_surf];
					total_resistive_loss += (sqr(curr0 + curr1) - curr0 * curr1) * loss_x;
				}
				if(edgeh1.m_conductor != INDEX_NONE) {
					real_t loss_x = 1.0 / 3.0 * delta_x * m_conductor_properties[edgeh1.m_conductor].m_surface_resistivity;
					real_t curr0 = solution_surf[node10.m_var_surf], curr1 = solution_surf[node11.m_var_surf];
					total_resistive_loss += (sqr(curr0 + curr1) - curr0 * curr1) * loss_x;
				}
				if(edgev0.m_conductor != INDEX_NONE) {
					real_t loss_y = 1.0 / 3.0 * delta_y * m_conductor_properties[edgev0.m_conductor].m_surface_resistivity;
					real_t curr0 = solution_surf[node00.m_var_surf], curr1 = solution_surf[node10.m_var_surf];
					total_resistive_loss += (sqr(curr0 + curr1) - curr0 * curr1) * loss_y;
				}
				if(edgev1.m_conductor != INDEX_NONE) {
					real_t loss_y = 1.0 / 3.0 * delta_y * m_conductor_properties[edgev1.m_conductor].m_surface_resistivity;
					real_t curr0 = solution_surf[node01.m_var_surf], curr1 = solution_surf[node11.m_var_surf];
					total_resistive_loss += (sqr(curr0 + curr1) - curr0 * curr1) * loss_y;
				}

			}
		}
		resistive_losses[i] = total_resistive_loss;
	}

}

void GridMesh2D::GetCellValues(std::vector<real_t> &cell_values, size_t mode, MeshImageType type) {
	assert(m_initialized);
	assert(mode < m_mode_count);
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

void GridMesh2D::GetNodeValues(std::vector<real_t> &node_values, size_t mode, MeshImageType type) {
	assert(m_initialized && m_solved);
	assert(mode < m_mode_count);
	assert(type == MESHIMAGETYPE_FIELD_E || type == MESHIMAGETYPE_FIELD_H || type == MESHIMAGETYPE_CURRENT);
	node_values.resize(m_nodes.size());
	if(type == MESHIMAGETYPE_FIELD_E || type == MESHIMAGETYPE_FIELD_H) {
		CholmodDenseMatrix &solution = (type == MESHIMAGETYPE_FIELD_E)? m_cholmod_solution_e : m_cholmod_solution_h;
		real_t *solution_values = solution.GetRealData() + solution.GetStride() * mode;
		real_t *fixed_values = m_modes.data() + m_vars_fixed * mode;
		for(size_t i = 0; i < m_nodes.size(); ++i) {
			size_t var = m_nodes[i].m_var;
			node_values[i] = (var < INDEX_OFFSET)? solution_values[var] : fixed_values[var - INDEX_OFFSET];
		}
	} else {
		real_t *solution_values = m_cholmod_solution_surf.GetRealData() + m_cholmod_solution_surf.GetStride() * mode;
		real_t max_value = 0.0;
		for(size_t i = 0; i < m_vars_surf; ++i) {
			max_value = std::max(max_value, fabs(solution_values[i]));
		}
		real_t scale = 1.0 / max_value;
		for(size_t i = 0; i < m_nodes.size(); ++i) {
			size_t var = m_nodes[i].m_var_surf;
			node_values[i] = (var == INDEX_NONE)? 0.0 : solution_values[var] * scale;
		}
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

	if(x2 - x1 < std::min(step1, step2)) {
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
