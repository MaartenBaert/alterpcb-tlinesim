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

#include "MatrixMarket.h"
#include "MiscMath.h"
#include "StringHelper.h"

#include <cfloat>

#include <algorithm>
#include <chrono>
#include <iostream>

#define SIMULATION_VERBOSE 1
#define SIMULATION_SAVE_MATRIXMARKET 0

template<class EigenSparseMatrix>
void EigenSparseFree(EigenSparseMatrix &matrix) {
	matrix.resize(0, 0);
	matrix.data().squeeze();
}

template<class SparseMatrix, typename F>
inline void BuildMatrix_Symm1(SparseMatrix &matrix, size_t i0, size_t i1, F cs, F cx) {
	matrix.Insert(i0, i0, cs * 0.5);
	matrix.Insert(i1, i1, cs * 0.5);
	matrix.Insert(i0, i1, cx);
}

template<class SparseMatrix, typename F>
inline void BuildMatrix_Symm2(
		SparseMatrix &matrix, size_t i00, size_t i01, size_t i10, size_t i11,
		F cs, F cx, F cy, F cd) {
	matrix.Insert(i00, i00, cs * 0.5);
	matrix.Insert(i01, i01, cs * 0.5);
	matrix.Insert(i10, i10, cs * 0.5);
	matrix.Insert(i11, i11, cs * 0.5);
	matrix.Insert(i00, i01, cx);
	matrix.Insert(i10, i11, cx);
	matrix.Insert(i00, i10, cy);
	matrix.Insert(i01, i11, cy);
	matrix.Insert(i00, i11, cd);
	matrix.Insert(i01, i10, cd);
}

template<class SparseMatrix, typename F>
inline void BuildMatrix_Asymm2_Sub(
		SparseMatrix &matrix, size_t i, size_t j00, size_t j01, size_t j10, size_t j11,
		F cs, F cx, F cy, F cd) {
	if(i != INDEX_NONE) {
		matrix.Insert(i, j00, cs);
		matrix.Insert(i, j01, cx);
		matrix.Insert(i, j10, cy);
		matrix.Insert(i, j11, cd);
	}
}

template<typename F>
inline void BuildMatrix_Asymm2(
		SparseBlockMatrixC<F> &matrix, size_t i00, size_t i01, size_t i10, size_t i11,
		size_t j00, size_t j01, size_t j10, size_t j11, F cs, F cx, F cy, F cd) {
	BuildMatrix_Asymm2_Sub(matrix, i00, j00, j01, j10, j11, cs, cx, cy, cd);
	BuildMatrix_Asymm2_Sub(matrix, i01, j01, j00, j11, j10, cs, cx, cy, cd);
	BuildMatrix_Asymm2_Sub(matrix, i10, j10, j11, j00, j01, cs, cx, cy, cd);
	BuildMatrix_Asymm2_Sub(matrix, i11, j11, j10, j01, j00, cs, cx, cy, cd);
}

GridMesh2D::GridMesh2D(const Box2D &world_box, const Box2D &world_focus, real_t grid_inc, real_t grid_epsilon) {
	if(!FinitePositive(grid_inc))
		throw std::runtime_error("GridMesh2D error: grid_inc must be positive.");
	if(!FinitePositive(grid_epsilon))
		throw std::runtime_error("GridMesh2D error: grid_epsilon must be positive.");
	m_world_box = world_box.Normalized();
	m_world_focus = world_focus.Normalized();
	m_grid_inc = grid_inc;
	m_grid_epsilon = grid_epsilon;
	m_pml_box = m_world_box;
	m_pml_step = Box2D(REAL_MAX, REAL_MAX, REAL_MAX, REAL_MAX);
	m_pml_attenuation = 0.0;
	m_vars_free = 0;
	m_vars_fixed = 0;
	m_vars_surf = 0;
}

GridMesh2D::~GridMesh2D() {
	// nothing
}

void GridMesh2D::SetPML(const Box2D &box, real_t step, real_t attenuation) {
	SetPML(box, Box2D(step, step, step, step), attenuation);
}

void GridMesh2D::SetPML(const Box2D &box, const Box2D &step, real_t attenuation) {
	m_pml_box = box.Clipped(m_world_box).Normalized();
	m_pml_step = step;
	m_pml_attenuation = attenuation;
}

size_t GridMesh2D::AddPort(GridMesh2D::PortType type, bool infinite_area) {
	if(IsInitialized())
		throw std::runtime_error("GridMesh2D error: Can't add port after initialization.");
	m_ports.emplace_back(type, infinite_area);
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
	m_conductors.emplace_back(box.Clipped(m_world_box).Normalized(), step, material, port);
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
	m_dielectrics.emplace_back(box.Clipped(m_world_box).Normalized(), step, material);
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
					real_t v0 = v00 + (v01 - v00) * frac_x[i];
					real_t v1 = v10 + (v11 - v10) * frac_x[i];
					row_value[i] = v0 + (v1 - v0) * frac_y[j];
				}
			}
			break;
		}
		case MESHIMAGETYPE_EPOT:
		case MESHIMAGETYPE_MPOT: {
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
		/*case MESHIMAGETYPE_CURRENT: {
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
		}*/
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
	std::cerr << "GridMesh2D stats:"
			  << " grid=" << m_grid_x.size() << "x" << m_grid_y.size()
			  << " nodes=" << m_nodes.size()
			  << " cells=" << m_cells.size()
			  << " vars=" << m_vars_free << "+" << m_vars_fixed
			  << std::endl;
	std::cerr << "GridMesh2D init time:"
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
	SolveMatrices();
#if SIMULATION_VERBOSE
	auto t3 = std::chrono::high_resolution_clock::now();
#endif

#if SIMULATION_VERBOSE
	std::cerr << "GridMesh2D solve time:"
			  << " build=" << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count() << "us"
			  << " solve=" << std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count() << "us"
			  << std::endl;
#endif

}

void GridMesh2D::DoCleanup() {
	EigenSparseFree(m_matrix_epot[0]);
	EigenSparseFree(m_matrix_epot[1]);
	EigenSparseFree(m_matrix_epot[2]);
	EigenSparseFree(m_matrix_mpot[0]);
	EigenSparseFree(m_matrix_mpot[1]);
	EigenSparseFree(m_matrix_mpot[2]);
	EigenSparseFree(m_matrix_surf_resid[0]);
	EigenSparseFree(m_matrix_surf_resid[1]);
	EigenSparseFree(m_matrix_surf_curr);
	EigenSparseFree(m_matrix_surf_loss);
	// TODO: m_eigen_chol
	// TODO: m_eigen_chol_surf;
	m_eigen_rhs.resize(0, 0);
	m_eigen_rhs_surf.resize(0, 0);
}

size_t GridMesh2D::GetFixedVariableCount() {
	assert(IsInitialized());
	return m_vars_fixed;
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

	// refine grid
	GridRefine(m_grid_x, grid_x, m_grid_inc, m_grid_epsilon);
	GridRefine(m_grid_y, grid_y, m_grid_inc, m_grid_epsilon);

	if(m_grid_x.size() < 2 || m_grid_y.size() < 2)
		throw std::runtime_error("GridMesh2D error: The mesh must have at least 2 grid lines.");

	// calculate midpoints
	GridMidpoints(m_midpoints_x, m_grid_x);
	GridMidpoints(m_midpoints_y, m_grid_y);

}

void GridMesh2D::InitCells() {
	assert(!IsInitialized());

	// generate nodes and cells
	m_nodes.resize(m_grid_x.size() * m_grid_y.size());
	m_edges_h.resize((m_grid_x.size() - 1) * m_grid_y.size());
	m_edges_v.resize(m_grid_x.size() * (m_grid_y.size() - 1));
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

}

void GridMesh2D::InitVariables() {
	assert(!IsInitialized());

	// assign variables to ports
	for(size_t i = 0; i < m_ports.size(); ++i) {
		Port &port = m_ports[i];
		switch(port.m_type) {
			case PORTTYPE_FIXED: {
				port.m_var = INDEX_OFFSET + m_vars_fixed++;
				break;
			}
			case PORTTYPE_FLOATING: {
				port.m_var = m_vars_free++;
				break;
			}
		}
	}

	// assign variables to nodes
	for(size_t iy = 0; iy < m_grid_y.size(); ++iy) {
		for(size_t ix = 0; ix < m_grid_x.size(); ++ix) {
			Node &node = GetNode(ix, iy);
			if(node.m_port == INDEX_NONE) {
				node.m_var = m_vars_free++;
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
	if(m_vars_free == 0)
		throw std::runtime_error("GridMesh2D error: The mesh has no free variables.");
	if(m_vars_fixed == 0)
		throw std::runtime_error("GridMesh2D error: The mesh has no fixed variables.");
	if(m_vars_surf == 0)
		throw std::runtime_error("GridMesh2D error: The mesh has no surface variables.");

}

void GridMesh2D::BuildMatrices() {
	assert(IsInitialized() && !IsSolved());

	//real_t omega = 2.0 * M_PI * GetFrequency();

	// load conductor properties
	m_conductor_properties.clear();
	m_conductor_properties.resize(m_conductors.size());
	for(size_t i = 0; i < m_conductors.size(); ++i) {
		GetConductorProperties(m_conductor_properties[i], m_conductors[i].m_material, GetFrequency());
	}

	// load dielectric properties
	m_dielectric_properties.clear();
	m_dielectric_properties.resize(m_dielectrics.size());
	for(size_t i = 0; i < m_dielectrics.size(); ++i) {
		GetDielectricProperties(m_dielectric_properties[i], m_dielectrics[i].m_material, GetFrequency());
	}

	// prepare PML
	complex_t pml_mult_x1 = m_pml_attenuation / ((m_pml_box.x1 - m_world_box.x1) * complex_t(1.0, -2.0 * GetFrequency() * (m_pml_box.x1 - m_world_box.x1) / SPEED_OF_LIGHT));
	complex_t pml_mult_x2 = m_pml_attenuation / ((m_world_box.x2 - m_pml_box.x2) * complex_t(1.0, -2.0 * GetFrequency() * (m_world_box.x2 - m_pml_box.x2) / SPEED_OF_LIGHT));
	complex_t pml_mult_y1 = m_pml_attenuation / ((m_pml_box.y1 - m_world_box.y1) * complex_t(1.0, -2.0 * GetFrequency() * (m_pml_box.y1 - m_world_box.y1) / SPEED_OF_LIGHT));
	complex_t pml_mult_y2 = m_pml_attenuation / ((m_world_box.y2 - m_pml_box.y2) * complex_t(1.0, -2.0 * GetFrequency() * (m_world_box.y2 - m_pml_box.y2) / SPEED_OF_LIGHT));

	// allocate matrices
	std::vector<real_t> port_dc_conductances(m_ports.size(), 0.0);
	SparseBlockMatrixCSL<complex_t> matrix_epot, matrix_mpot;
	SparseBlockMatrixC<real_t> matrix_surf_resid;
	SparseMatrixCSL<real_t> matrix_surf_curr, matrix_surf_loss;
	matrix_epot.Reset(m_vars_free, m_vars_fixed, m_vars_free, m_vars_fixed, 5);
	matrix_mpot.Reset(m_vars_free, m_vars_fixed, m_vars_free, m_vars_fixed, 5);
	matrix_surf_resid.Reset(m_vars_surf, 0, m_vars_free, m_vars_fixed, 5);
	matrix_surf_curr.Reset(m_vars_surf, m_vars_surf, 2);
	matrix_surf_loss.Reset(m_vars_surf, m_vars_surf, 2);

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
				real_t conductivity = m_conductor_properties[cell.m_conductor].m_conductivity;
				port_dc_conductances[port] += conductivity * delta_x * delta_y;
				continue;
			}

			// get neighboring nodes
			Node &node00 = GetNode(ix    , iy    );
			Node &node01 = GetNode(ix + 1, iy    );
			Node &node10 = GetNode(ix    , iy + 1);
			Node &node11 = GetNode(ix + 1, iy + 1);

			// get dielectric properties
			complex_t permittivity_x(VACUUM_PERMITTIVITY, 0.0), permittivity_y(VACUUM_PERMITTIVITY, 0.0);
			if(cell.m_dielectric != INDEX_NONE) {
				permittivity_x = m_dielectric_properties[cell.m_dielectric].m_permittivity_x;
				permittivity_y = m_dielectric_properties[cell.m_dielectric].m_permittivity_y;
			}

			// get PML properties
			real_t center_x = 0.5 * (m_grid_x[ix] + m_grid_x[ix + 1]);
			real_t center_y = 0.5 * (m_grid_y[iy] + m_grid_y[iy + 1]);
			complex_t pml_sx = 1.0, pml_sy = 1.0;
			if(center_x < m_pml_box.x1) pml_sx += (m_pml_box.x1 - center_x) * pml_mult_x1;
			if(center_x > m_pml_box.x2) pml_sx += (center_x - m_pml_box.x2) * pml_mult_x2;
			if(center_y < m_pml_box.y1) pml_sy += (m_pml_box.y1 - center_y) * pml_mult_y1;
			if(center_y > m_pml_box.y2) pml_sy += (center_y - m_pml_box.y2) * pml_mult_y2;
			complex_t delta_pml_x = delta_x * pml_sx;
			complex_t delta_pml_y = delta_y * pml_sy;

			// calculate scale factors
			complex_t scale_x_epot = 1.0 / 6.0 * delta_pml_y / delta_pml_x * permittivity_x;
			complex_t scale_y_epot = 1.0 / 6.0 * delta_pml_x / delta_pml_y * permittivity_y;
			complex_t scale_x_mpot = 1.0 / 6.0 * delta_pml_y / delta_pml_x / VACUUM_PERMEABILITY;
			complex_t scale_y_mpot = 1.0 / 6.0 * delta_pml_x / delta_pml_y / VACUUM_PERMEABILITY;

			// calculate electric potential coefficients
			complex_t coef_s_epot = 2.0 * (scale_x_epot + scale_y_epot);
			complex_t coef_x_epot = scale_y_epot - 2.0 * scale_x_epot;
			complex_t coef_y_epot = scale_x_epot - 2.0 * scale_y_epot;
			complex_t coef_d_epot = -(scale_x_epot + scale_y_epot);

			// calculate magnetic potential coefficients
			complex_t coef_s_mpot = 2.0 * (scale_x_mpot + scale_y_mpot);
			complex_t coef_x_mpot = scale_y_mpot - 2.0 * scale_x_mpot;
			complex_t coef_y_mpot = scale_x_mpot - 2.0 * scale_y_mpot;
			complex_t coef_d_mpot = -(scale_x_mpot + scale_y_mpot);

			// add to potential matrices
			BuildMatrix_Symm2(matrix_epot, node00.m_var, node01.m_var, node10.m_var, node11.m_var,
							  coef_s_epot, coef_x_epot, coef_y_epot, coef_d_epot);
			BuildMatrix_Symm2(matrix_mpot, node00.m_var, node01.m_var, node10.m_var, node11.m_var,
							  coef_s_mpot, coef_x_mpot, coef_y_mpot, coef_d_mpot);

			// add to surface residual matrix
			BuildMatrix_Asymm2(matrix_surf_resid,
							   node00.m_var_surf, node01.m_var_surf, node10.m_var_surf, node11.m_var_surf,
							   node00.m_var, node01.m_var, node10.m_var, node11.m_var,
							   coef_s_mpot.real(), coef_x_mpot.real(), coef_y_mpot.real(), coef_d_mpot.real());

			// process conductor surfaces
			Edge &edgeh0 = GetEdgeH(ix, iy);
			Edge &edgeh1 = GetEdgeH(ix, iy + 1);
			Edge &edgev0 = GetEdgeV(ix, iy);
			Edge &edgev1 = GetEdgeV(ix + 1, iy);
			if(edgeh0.m_conductor != INDEX_NONE) {
				assert(node00.m_var_surf != INDEX_NONE);
				assert(node01.m_var_surf != INDEX_NONE);
				real_t coef_curr = 1.0 / 6.0 * delta_x;
				real_t coef_loss = coef_curr / m_conductor_properties[edgeh0.m_conductor].m_surface_conductivity;
				BuildMatrix_Symm1(matrix_surf_curr, node00.m_var_surf, node01.m_var_surf, 2.0 * coef_curr, coef_curr);
				BuildMatrix_Symm1(matrix_surf_loss, node00.m_var_surf, node01.m_var_surf, 2.0 * coef_loss, coef_loss);
			}
			if(edgeh1.m_conductor != INDEX_NONE) {
				assert(node10.m_var_surf != INDEX_NONE);
				assert(node11.m_var_surf != INDEX_NONE);
				real_t coef_curr = 1.0 / 6.0 * delta_x;
				real_t coef_loss = coef_curr / m_conductor_properties[edgeh1.m_conductor].m_surface_conductivity;
				BuildMatrix_Symm1(matrix_surf_curr, node10.m_var_surf, node11.m_var_surf, 2.0 * coef_curr, coef_curr);
				BuildMatrix_Symm1(matrix_surf_loss, node10.m_var_surf, node11.m_var_surf, 2.0 * coef_loss, coef_loss);
			}
			if(edgev0.m_conductor != INDEX_NONE) {
				assert(node00.m_var_surf != INDEX_NONE);
				assert(node10.m_var_surf != INDEX_NONE);
				real_t coef_curr = 1.0 / 6.0 * delta_y;
				real_t coef_loss = coef_curr / m_conductor_properties[edgev0.m_conductor].m_surface_conductivity;
				BuildMatrix_Symm1(matrix_surf_curr, node00.m_var_surf, node10.m_var_surf, 2.0 * coef_curr, coef_curr);
				BuildMatrix_Symm1(matrix_surf_loss, node00.m_var_surf, node10.m_var_surf, 2.0 * coef_loss, coef_loss);
			}
			if(edgev1.m_conductor != INDEX_NONE) {
				assert(node01.m_var_surf != INDEX_NONE);
				assert(node11.m_var_surf != INDEX_NONE);
				real_t coef_curr = 1.0 / 6.0 * delta_y;
				real_t coef_loss = coef_curr / m_conductor_properties[edgev1.m_conductor].m_surface_conductivity;
				BuildMatrix_Symm1(matrix_surf_curr, node01.m_var_surf, node11.m_var_surf, 2.0 * coef_curr, coef_curr);
				BuildMatrix_Symm1(matrix_surf_loss, node01.m_var_surf, node11.m_var_surf, 2.0 * coef_loss, coef_loss);
			}

		}
	}

	// convert port DC conductance to resistance
	m_vector_dc_resistances.resize((Eigen::Index) m_vars_fixed);
	for(size_t i = 0; i < m_ports.size(); ++i) {
		size_t var = m_ports[i].m_var;
		if(var >= INDEX_OFFSET) {
			m_vector_dc_resistances[(Eigen::Index) (var - INDEX_OFFSET)] = (m_ports[i].m_infinite_area)? 0.0 : 1.0 / port_dc_conductances[i];
		}
	}

	// convert to Eigen sparse matrices
	matrix_epot.GetMatrixA().ToEigen(m_matrix_epot[0]);
	matrix_epot.GetMatrixBC().ToEigen(m_matrix_epot[1]);
	matrix_epot.GetMatrixD().ToEigen(m_matrix_epot[2]);
	matrix_mpot.GetMatrixA().ToEigen(m_matrix_mpot[0]);
	matrix_mpot.GetMatrixBC().ToEigen(m_matrix_mpot[1]);
	matrix_mpot.GetMatrixD().ToEigen(m_matrix_mpot[2]);
	matrix_surf_resid.GetMatrixA().ToEigen(m_matrix_surf_resid[0]);
	matrix_surf_resid.GetMatrixB().ToEigen(m_matrix_surf_resid[1]);
	matrix_surf_curr.ToEigen(m_matrix_surf_curr);
	matrix_surf_loss.ToEigen(m_matrix_surf_loss);

#if SIMULATION_SAVE_MATRIXMARKET
	MatrixMarket::Save("matrix_epot.mtx", m_matrix_epot[0], true);
	MatrixMarket::Save("matrix_mpot.mtx", m_matrix_mpot[0], true);
#endif

}

void GridMesh2D::SolveMatrices() {
	assert(IsInitialized() && !IsSolved());

	real_t omega = 2.0 * M_PI * GetFrequency();

	// factorize electric potential matrix
	if(m_eigen_chol.permutationP().size() == 0) {
		m_eigen_chol.analyzePattern(m_matrix_epot[0].real());
	}
	m_eigen_chol.factorize(m_matrix_epot[0].real());
	if(m_eigen_chol.info() != Eigen::Success)
		throw std::runtime_error("Sparse matrix factorization failed!");

	// solve electric potential matrix
	m_eigen_rhs = -m_matrix_epot[1].real().transpose() * GetModes();
	m_eigen_solution_epot = m_eigen_chol.solve(m_eigen_rhs);

	// factorize magnetic potential matrix
	if(m_eigen_chol.permutationP().size() == 0) {
		m_eigen_chol.analyzePattern(m_matrix_mpot[0].real());
	}
	m_eigen_chol.factorize(m_matrix_mpot[0].real());
	if(m_eigen_chol.info() != Eigen::Success)
		throw std::runtime_error("Sparse matrix factorization failed!");

	// solve magnetic potential matrix
	m_eigen_rhs = -m_matrix_mpot[1].real().transpose() * GetModes();
	m_eigen_solution_mpot = m_eigen_chol.solve(m_eigen_rhs);

	// calculate residuals
	Eigen::MatrixXr residual_epot = m_matrix_epot[1].real() * m_eigen_solution_epot + m_matrix_epot[2].real().selfadjointView<Eigen::Lower>() * GetModes();
	Eigen::MatrixXr residual_mpot = m_matrix_mpot[1].real() * m_eigen_solution_mpot + m_matrix_mpot[2].real().selfadjointView<Eigen::Lower>() * GetModes();
	Eigen::MatrixXr charge_matrix = GetModes().transpose() * residual_epot;
	Eigen::MatrixXr current_matrix = GetModes().transpose() * residual_mpot;

	// calculate losses
	Eigen::MatrixXr electric_loss_matrix =
			m_eigen_solution_epot.transpose() * (m_matrix_epot[0].imag().selfadjointView<Eigen::Lower>() * m_eigen_solution_epot + m_matrix_epot[1].imag().transpose() * GetModes()) +
			GetModes().transpose() * (m_matrix_epot[1].imag() * m_eigen_solution_epot + m_matrix_epot[2].imag().selfadjointView<Eigen::Lower>() * GetModes());
	Eigen::MatrixXr magnetic_loss_matrix =
			m_eigen_solution_mpot.transpose() * (m_matrix_mpot[0].imag().selfadjointView<Eigen::Lower>() * m_eigen_solution_mpot + m_matrix_mpot[1].imag().transpose() * GetModes()) +
			GetModes().transpose() * (m_matrix_mpot[1].imag() * m_eigen_solution_mpot + m_matrix_mpot[2].imag().selfadjointView<Eigen::Lower>() * GetModes());

	// factorize current matrix
	if(m_eigen_chol_surf.permutationP().size() == 0) {
		m_eigen_chol_surf.analyzePattern(m_matrix_surf_curr);
	}
	m_eigen_chol_surf.factorize(m_matrix_surf_curr);
	if(m_eigen_chol_surf.info() != Eigen::Success)
		throw std::runtime_error("Sparse matrix factorization failed!");

	// calculate surface currents
	m_eigen_rhs_surf = m_matrix_surf_resid[0] * m_eigen_solution_mpot + m_matrix_surf_resid[1] * GetModes();
	m_eigen_solution_surf = m_eigen_chol_surf.solve(m_eigen_rhs_surf);

	// calculate surface losses
	Eigen::MatrixXr surface_loss_matrix = m_eigen_solution_surf.transpose() * m_matrix_surf_loss.selfadjointView<Eigen::Lower>() * m_eigen_solution_surf;

	// calculate DC losses and combine with surface losses
	Eigen::MatrixXr dc_loss_matrix = residual_mpot.transpose() * m_vector_dc_resistances.asDiagonal() * residual_mpot;
	Eigen::MatrixXr combined_loss_matrix = surface_loss_matrix.cwiseMax(dc_loss_matrix);

	std::cerr << "charges =\n" << charge_matrix << std::endl;
	std::cerr << "currents =\n" << current_matrix << std::endl;
	std::cerr << "electric_loss_matrix =\n" << electric_loss_matrix << std::endl;
	std::cerr << "magnetic_loss_matrix =\n" << magnetic_loss_matrix << std::endl;
	std::cerr << "surface_loss_matrix =\n" << surface_loss_matrix << std::endl;
	std::cerr << "dc_loss_matrix =\n" << dc_loss_matrix << std::endl;
	std::cerr << "combined_loss_matrix =\n" << combined_loss_matrix << std::endl;
	std::cerr << std::endl;

	// calculate L, C, R and G
	m_inductance_matrix = current_matrix.inverse();
	m_capacitance_matrix = charge_matrix;
	m_resistance_matrix = m_inductance_matrix.transpose() * (-omega * magnetic_loss_matrix + combined_loss_matrix) * m_inductance_matrix;
	m_conductance_matrix = -omega * electric_loss_matrix;

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

void GridMesh2D::GetNodeValues(std::vector<real_t> &node_values, size_t mode, MeshImageType type) {
	assert(IsInitialized() && IsSolved());
	assert(mode < GetModeCount());
	assert(type == MESHIMAGETYPE_EPOT || type == MESHIMAGETYPE_MPOT);
	node_values.clear();
	node_values.resize(m_nodes.size(), 0.0);
	Eigen::MatrixXr &solution = (type == MESHIMAGETYPE_EPOT)? m_eigen_solution_epot : m_eigen_solution_mpot;
	real_t *solution_values = solution.data() + solution.outerStride() * (ptrdiff_t) mode;
	const real_t *fixed_values = GetModes().data() + GetModes().outerStride() * (ptrdiff_t) mode;
	real_t max_value = 0.0;
	for(size_t i = 0; i < m_vars_fixed; ++i) {
		max_value = std::max(max_value, fabs(fixed_values[i]));
	}
	real_t scale = 1.0 / max_value;
	for(size_t i = 0; i < m_nodes.size(); ++i) {
		size_t var = m_nodes[i].m_var;
		node_values[i] = (var < INDEX_OFFSET)? solution_values[var] * scale : fixed_values[var - INDEX_OFFSET] * scale;
	}
}

void GridMesh2D::GetCellNodeValues(std::vector<std::array<real_t, 4>> &cellnode_values, size_t mode, MeshImageType type) {
	assert(IsInitialized() && IsSolved());
	assert(mode < GetModeCount());
	assert(type == MESHIMAGETYPE_ENERGY || type == MESHIMAGETYPE_CURRENT);
	cellnode_values.clear();
	cellnode_values.resize(m_cells.size());
	if(type == MESHIMAGETYPE_ENERGY) {
		real_t *solution_epot = m_eigen_solution_epot.data() + m_eigen_solution_epot.outerStride() * (ptrdiff_t) mode;
		real_t *solution_mpot = m_eigen_solution_mpot.data() + m_eigen_solution_mpot.outerStride() * (ptrdiff_t) mode;
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
				real_t node00_value_epot = (node00.m_var < INDEX_OFFSET)? solution_epot[node00.m_var] : fixed_values[node00.m_var - INDEX_OFFSET];
				real_t node01_value_epot = (node01.m_var < INDEX_OFFSET)? solution_epot[node01.m_var] : fixed_values[node01.m_var - INDEX_OFFSET];
				real_t node10_value_epot = (node10.m_var < INDEX_OFFSET)? solution_epot[node10.m_var] : fixed_values[node10.m_var - INDEX_OFFSET];
				real_t node11_value_epot = (node11.m_var < INDEX_OFFSET)? solution_epot[node11.m_var] : fixed_values[node11.m_var - INDEX_OFFSET];
				real_t node00_value_mpot = (node00.m_var < INDEX_OFFSET)? solution_mpot[node00.m_var] : fixed_values[node00.m_var - INDEX_OFFSET];
				real_t node01_value_mpot = (node01.m_var < INDEX_OFFSET)? solution_mpot[node01.m_var] : fixed_values[node01.m_var - INDEX_OFFSET];
				real_t node10_value_mpot = (node10.m_var < INDEX_OFFSET)? solution_mpot[node10.m_var] : fixed_values[node10.m_var - INDEX_OFFSET];
				real_t node11_value_mpot = (node11.m_var < INDEX_OFFSET)? solution_mpot[node11.m_var] : fixed_values[node11.m_var - INDEX_OFFSET];
				real_t ex0 = node01_value_epot - node00_value_epot, ex1 = node11_value_epot - node10_value_epot;
				real_t ey0 = node10_value_epot - node00_value_epot, ey1 = node11_value_epot - node01_value_epot;
				real_t hx0 = node01_value_mpot - node00_value_mpot, hx1 = node11_value_mpot - node10_value_mpot;
				real_t hy0 = node10_value_mpot - node00_value_mpot, hy1 = node11_value_mpot - node01_value_mpot;
				real_t scale_x = 1.0 / square(m_grid_x[ix + 1] - m_grid_x[ix]);
				real_t scale_y = 1.0 / square(m_grid_y[iy + 1] - m_grid_y[iy]);
				size_t dielectric_index = (cell.m_dielectric == INDEX_NONE)? 0 : cell.m_dielectric + 1;
				dielectric_values[dielectric_index][node00_index] += ex0 * hx0 * scale_x + ey0 * hy0 * scale_y;
				dielectric_values[dielectric_index][node01_index] += ex0 * hx0 * scale_x + ey1 * hy1 * scale_y;
				dielectric_values[dielectric_index][node10_index] += ex1 * hx1 * scale_x + ey0 * hy0 * scale_y;
				dielectric_values[dielectric_index][node11_index] += ex1 * hx1 * scale_x + ey1 * hy1 * scale_y;
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
	} else {
		real_t *solution_values = m_eigen_solution_surf.data() + m_eigen_solution_surf.outerStride() * (ptrdiff_t) mode;
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
		}
	}
}

void GridMesh2D::GridAddBox(std::vector<GridLine> &grid_x, std::vector<GridLine> &grid_y, const Box2D &box, const Box2D &step) {
	grid_x.emplace_back(box.x1, step.x1);
	grid_x.emplace_back(box.x2, step.x2);
	grid_y.emplace_back(box.y1, step.y1);
	grid_y.emplace_back(box.y2, step.y2);
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
		size_t n1 = (size_t) std::max<ptrdiff_t>(1, rintp(log2(st1 / (1.0 + (st1 - 1.0) * frac)) / log2frac + 1.5));
		real_t vmin1 = exp2(log2frac * (real_t) n1);
		for(size_t i = n1 - 1; i > 0; --i) {
			result.push_back(x1 + delta1 * (exp2(log2frac * (real_t) i) - vmin1) / (1.0 - vmin1));
		}
		result.push_back(x1 + delta1);
	}

	// second half
	if(delta2 != 0.0) {
		real_t st2 = step2 / delta2;
		size_t n2 = (size_t) std::max<ptrdiff_t>(1, rintp(log2(st2 / (1.0 + (st2 - 1.0) * frac)) / log2frac + 1.5));
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
		frac[i] = clamp((value - grid[j + 1]) * scale + 0.5, 0.0, 1.0);
	}
}
