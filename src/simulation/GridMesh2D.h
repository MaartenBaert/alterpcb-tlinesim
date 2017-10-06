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
#include "CholmodSolver.h"
#include "GenericMesh.h"
#include "MaterialDatabase.h"
#include "SparseMatrix.h"
#include "Vector.h"

class GridMesh2D : public GenericMesh {

public:
	enum PortType {
		PORTTYPE_FIXED,
		PORTTYPE_FLOATING,
	};

private:
	struct Port {
		PortType m_type;
		index_t m_var;
		inline Port(PortType type) : m_type(type), m_var(INDEX_NONE) {}
	};
	struct Conductor {
		Box2D m_box;
		real_t m_step_left, m_step_right, m_step_top, m_step_bottom;
		const MaterialConductor *m_material;
		index_t m_port;
		inline Conductor(const Box2D &box, real_t step_left, real_t step_right, real_t step_top, real_t step_bottom, const MaterialConductor *material, index_t port)
			: m_box(box), m_step_left(step_left), m_step_right(step_right), m_step_top(step_top), m_step_bottom(step_bottom), m_material(material), m_port(port) {}
	};
	struct Dielectric {
		Box2D m_box;
		real_t m_step_left, m_step_right, m_step_top, m_step_bottom;
		const MaterialDielectric *m_material;
		inline Dielectric(const Box2D &box, real_t step_left, real_t step_right, real_t step_top, real_t step_bottom, const MaterialDielectric *material)
			: m_box(box), m_step_left(step_left), m_step_right(step_right), m_step_top(step_top), m_step_bottom(step_bottom), m_material(material) {}
	};
	struct GridLine {
		real_t m_value,  m_step;
		inline GridLine(real_t value, real_t step) : m_value(value), m_step(step) {}
		inline bool operator<(const GridLine &other) const { return m_value < other.m_value; }
	};
	struct Node {
		index_t m_port;
		index_t m_var;
		inline Node() : m_port(INDEX_NONE), m_var(INDEX_NONE) {}
	};
	struct Cell {
		index_t m_conductor;
		index_t m_dielectric;
		inline Cell() : m_conductor(INDEX_NONE), m_dielectric(INDEX_NONE) {}
	};

private:
	Box2D m_world_box, m_world_focus;
	real_t m_grid_inc, m_grid_epsilon;

	std::vector<Port> m_ports;
	std::vector<Conductor> m_conductors;
	std::vector<Dielectric> m_dielectrics;

	bool m_initialized;
	std::vector<real_t> m_grid_x, m_grid_y, m_midpoints_x, m_midpoints_y;
	std::vector<Node> m_nodes;
	std::vector<Cell> m_cells;
	index_t m_vars_real, m_vars_fixed;

	bool m_solved;
	real_t m_solution_frequency;
	index_t m_mode_count;
	std::vector<real_t> m_modes;
	SymmetricSparseMatrix<real_t> m_matrix_e, m_matrix_h;
	CholmodSparseMatrix m_cholmod_sparse;
	CholmodFactorization m_cholmod_factor;
	CholmodDenseMatrix m_cholmod_rhs, m_cholmod_ws1, m_cholmod_ws2;
	CholmodDenseMatrix m_cholmod_solution_e, m_cholmod_solution_h;

public:
	GridMesh2D(const Box2D &world_box, const Box2D &world_focus, real_t grid_inc, real_t grid_epsilon);
	virtual ~GridMesh2D();

	// noncopyable
	GridMesh2D(const GridMesh2D&) = delete;
	GridMesh2D& operator=(const GridMesh2D&) = delete;

	index_t AddPort(PortType type);
	void AddConductor(const Box2D &box, real_t step, const MaterialConductor *material, index_t port);
	void AddConductor(const Box2D &box, real_t step_left, real_t step_right, real_t step_top, real_t step_bottom, const MaterialConductor *material, index_t port);
	void AddDielectric(const Box2D &box, real_t step, const MaterialDielectric *material);
	void AddDielectric(const Box2D &box, real_t step_left, real_t step_right, real_t step_top, real_t step_bottom, const MaterialDielectric *material);

	virtual void Initialize() override;
	virtual void Solve(std::vector<real_t> &charges, std::vector<real_t> &currents, const std::vector<real_t> &modes, index_t mode_count, real_t frequency) override;
	virtual void Cleanup() override;

	virtual Box2D GetWorldBox2D() override;
	virtual Box2D GetWorldFocus2D() override;
	virtual bool IsInitialized() override;
	virtual bool IsSolved() override;
	virtual index_t GetModeCount() override;
	virtual bool GetImage2D(std::vector<real_t> &image_value, std::vector<Vector2D> &image_gradient, size_t width, size_t height, const Box2D &view, MeshImageType type, index_t mode) override;

private:
	void InitGrid();
	void InitCells();
	void InitVariables();
	void GenerateMatrices();
	void SolveModes(std::vector<real_t> &charges, std::vector<real_t> &currents);

	void GetCellValues(std::vector<real_t> &cell_values, index_t mode, MeshImageType type);
	void GetNodeValues(std::vector<real_t> &node_values, index_t mode, MeshImageType type);

private:
	inline index_t GetNodeIndex(index_t ix, index_t iy) { assert(ix < m_grid_x.size()); assert(iy < m_grid_y.size()); return ix + iy * m_grid_x.size(); }
	inline index_t GetCellIndex(index_t ix, index_t iy) { assert(ix < m_grid_x.size() - 1); assert(iy < m_grid_y.size() - 1); return ix + iy * (m_grid_x.size() - 1); }
	inline Node& GetNode(index_t ix, index_t iy) { return m_nodes[GetNodeIndex(ix, iy)]; }
	inline Cell& GetCell(index_t ix, index_t iy) { return m_cells[GetCellIndex(ix, iy)]; }

private:
	static void GridAddBox(std::vector<GridLine> &grid_x, std::vector<GridLine> &grid_y, const Box2D &box, real_t step_left, real_t step_right, real_t step_top, real_t step_bottom);
	static void GridRefine(std::vector<real_t> &m_cholmod_result, std::vector<GridLine> &grid, real_t inc, real_t epsilon);
	static void GridRefine2(std::vector<real_t> &m_cholmod_result, real_t x1, real_t x2, real_t step1, real_t step2, real_t inc);
	static void GridMidpoints(std::vector<real_t> &m_cholmod_result, const std::vector<real_t> &grid);
	static void PrepareNodeImage(std::vector<index_t> &index, std::vector<real_t> &frac, const std::vector<real_t> &grid, real_t value1, real_t value2, index_t size);
	static void PrepareCellImage(std::vector<index_t> &index, std::vector<real_t> &frac, const std::vector<real_t> &grid, const std::vector<real_t> &midpoints, real_t value1, real_t value2, index_t size);

};
