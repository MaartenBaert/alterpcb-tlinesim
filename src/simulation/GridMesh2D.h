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
#include "Eigen.h"
#include "EigenSparse.h"
#include "GenericMesh.h"
#include "MaterialDatabase.h"
#include "SparseMatrix.h"
#include "Vector.h"

class GridMesh2D : public GenericMesh {

public:
	static constexpr real_t DEFAULT_LAMBDA_FACTOR = 0.2;
	static constexpr real_t DEFAULT_GRID_INC = 0.30;
	static constexpr real_t DEFAULT_GRID_STEP = 0.01;
	//static constexpr real_t DEFAULT_GRID_INC = 0.02;
	//static constexpr real_t DEFAULT_GRID_STEP = 0.05;

private:
	struct Port {
		Vector2D m_anchor;
		bool m_infinite_area;
		size_t m_var_static_e, m_var_static_m, m_var_full_e;
		inline Port(const Vector2D &anchor, bool infinite_area) : m_anchor(anchor), m_infinite_area(infinite_area),
				m_var_static_e(INDEX_NONE), m_var_static_m(INDEX_NONE), m_var_full_e(INDEX_NONE) {}
	};
	struct Conductor {
		Box2D m_box, m_step;
		const MaterialConductor *m_material;
		size_t m_port;
		inline Conductor(const Box2D &box, const Box2D &step, const MaterialConductor *material, size_t port)
			: m_box(box), m_step(step), m_material(material), m_port(port) {}
	};
	struct Dielectric {
		Box2D m_box, m_step;
		const MaterialDielectric *m_material;
		inline Dielectric(const Box2D &box, const Box2D &step, const MaterialDielectric *material)
			: m_box(box), m_step(step), m_material(material) {}
	};
	struct GridLine {
		real_t m_value, m_step;
		inline GridLine(real_t value, real_t step) : m_value(value), m_step(step) {}
		inline bool operator<(const GridLine &other) const { return m_value < other.m_value; }
	};
	struct Node {
		size_t m_conductor;
		bool m_surface;
		size_t m_var_static_e1, m_var_static_m1;
		size_t m_var_full_e1, m_var_full_m1;
		inline Node() : m_conductor(INDEX_NONE), m_surface(false), m_var_static_e1(INDEX_NONE), m_var_static_m1(INDEX_NONE),
				m_var_full_e1(INDEX_NONE), m_var_full_m1(INDEX_NONE) {}
	};
	struct Edge {
		size_t m_conductor;
		bool m_surface;
		size_t m_var_static_e2, m_var_static_m2;
		size_t m_var_full_e2, m_var_full_m2;
		size_t m_var_full_mt1, m_var_full_mt2;
		inline Edge() : m_conductor(INDEX_NONE), m_surface(false), m_var_static_e2(INDEX_NONE), m_var_static_m2(INDEX_NONE),
				m_var_full_e2(INDEX_NONE), m_var_full_m2(INDEX_NONE), m_var_full_mt1(INDEX_NONE), m_var_full_mt2(INDEX_NONE) {}
	};
	struct Cell {
		size_t m_conductor;
		size_t m_dielectric;
		inline Cell() : m_conductor(INDEX_NONE), m_dielectric(INDEX_NONE) {}
	};

public:
	struct NodeField {
		complex_t e1, m1;
	};
	struct EdgeField {
		complex_t e2, m2;
		complex_t mt1, mt2;
	};
	struct SolutionField {
		complex_t m_propagation_constant;
		complex_t m_effective_index;
		std::vector<NodeField> m_field_nodes;
		std::vector<EdgeField> m_field_edges_x, m_field_edges_y;
	};

private:
	SolverType m_solver_type;
	ElementType m_element_type;
	Box2D m_world_box, m_world_focus;
	real_t m_max_frequency, m_lambda_factor;
	real_t m_grid_inc, m_grid_epsilon;

	Box2D m_pml_box, m_pml_step;
	real_t m_pml_attenuation;

	std::vector<Port> m_ports;
	std::vector<Conductor> m_conductors;
	std::vector<Dielectric> m_dielectrics;
	std::vector<Box2D> m_integration_lines;

	std::vector<real_t> m_grid_x, m_grid_y, m_midpoints_x, m_midpoints_y;
	std::vector<Node> m_nodes;
	std::vector<Edge> m_edges_x, m_edges_y;
	std::vector<Cell> m_cells;
	size_t m_vars_static_e_free, m_vars_static_e_fixed;
	size_t m_vars_static_m_free, m_vars_static_m_fixed;
	size_t m_vars_full_em;

	std::vector<MaterialConductor::Properties> m_conductor_properties;
	std::vector<MaterialDielectric::Properties> m_dielectric_properties;
	Eigen::VectorXr m_vector_dc_resistances;
	Eigen::SparseMatrix<complex_t> m_matrix_static_e[4], m_matrix_static_m[4];
	Eigen::SparseMatrix<complex_t> m_matrix_full_em[2];

	Eigen::SparseLU<Eigen::SparseMatrix<complex_t>> m_lu_static_e, m_lu_static_m;
	Eigen::MatrixXc m_solution_static_e, m_solution_static_m, m_solution_full_em;

	Eigen::SparseLU<Eigen::SparseMatrix<complex_t>> m_lu_full_em;

	std::vector<SolutionField> m_solution_fields;

public:
	GridMesh2D(SolverType solver_type, ElementType element_type, const Box2D &world_box, const Box2D &world_focus, real_t max_frequency, real_t lambda_factor, real_t grid_inc, real_t grid_epsilon);
	virtual ~GridMesh2D() override;

	// noncopyable
	GridMesh2D(const GridMesh2D&) = delete;
	GridMesh2D& operator=(const GridMesh2D&) = delete;

	void SetPML(const Box2D &box, real_t step, real_t attenuation);
	void SetPML(const Box2D &box, const Box2D &step, real_t attenuation);

	size_t AddPort(const Vector2D &anchor, bool infinite_area);
	void AddConductor(const Box2D &box, real_t step, const MaterialConductor *material, size_t port);
	void AddConductor(const Box2D &box, const Box2D &step, const MaterialConductor *material, size_t port);
	void AddDielectric(const Box2D &box, real_t step, const MaterialDielectric *material);
	void AddDielectric(const Box2D &box, const Box2D &step, const MaterialDielectric *material);
	void AddIntegrationLine(const Box2D &box);

	virtual Box2D GetWorldBox2D() override;
	virtual Box2D GetWorldFocus2D() override;
	virtual void GetImage2D(std::vector<real_t> &image_value, std::vector<Vector2D> &image_gradient, size_t width, size_t height, const Box2D &view, MeshImageType type, size_t mode) override;

protected:
	virtual void DoInitialize() override;
	virtual void DoSolve() override;
	virtual void DoCleanup() override;
	virtual size_t GetPortCount() override;

private:
	void InitGrid();
	void InitCells();
	void InitVariables();
	void BuildMatrices();
	void SolveStaticModes();
	void SolveStaticEigenModes();
	void SolveFullEigenModes();

	void GetCellField(size_t mode, size_t ix, size_t iy, real_t fx, real_t fy, complex_t &fex, complex_t &fey, complex_t &fez, complex_t &fmx, complex_t &fmy, complex_t &fmz);

	// TODO: remove
	void GetCellValues(std::vector<real_t> &cell_values, size_t mode, MeshImageType type);
	void GetNodeValues(std::vector<real_t> &node_values, size_t mode, MeshImageType type);
	void GetCellNodeValues(std::vector<std::array<real_t, 4>> &cellnode_values, size_t mode, MeshImageType type);

private:
	inline size_t GetNodeIndex(size_t ix, size_t iy) { assert(ix < m_grid_x.size()); assert(iy < m_grid_y.size()); return ix + iy * m_grid_x.size(); }
	inline size_t GetEdgeXIndex(size_t ix, size_t iy) { assert(ix < m_grid_x.size() - 1); assert(iy < m_grid_y.size()); return ix + iy * (m_grid_x.size() - 1); }
	inline size_t GetEdgeYIndex(size_t ix, size_t iy) { assert(ix < m_grid_x.size()); assert(iy < m_grid_y.size() - 1); return ix + iy * m_grid_x.size(); }
	inline size_t GetCellIndex(size_t ix, size_t iy) { assert(ix < m_grid_x.size() - 1); assert(iy < m_grid_y.size() - 1); return ix + iy * (m_grid_x.size() - 1); }
	inline Node& GetNode(size_t ix, size_t iy) { return m_nodes[GetNodeIndex(ix, iy)]; }
	inline Edge& GetEdgeX(size_t ix, size_t iy) { return m_edges_x[GetEdgeXIndex(ix, iy)]; }
	inline Edge& GetEdgeY(size_t ix, size_t iy) { return m_edges_y[GetEdgeYIndex(ix, iy)]; }
	inline Cell& GetCell(size_t ix, size_t iy) { return m_cells[GetCellIndex(ix, iy)]; }

private:
	static void GridAddBox(std::vector<GridLine> &grid_x, std::vector<GridLine> &grid_y, const Box2D &box, const Box2D &step);
	static void GridRefine(std::vector<real_t> &result, std::vector<GridLine> &grid, const std::vector<real_t> &max_step, real_t inc, real_t epsilon);
	static void GridRefine2(std::vector<real_t> &result, real_t x1, real_t x2, real_t step1, real_t step2, real_t max_step, real_t inc);
	static void GridMidpoints(std::vector<real_t> &result, const std::vector<real_t> &grid);
	static void GridMidpoints(std::vector<real_t> &result, const std::vector<GridLine> &grid);
	static void PrepareNodeImage(std::vector<size_t> &index, std::vector<real_t> &frac, const std::vector<real_t> &grid, real_t value1, real_t value2, size_t size);
	static void PrepareCellImage(std::vector<size_t> &index, std::vector<real_t> &frac, const std::vector<real_t> &grid, const std::vector<real_t> &midpoints, real_t value1, real_t value2, size_t size);

};
