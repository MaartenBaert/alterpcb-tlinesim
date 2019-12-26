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

#include "TLineTypes.h"

#include "GridMesh2D.h"
#include "Json.h"
#include "MaterialDatabase.h"

#include <iostream> // TODO: remove

void TLine_Stripline_Single(TLineContext &context) {

	VDataDictReader root(context.m_parameters);

	real_t track_width = root.GetMember("track_width").AsFloat() * 1e-3;
	real_t track_thickness = root.GetMember("track_thickness").AsFloat() * 1e-3;
	const MaterialConductor *track_material = FindConductor(root, "track_material", context.m_material_database);
	real_t substrate_thickness_1 = root.GetMember("substrate_thickness_1").AsFloat() * 1e-3;
	const MaterialDielectric *substrate_material_1 = FindDielectric(root, "substrate_material_1", context.m_material_database);
	real_t substrate_thickness_2 = root.GetMember("substrate_thickness_2").AsFloat() * 1e-3;
	const MaterialDielectric *substrate_material_2 = FindDielectric(root, "substrate_material_2", context.m_material_database);
	bool reverse_buildup = root.GetMember("reverse_buildup").AsBool();

	real_t substrate_thickness_total = substrate_thickness_1 + substrate_thickness_2;
	real_t space_x = (track_width + track_thickness + substrate_thickness_1 + substrate_thickness_2) * 10.0;
	Box2D track_box = {
		-0.5 * track_width,
		0.5 * track_width,
		substrate_thickness_2 - ((reverse_buildup)? track_thickness : 0.0),
		substrate_thickness_2 + ((reverse_buildup)? 0.0 : track_thickness),
	};
	Box2D world_box = {
		track_box.x1 - space_x,
		track_box.x2 + space_x,
		0.0,
		substrate_thickness_total,
	};
	Box2D world_focus = {
		track_box.x1,
		track_box.x2,
		world_box.y1,
		world_box.y2,
	};
	Box2D ground1_box = {
		world_box.x1,
		world_box.x2,
		world_box.y1,
		world_box.y1,
	};
	Box2D ground2_box = {
		world_box.x1,
		world_box.x2,
		world_box.y2,
		world_box.y2,
	};
	Box2D substrate1_box = {
		world_box.x1,
		world_box.x2,
		substrate_thickness_2,
		substrate_thickness_total,
	};
	Box2D substrate2_box = {
		world_box.x1,
		world_box.x2,
		0.0,
		substrate_thickness_2,
	};

	real_t critical_dimension = vmin(track_width, substrate_thickness_1, substrate_thickness_2);
	real_t step0 = REAL_MAX, step1 = critical_dimension * GridMesh2D::DEFAULT_GRID_STEP / context.m_mesh_detail;

	std::unique_ptr<GridMesh2D> mesh(new GridMesh2D(world_box, world_focus, GridMesh2D::DEFAULT_GRID_INC / context.m_mesh_detail, critical_dimension * 1.0e-6));

	size_t port_ground = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED, true);
	size_t port_signal = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED, false);

	mesh->AddConductor(ground1_box, step0, track_material, port_ground);
	mesh->AddConductor(ground2_box, step0, track_material, port_ground);
	mesh->AddConductor(track_box, step1, track_material, port_signal);
	mesh->AddDielectric(substrate1_box, step0, substrate_material_1);
	mesh->AddDielectric(substrate2_box, step0, substrate_material_2);

	context.m_output_mesh = std::move(mesh);

	Eigen::MatrixXr modes(2, 1);
	modes.col(0) << 0.0, 1.0;
	TLineSolveModes(context, modes);

}

void TLine_Stripline_Single_Shielded(TLineContext &context) {

	VDataDictReader root(context.m_parameters);

	real_t track_width = root.GetMember("track_width").AsFloat() * 1e-3;
	real_t ground_spacing = root.GetMember("ground_spacing").AsFloat() * 1e-3;
	real_t track_thickness = root.GetMember("track_thickness").AsFloat() * 1e-3;
	const MaterialConductor *track_material = FindConductor(root, "track_material", context.m_material_database);
	real_t substrate_thickness_1 = root.GetMember("substrate_thickness_1").AsFloat() * 1e-3;
	const MaterialDielectric *substrate_material_1 = FindDielectric(root, "substrate_material_1", context.m_material_database);
	real_t substrate_thickness_2 = root.GetMember("substrate_thickness_2").AsFloat() * 1e-3;
	const MaterialDielectric *substrate_material_2 = FindDielectric(root, "substrate_material_2", context.m_material_database);
	bool reverse_buildup = root.GetMember("reverse_buildup").AsBool();

	real_t substrate_thickness_total = substrate_thickness_1 + substrate_thickness_2;
	real_t space_x = (track_width + track_thickness + substrate_thickness_1 + substrate_thickness_2) * 10.0;
	Box2D track_box = {
		-0.5 * track_width,
		0.5 * track_width,
		substrate_thickness_2 - ((reverse_buildup)? track_thickness : 0.0),
		substrate_thickness_2 + ((reverse_buildup)? 0.0 : track_thickness),
	};
	Box2D ground3_box = {
		track_box.x1 - ground_spacing - space_x,
		track_box.x1 - ground_spacing,
		substrate_thickness_2 - ((reverse_buildup)? track_thickness : 0.0),
		substrate_thickness_2 + ((reverse_buildup)? 0.0 : track_thickness),
	};
	Box2D ground4_box = ground3_box.MirroredX();
	Box2D world_box = {
		track_box.x1 - space_x,
		track_box.x2 + space_x,
		0.0,
		substrate_thickness_total,
	};
	Box2D world_focus = {
		ground3_box.x2,
		ground4_box.x1,
		world_box.y1,
		world_box.y2,
	};
	Box2D ground1_box = {
		world_box.x1,
		world_box.x2,
		world_box.y1,
		world_box.y1,
	};
	Box2D ground2_box = {
		world_box.x1,
		world_box.x2,
		world_box.y2,
		world_box.y2,
	};
	Box2D substrate1_box = {
		world_box.x1,
		world_box.x2,
		substrate_thickness_2,
		substrate_thickness_total,
	};
	Box2D substrate2_box = {
		world_box.x1,
		world_box.x2,
		0.0,
		substrate_thickness_2,
	};

	real_t critical_dimension = vmin(track_width, ground_spacing, substrate_thickness_1, substrate_thickness_2);
	real_t step0 = REAL_MAX, step1 = critical_dimension * GridMesh2D::DEFAULT_GRID_STEP / context.m_mesh_detail;

	std::unique_ptr<GridMesh2D> mesh(new GridMesh2D(world_box, world_focus, GridMesh2D::DEFAULT_GRID_INC / context.m_mesh_detail, critical_dimension * 1.0e-6));

	size_t port_ground = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED, true);
	size_t port_signal = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED, false);

	mesh->AddConductor(ground1_box, step0, track_material, port_ground);
	mesh->AddConductor(ground2_box, step0, track_material, port_ground);
	mesh->AddConductor(track_box, step1, track_material, port_signal);
	mesh->AddConductor(ground3_box, Box2D(step0, step1, step1, step1), track_material, port_ground);
	mesh->AddConductor(ground4_box, Box2D(step1, step0, step1, step1), track_material, port_ground);
	mesh->AddDielectric(substrate1_box, step0, substrate_material_1);
	mesh->AddDielectric(substrate2_box, step0, substrate_material_2);

	context.m_output_mesh = std::move(mesh);

	Eigen::MatrixXr modes(2, 1);
	modes.col(0) << 0.0, 1.0;
	TLineSolveModes(context, modes);

}

void TLine_Stripline_Differential(TLineContext &context) {

	VDataDictReader root(context.m_parameters);

	real_t track_width = root.GetMember("track_width").AsFloat() * 1e-3;
	real_t track_spacing = root.GetMember("track_spacing").AsFloat() * 1e-3;
	real_t track_thickness = root.GetMember("track_thickness").AsFloat() * 1e-3;
	const MaterialConductor *track_material = FindConductor(root, "track_material", context.m_material_database);
	real_t substrate_thickness_1 = root.GetMember("substrate_thickness_1").AsFloat() * 1e-3;
	const MaterialDielectric *substrate_material_1 = FindDielectric(root, "substrate_material_1", context.m_material_database);
	real_t substrate_thickness_2 = root.GetMember("substrate_thickness_2").AsFloat() * 1e-3;
	const MaterialDielectric *substrate_material_2 = FindDielectric(root, "substrate_material_2", context.m_material_database);
	bool reverse_buildup = root.GetMember("reverse_buildup").AsBool();

	real_t substrate_thickness_total = substrate_thickness_1 + substrate_thickness_2;
	real_t space_x = (track_width * 2 + track_spacing + track_thickness + substrate_thickness_1 + substrate_thickness_2) * 10.0;
	Box2D track1_box = {
		-0.5 * track_spacing - track_width,
		-0.5 * track_spacing,
		substrate_thickness_2 - ((reverse_buildup)? track_thickness : 0.0),
		substrate_thickness_2 + ((reverse_buildup)? 0.0 : track_thickness),
	};
	Box2D track2_box = track1_box.MirroredX();
	Box2D world_box = {
		track1_box.x1 - space_x,
		track2_box.x2 + space_x,
		0.0,
		substrate_thickness_total,
	};
	Box2D world_focus = {
		track1_box.x1,
		track2_box.x2,
		world_box.y1,
		world_box.y2,
	};
	Box2D ground1_box = {
		world_box.x1,
		world_box.x2,
		world_box.y1,
		world_box.y1,
	};
	Box2D ground2_box = {
		world_box.x1,
		world_box.x2,
		world_box.y2,
		world_box.y2,
	};
	Box2D substrate1_box = {
		world_box.x1,
		world_box.x2,
		substrate_thickness_2,
		substrate_thickness_total,
	};
	Box2D substrate2_box = {
		world_box.x1,
		world_box.x2,
		0.0,
		substrate_thickness_2,
	};

	real_t critical_dimension = vmin(track_width, track_spacing, substrate_thickness_1, substrate_thickness_2);
	real_t step0 = REAL_MAX, step1 = critical_dimension * GridMesh2D::DEFAULT_GRID_STEP / context.m_mesh_detail;

	std::unique_ptr<GridMesh2D> mesh(new GridMesh2D(world_box, world_focus, GridMesh2D::DEFAULT_GRID_INC / context.m_mesh_detail, critical_dimension * 1.0e-6));

	size_t port_ground = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED, true);
	size_t port_signal1 = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED, false);
	size_t port_signal2 = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED, false);

	mesh->AddConductor(ground1_box, step0, track_material, port_ground);
	mesh->AddConductor(ground2_box, step0, track_material, port_ground);
	mesh->AddConductor(track1_box, step1, track_material, port_signal1);
	mesh->AddConductor(track2_box, step1, track_material, port_signal2);
	mesh->AddDielectric(substrate1_box, step0, substrate_material_1);
	mesh->AddDielectric(substrate2_box, step0, substrate_material_2);

	context.m_output_mesh = std::move(mesh);

	Eigen::MatrixXr modes(3, 2);
	modes.col(0) << 0.0, 0.5, -0.5;
	modes.col(1) << 0.0, 1.0, 1.0;
	TLineSolveModes(context, modes);

}

void TLine_Stripline_Differential_Shielded(TLineContext &context) {

	VDataDictReader root(context.m_parameters);

	real_t track_width = root.GetMember("track_width").AsFloat() * 1e-3;
	real_t track_spacing = root.GetMember("track_spacing").AsFloat() * 1e-3;
	real_t ground_spacing = root.GetMember("ground_spacing").AsFloat() * 1e-3;
	real_t track_thickness = root.GetMember("track_thickness").AsFloat() * 1e-3;
	const MaterialConductor *track_material = FindConductor(root, "track_material", context.m_material_database);
	real_t substrate_thickness_1 = root.GetMember("substrate_thickness_1").AsFloat() * 1e-3;
	const MaterialDielectric *substrate_material_1 = FindDielectric(root, "substrate_material_1", context.m_material_database);
	real_t substrate_thickness_2 = root.GetMember("substrate_thickness_2").AsFloat() * 1e-3;
	const MaterialDielectric *substrate_material_2 = FindDielectric(root, "substrate_material_2", context.m_material_database);
	bool reverse_buildup = root.GetMember("reverse_buildup").AsBool();

	real_t substrate_thickness_total = substrate_thickness_1 + substrate_thickness_2;
	real_t space_x = (track_width * 2 + track_spacing + track_thickness + substrate_thickness_1 + substrate_thickness_2) * 10.0;
	Box2D track1_box = {
		-0.5 * track_spacing - track_width,
		-0.5 * track_spacing,
		substrate_thickness_2 - ((reverse_buildup)? track_thickness : 0.0),
		substrate_thickness_2 + ((reverse_buildup)? 0.0 : track_thickness),
	};
	Box2D track2_box = track1_box.MirroredX();
	Box2D ground3_box = {
		track1_box.x1 - ground_spacing - space_x,
		track1_box.x1 - ground_spacing,
		substrate_thickness_2 - ((reverse_buildup)? track_thickness : 0.0),
		substrate_thickness_2 + ((reverse_buildup)? 0.0 : track_thickness),
	};
	Box2D ground4_box = ground3_box.MirroredX();
	Box2D world_box = {
		track1_box.x1 - space_x,
		track2_box.x2 + space_x,
		0.0,
		substrate_thickness_total,
	};
	Box2D world_focus = {
		ground3_box.x2,
		ground4_box.x1,
		world_box.y1,
		world_box.y2,
	};
	Box2D ground1_box = {
		world_box.x1,
		world_box.x2,
		world_box.y1,
		world_box.y1,
	};
	Box2D ground2_box = {
		world_box.x1,
		world_box.x2,
		world_box.y2,
		world_box.y2,
	};
	Box2D substrate1_box = {
		world_box.x1,
		world_box.x2,
		substrate_thickness_2,
		substrate_thickness_total,
	};
	Box2D substrate2_box = {
		world_box.x1,
		world_box.x2,
		0.0,
		substrate_thickness_2,
	};

	real_t critical_dimension = vmin(track_width, track_spacing, ground_spacing, substrate_thickness_1, substrate_thickness_2);
	real_t step0 = REAL_MAX, step1 = critical_dimension * GridMesh2D::DEFAULT_GRID_STEP / context.m_mesh_detail;

	std::unique_ptr<GridMesh2D> mesh(new GridMesh2D(world_box, world_focus, GridMesh2D::DEFAULT_GRID_INC / context.m_mesh_detail, critical_dimension * 1.0e-6));

	size_t port_ground = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED, true);
	size_t port_signal1 = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED, false);
	size_t port_signal2 = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED, false);

	mesh->AddConductor(ground1_box, step0, track_material, port_ground);
	mesh->AddConductor(ground2_box, step0, track_material, port_ground);
	mesh->AddConductor(track1_box, step1, track_material, port_signal1);
	mesh->AddConductor(track2_box, step1, track_material, port_signal2);
	mesh->AddConductor(ground3_box, Box2D(step0, step1, step1, step1), track_material, port_ground);
	mesh->AddConductor(ground4_box, Box2D(step1, step0, step1, step1), track_material, port_ground);
	mesh->AddDielectric(substrate1_box, step0, substrate_material_1);
	mesh->AddDielectric(substrate2_box, step0, substrate_material_2);

	context.m_output_mesh = std::move(mesh);

	Eigen::MatrixXr modes(3, 2);
	modes.col(0) << 0.0, 0.5, -0.5;
	modes.col(1) << 0.0, 1.0, 1.0;
	TLineSolveModes(context, modes);

}

void RegisterTLine_Stripline() {

	VData default_track_width = Json::FromString("1.0");
	VData default_track_spacing = Json::FromString("1.0");
	VData default_ground_spacing = Json::FromString("1.0");
	VData default_track_thickness = Json::FromString("0.035");
	VData default_substrate_thickness = Json::FromString("0.5");
	VData default_track_material = "Copper";
	VData default_substrate_material = "Isola DE104";

	g_tline_types.push_back(TLineType{
		"Stripline (single)",
		"A single track embedded in a substrate, with ground planes above and below. "
		"Striplines have exceptionally low crosstalk and radiation loss since they are completely shielded. "
		"Ground vias should be placed on both sides of the track at regular intervals (less than 1/10th of the wavelength) to ensure correct behavior. "
		"Isolation can be improved by adding more vias.",
		{
			{"Track Width"          , TLINE_PARAMETERTYPE_REAL               , default_track_width        , true , 0},
			{"Track Thickness"      , TLINE_PARAMETERTYPE_REAL               , default_track_thickness    , true , 0},
			{"Track Material"       , TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR , default_track_material     , false, 0},
			{"Reverse Buildup"      , TLINE_PARAMETERTYPE_BOOL               , false                      , false, 1},
			{"Substrate Thickness 1", TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness, true , 0},
			{"Substrate Material 1" , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material , false, 1},
			{"Substrate Thickness 2", TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness, true , 0},
			{"Substrate Material 2" , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material , false, 0},
		},
		{"Single-ended"},
		&TLine_Stripline_Single,
	});

	g_tline_types.push_back(TLineType{
		"Stripline (single, shielded)",
		"A single track embedded in a substrate, with ground planes above, below and on both sides. "
		"Striplines have exceptionally low crosstalk and radiation loss since they are completely shielded. "
		"Ground vias should be placed on both sides of the track at regular intervals (less than 1/10th of the wavelength) to ensure correct behavior. "
		"Isolation can be improved by adding more vias.",
		{
			{"Track Width"          , TLINE_PARAMETERTYPE_REAL               , default_track_width        , true , 0},
			{"Ground Spacing"       , TLINE_PARAMETERTYPE_REAL               , default_ground_spacing     , true , 0},
			{"Track Thickness"      , TLINE_PARAMETERTYPE_REAL               , default_track_thickness    , true , 0},
			{"Track Material"       , TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR , default_track_material     , false, 0},
			{"Reverse Buildup"      , TLINE_PARAMETERTYPE_BOOL               , false                      , false, 1},
			{"Substrate Thickness 1", TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness, true , 0},
			{"Substrate Material 1" , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material , false, 1},
			{"Substrate Thickness 2", TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness, true , 0},
			{"Substrate Material 2" , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material , false, 0},
		},
		{"Single-ended"},
		&TLine_Stripline_Single_Shielded,
	});

	g_tline_types.push_back(TLineType{
		"Stripline (differential)",
		"A differential pair embedded in a substrate, with ground planes above and below. "
		"Striplines have exceptionally low crosstalk and radiation loss since they are completely shielded. "
		"Ground vias should be placed on both sides of the tracks at regular intervals (less than 1/10th of the wavelength) to ensure correct behavior. "
		"Isolation can be improved by adding more vias.",
		{
			{"Track Width"          , TLINE_PARAMETERTYPE_REAL               , default_track_width        , true , 0},
			{"Track Spacing"        , TLINE_PARAMETERTYPE_REAL               , default_track_spacing      , true , 0},
			{"Track Thickness"      , TLINE_PARAMETERTYPE_REAL               , default_track_thickness    , true , 0},
			{"Track Material"       , TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR , default_track_material     , false, 0},
			{"Reverse Buildup"      , TLINE_PARAMETERTYPE_BOOL               , false                      , false, 1},
			{"Substrate Thickness 1", TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness, true , 0},
			{"Substrate Material 1" , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material , false, 1},
			{"Substrate Thickness 2", TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness, true , 0},
			{"Substrate Material 2" , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material , false, 0},
		},
		{"Differential", "Common-mode"},
		&TLine_Stripline_Differential,
	});

	g_tline_types.push_back(TLineType{
		"Stripline (differential, shielded)",
		"A differential pair embedded in a substrate, with ground planes above, below and on both sides. "
		"Striplines have exceptionally low crosstalk and radiation loss since they are completely shielded. "
		"Ground vias should be placed on both sides of the tracks at regular intervals (less than 1/10th of the wavelength) to ensure correct behavior. "
		"Isolation can be improved by adding more vias.",
		{
			{"Track Width"          , TLINE_PARAMETERTYPE_REAL               , default_track_width        , true , 0},
			{"Track Spacing"        , TLINE_PARAMETERTYPE_REAL               , default_track_spacing      , true , 0},
			{"Ground Spacing"       , TLINE_PARAMETERTYPE_REAL               , default_ground_spacing     , true , 0},
			{"Track Thickness"      , TLINE_PARAMETERTYPE_REAL               , default_track_thickness    , true , 0},
			{"Track Material"       , TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR , default_track_material     , false, 0},
			{"Reverse Buildup"      , TLINE_PARAMETERTYPE_BOOL               , false                      , false, 1},
			{"Substrate Thickness 1", TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness, true , 0},
			{"Substrate Material 1" , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material , false, 1},
			{"Substrate Thickness 2", TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness, true , 0},
			{"Substrate Material 2" , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material , false, 0},
		},
		{"Differential", "Common-mode"},
		&TLine_Stripline_Differential_Shielded,
	});

}
