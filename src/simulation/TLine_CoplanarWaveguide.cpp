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

void TLine_CoplanarWaveguide_Single(TLineContext &context) {

	VDataDictReader root(context.m_parameters);

	real_t track_width = root.GetMember("track_width").AsFloat() * 1e-3;
	real_t ground_spacing = root.GetMember("ground_spacing").AsFloat() * 1e-3;
	real_t track_thickness = root.GetMember("track_thickness").AsFloat() * 1e-3;
	const MaterialConductor *track_material = FindConductor(root, "track_material", context.m_material_database);
	real_t substrate_thickness = root.GetMember("substrate_thickness").AsFloat() * 1e-3;
	const MaterialDielectric *substrate_material = FindDielectric(root, "substrate_material", context.m_material_database);
	real_t solder_mask_thickness_1 = root.GetMember("solder_mask_thickness_1").AsFloat() * 1e-3;
	real_t solder_mask_thickness_2 = root.GetMember("solder_mask_thickness_2").AsFloat() * 1e-3;
	const MaterialDielectric *solder_mask_material = FindDielectric(root, "solder_mask_material", context.m_material_database);

	real_t space_x = (track_width + track_thickness + substrate_thickness + std::max(solder_mask_thickness_1, solder_mask_thickness_2)) * 15.0;
	real_t space_y = (track_width + track_thickness + substrate_thickness + std::max(solder_mask_thickness_1, solder_mask_thickness_2)) * 25.0;
	Box2D track_box = {
		-0.5 * track_width,
		substrate_thickness,
		0.5 * track_width,
		substrate_thickness + track_thickness,
	};
	Box2D ground1_box = {
		track_box.x1 - ground_spacing - space_x,
		substrate_thickness,
		track_box.x1 - ground_spacing,
		substrate_thickness + track_thickness,
	};
	Box2D ground2_box = ground1_box.MirroredX();
	Box2D world_box = {
		ground1_box.x1,
		0.0,
		ground2_box.x2,
		track_box.y2 + space_y,
	};
	Box2D world_focus = {
		ground1_box.x2,
		0.0,
		ground2_box.x1,
		track_box.y2,
	};
	Box2D ground_box = {world_box.x1, 0.0, world_box.x2, 0.0};
	Box2D substrate_box = {world_box.x1, 0.0, world_box.x2, substrate_thickness};
	Box2D solder_mask_box1 = {substrate_box.x2, substrate_box.y2, substrate_box.x1, substrate_box.y2 + solder_mask_thickness_2};
	Box2D solder_mask_box2 = {track_box.x1 - solder_mask_thickness_1, track_box.y1, track_box.x2 + solder_mask_thickness_1, track_box.y2 + solder_mask_thickness_1};
	Box2D solder_mask_box3 = {ground1_box.x1, ground1_box.y1, ground1_box.x2 + solder_mask_thickness_1, ground1_box.y2 + solder_mask_thickness_1};
	Box2D solder_mask_box4 = {ground2_box.x1 - solder_mask_thickness_1, ground2_box.y1, ground2_box.x2, ground2_box.y2 + solder_mask_thickness_1};

	real_t step0 = REAL_MAX, step1 = std::min(ground_spacing, substrate_thickness) * GridMesh2D::DEFAULT_GRID_STEP;

	std::unique_ptr<GridMesh2D> mesh(new GridMesh2D(world_box, world_focus, GridMesh2D::DEFAULT_GRID_INC, substrate_thickness * 1.0e-6));

	size_t port_ground = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);
	size_t port_signal = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);

	mesh->AddConductor(ground_box, step0, track_material, port_ground);
	mesh->AddConductor(track_box, step1, track_material, port_signal);
	mesh->AddConductor(ground1_box, step0, step1, step1, step1, track_material, port_ground);
	mesh->AddConductor(ground2_box, step1, step0, step1, step1, track_material, port_ground);
	mesh->AddDielectric(substrate_box, step0, substrate_material);
	mesh->AddDielectric(solder_mask_box1, step0, solder_mask_material);
	mesh->AddDielectric(solder_mask_box2, step0, solder_mask_material);
	mesh->AddDielectric(solder_mask_box3, step0, solder_mask_material);
	mesh->AddDielectric(solder_mask_box4, step0, solder_mask_material);

	context.m_output_mesh = std::move(mesh);

	Eigen::MatrixXr modes(2, 1);
	modes.col(0) << 0.0, 1.0;
	TLineSolveModes(context, modes);

}

void TLine_CoplanarWaveguide_Differential(TLineContext &context) {

	VDataDictReader root(context.m_parameters);

	real_t track_width = root.GetMember("track_width").AsFloat() * 1e-3;
	real_t track_spacing = root.GetMember("track_spacing").AsFloat() * 1e-3;
	real_t ground_spacing = root.GetMember("ground_spacing").AsFloat() * 1e-3;
	real_t track_thickness = root.GetMember("track_thickness").AsFloat() * 1e-3;
	const MaterialConductor *track_material = FindConductor(root, "track_material", context.m_material_database);
	real_t substrate_thickness = root.GetMember("substrate_thickness").AsFloat() * 1e-3;
	const MaterialDielectric *substrate_material = FindDielectric(root, "substrate_material", context.m_material_database);
	real_t solder_mask_thickness_1 = root.GetMember("solder_mask_thickness_1").AsFloat() * 1e-3;
	real_t solder_mask_thickness_2 = root.GetMember("solder_mask_thickness_2").AsFloat() * 1e-3;
	const MaterialDielectric *solder_mask_material = FindDielectric(root, "solder_mask_material", context.m_material_database);

	real_t space_x = (track_width * 2 + track_spacing + track_thickness + substrate_thickness + std::max(solder_mask_thickness_1, solder_mask_thickness_2)) * 15.0;
	real_t space_y = (track_width * 2 + track_spacing + track_thickness + substrate_thickness + std::max(solder_mask_thickness_1, solder_mask_thickness_2)) * 25.0;
	Box2D track1_box = {
		-0.5 * track_spacing - track_width,
		substrate_thickness,
		-0.5 * track_spacing,
		substrate_thickness + track_thickness,
	};
	Box2D track2_box = track1_box.MirroredX();
	Box2D ground1_box = {
		track1_box.x1 - ground_spacing - space_x,
		substrate_thickness,
		track1_box.x1 - ground_spacing,
		substrate_thickness + track_thickness,
	};
	Box2D ground2_box = ground1_box.MirroredX();
	Box2D world_box = {
		ground1_box.x1,
		0.0,
		ground2_box.x2,
		track1_box.y2 + space_y,
	};
	Box2D world_focus = {
		ground1_box.x2,
		0.0,
		ground2_box.x1,
		track1_box.y2,
	};
	Box2D ground_box = {world_box.x1, 0.0, world_box.x2, 0.0};
	Box2D substrate_box = {world_box.x1, 0.0, world_box.x2, substrate_thickness};
	Box2D solder_mask_box1 = {substrate_box.x2, substrate_box.y2, substrate_box.x1, substrate_box.y2 + solder_mask_thickness_2};
	Box2D solder_mask_box2 = {track1_box.x1 - solder_mask_thickness_1, track1_box.y1, track1_box.x2 + solder_mask_thickness_1, track1_box.y2 + solder_mask_thickness_1};
	Box2D solder_mask_box3 = {track2_box.x1 - solder_mask_thickness_1, track2_box.y1, track2_box.x2 + solder_mask_thickness_1, track2_box.y2 + solder_mask_thickness_1};
	Box2D solder_mask_box4 = {ground1_box.x1, ground1_box.y1, ground1_box.x2 + solder_mask_thickness_1, ground1_box.y2 + solder_mask_thickness_1};
	Box2D solder_mask_box5 = {ground2_box.x1 - solder_mask_thickness_1, ground2_box.y1, ground2_box.x2, ground2_box.y2 + solder_mask_thickness_1};

	real_t step0 = REAL_MAX, step1 = std::min(std::min(track_spacing, ground_spacing), substrate_thickness) * GridMesh2D::DEFAULT_GRID_STEP;

	std::unique_ptr<GridMesh2D> mesh(new GridMesh2D(world_box, world_focus, GridMesh2D::DEFAULT_GRID_INC, substrate_thickness * 1.0e-6));

	size_t port_ground = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);
	size_t port_signal1 = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);
	size_t port_signal2 = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);

	mesh->AddConductor(ground_box, step0, track_material, port_ground);
	mesh->AddConductor(track1_box, step1, track_material, port_signal1);
	mesh->AddConductor(track2_box, step1, track_material, port_signal2);
	mesh->AddConductor(ground1_box, step0, step1, step1, step1, track_material, port_ground);
	mesh->AddConductor(ground2_box, step1, step0, step1, step1, track_material, port_ground);
	mesh->AddDielectric(substrate_box, step0, substrate_material);
	mesh->AddDielectric(solder_mask_box1, step0, solder_mask_material);
	mesh->AddDielectric(solder_mask_box2, step0, solder_mask_material);
	mesh->AddDielectric(solder_mask_box3, step0, solder_mask_material);
	mesh->AddDielectric(solder_mask_box4, step0, solder_mask_material);
	mesh->AddDielectric(solder_mask_box5, step0, solder_mask_material);

	context.m_output_mesh = std::move(mesh);

	Eigen::MatrixXr modes(3, 2);
	modes.col(0) << 0.0, 0.5, -0.5;
	modes.col(1) << 0.0, 1.0, 1.0;
	TLineSolveModes(context, modes);

}

void RegisterTLine_CoplanarWaveguide() {

	VData default_track_width = Json::FromString("1.0");
	VData default_track_spacing = Json::FromString("1.0");
	VData default_ground_spacing = Json::FromString("1.0");
	VData default_track_thickness = Json::FromString("0.035");
	VData default_substrate_thickness = Json::FromString("1.6");
	VData default_solder_mask_thickness_1 = Json::FromString("0.015");
	VData default_solder_mask_thickness_2 = Json::FromString("0.025");
	VData default_track_material = "Copper";
	VData default_substrate_material = "Isola DE104";
	VData default_solder_mask_material = "Solder Mask";

	g_tline_types.push_back(TLineType{
		"Coplanar Waveguide (single)",
		"A single track above a ground plane, with ground planes on both sides. "
		"Coplanar waveguides are much less susceptible to crosstalk compared to microstrips. "
		"Ground vias should be placed on both sides of the track at regular intervals (less than 1/10th of the wavelength) to ensure correct behavior. "
		"Isolation can be improved by adding more vias.",
		{
			{"Track Width"            , TLINE_PARAMETERTYPE_REAL               , default_track_width            , true , 0},
			{"Ground Spacing"         , TLINE_PARAMETERTYPE_REAL               , default_ground_spacing         , true , 0},
			{"Track Thickness"        , TLINE_PARAMETERTYPE_REAL               , default_track_thickness        , true , 0},
			{"Track Material"         , TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR , default_track_material         , false, 1},
			{"Substrate Thickness"    , TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness    , true , 0},
			{"Substrate Material"     , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material     , false, 1},
			{"Solder Mask Thickness 1", TLINE_PARAMETERTYPE_REAL               , default_solder_mask_thickness_1, true , 0},
			{"Solder Mask Thickness 2", TLINE_PARAMETERTYPE_REAL               , default_solder_mask_thickness_2, true , 0},
			{"Solder Mask Material"   , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_solder_mask_material   , false, 0},
		},
		{"Single-ended"},
		&TLine_CoplanarWaveguide_Single,
	});

	g_tline_types.push_back(TLineType{
		"Coplanar Waveguide (differential)",
		"A differential pair above a ground plane, with ground planes on both sides. "
		"Coplanar waveguides are much less susceptible to crosstalk compared to microstrips. "
		"Ground vias should be placed on both sides of the tracks at regular intervals (less than 1/10th of the wavelength) to ensure correct behavior. "
		"Isolation can be improved by adding more vias.",
		{
			{"Track Width"            , TLINE_PARAMETERTYPE_REAL               , default_track_width            , true , 0},
			{"Track Spacing"          , TLINE_PARAMETERTYPE_REAL               , default_track_spacing          , true , 0},
			{"Ground Spacing"         , TLINE_PARAMETERTYPE_REAL               , default_ground_spacing         , true , 0},
			{"Track Thickness"        , TLINE_PARAMETERTYPE_REAL               , default_track_thickness        , true , 0},
			{"Track Material"         , TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR , default_track_material         , false, 1},
			{"Substrate Thickness"    , TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness    , true , 0},
			{"Substrate Material"     , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material     , false, 1},
			{"Solder Mask Thickness 1", TLINE_PARAMETERTYPE_REAL               , default_solder_mask_thickness_1, true , 0},
			{"Solder Mask Thickness 2", TLINE_PARAMETERTYPE_REAL               , default_solder_mask_thickness_2, true , 0},
			{"Solder Mask Material"   , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_solder_mask_material   , false, 0},
		},
		{"Differential", "Common-mode"},
		&TLine_CoplanarWaveguide_Differential,
	});

}
