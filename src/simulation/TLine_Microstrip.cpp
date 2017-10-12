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

#include <cfloat>

#include <iostream> // TODO: remove

void TLine_Microstrip_Single(TLineContext &context) {

	VDataDictReader root(context.m_parameters);

	real_t track_width = root.GetMember("track_width").AsFloat() * 1e-3;
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
	Box2D world_box = {
		track_box.x1 - space_x,
		0.0,
		track_box.x2 + space_x,
		track_box.y2 + space_y,
	};
	Box2D world_focus = {
		track_box.x1,
		0.0,
		track_box.x2,
		track_box.y2,
	};
	Box2D ground_box = {world_box.x1, 0.0, world_box.x2, 0.0};
	Box2D substrate_box = {world_box.x1, 0.0, world_box.x2, substrate_thickness};
	Box2D solder_mask_box1 = {substrate_box.x2, substrate_box.y2, substrate_box.x1, substrate_box.y2 + solder_mask_thickness_2};
	Box2D solder_mask_box2 = {track_box.x1 - solder_mask_thickness_1, track_box.y1, track_box.x2 + solder_mask_thickness_1, track_box.y2 + solder_mask_thickness_1};

	real_t step0 = REAL_MAX, step1 = substrate_thickness * 0.02;

	std::unique_ptr<GridMesh2D> mesh(new GridMesh2D(world_box, world_focus, 0.15, substrate_thickness * 1.0e-6));

	size_t port_ground = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);
	size_t port_signal = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);

	mesh->AddConductor(ground_box, step0, track_material, port_ground);
	mesh->AddConductor(track_box, step1, track_material, port_signal);
	mesh->AddDielectric(substrate_box, step0, substrate_material);
	mesh->AddDielectric(solder_mask_box1, step0, solder_mask_material);
	mesh->AddDielectric(solder_mask_box2, step0, solder_mask_material);

	context.m_output_mesh = std::move(mesh);
	TLineSolveModes(context, {0.0, 1.0}, {1.0});

}

void TLine_Microstrip_Differential(TLineContext &context) {

	VDataDictReader root(context.m_parameters);

	real_t track_width = root.GetMember("track_width").AsFloat() * 1e-3;
	real_t track_spacing = root.GetMember("track_spacing").AsFloat() * 1e-3;
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
	Box2D world_box = {
		track1_box.x1 - space_x,
		0.0,
		track2_box.x2 + space_x,
		track1_box.y2 + space_y,
	};
	Box2D world_focus = {
		track1_box.x1,
		0.0,
		track2_box.x2,
		track1_box.y2,
	};
	Box2D ground_box = {world_box.x1, 0.0, world_box.x2, 0.0};
	Box2D substrate_box = {world_box.x1, 0.0, world_box.x2, substrate_thickness};
	Box2D solder_mask_box1 = {substrate_box.x2, substrate_box.y2, substrate_box.x1, substrate_box.y2 + solder_mask_thickness_2};
	Box2D solder_mask_box2 = {track1_box.x1 - solder_mask_thickness_1, track1_box.y1, track1_box.x2 + solder_mask_thickness_1, track1_box.y2 + solder_mask_thickness_1};
	Box2D solder_mask_box3 = {track2_box.x1 - solder_mask_thickness_1, track2_box.y1, track2_box.x2 + solder_mask_thickness_1, track2_box.y2 + solder_mask_thickness_1};

	real_t step0 = REAL_MAX, step1 = std::min(track_spacing, substrate_thickness) * 0.02;

	std::unique_ptr<GridMesh2D> mesh(new GridMesh2D(world_box, world_focus, 0.15, substrate_thickness * 1.0e-6));

	size_t port_ground = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);
	size_t port_signal1 = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);
	size_t port_signal2 = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);

	mesh->AddConductor(ground_box, step0, track_material, port_ground);
	mesh->AddConductor(track1_box, step1, track_material, port_signal1);
	mesh->AddConductor(track2_box, step1, track_material, port_signal2);
	mesh->AddDielectric(substrate_box, step0, substrate_material);
	mesh->AddDielectric(solder_mask_box1, step0, solder_mask_material);
	mesh->AddDielectric(solder_mask_box2, step0, solder_mask_material);
	mesh->AddDielectric(solder_mask_box3, step0, solder_mask_material);

	context.m_output_mesh = std::move(mesh);
	TLineSolveModes(context, {0.0, 1.0, -1.0, 0.0, 1.0, 1.0}, {2.0, 1.0});

}

void RegisterTLine_Microstrip() {

	VData default_track_width = Json::FromString("1.0");
	VData default_track_spacing = Json::FromString("1.0");
	VData default_track_thickness = Json::FromString("0.035");
	VData default_substrate_thickness = Json::FromString("1.6");
	VData default_solder_mask_thickness_1 = Json::FromString("0.015");
	VData default_solder_mask_thickness_2 = Json::FromString("0.025");
	VData default_track_material = "Copper";
	VData default_substrate_material = "Isola DE104";
	VData default_solder_mask_material = "Solder Mask";

	g_tline_types.push_back(TLineType{
		"Microstrip (single)",
		"A single track above a ground plane. This is the simplest and most common PCB transmission line. "
		"Microstrips require very little space but are more susceptible to crosstalk than most other types of transmission lines.",
		{
			{"Track Width"            , TLINE_PARAMETERTYPE_REAL               , default_track_width            , true , 0},
			{"Track Thickness"        , TLINE_PARAMETERTYPE_REAL               , default_track_thickness        , true , 0},
			{"Track Material"         , TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR , default_track_material         , false, 1},
			{"Substrate Thickness"    , TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness    , true , 0},
			{"Substrate Material"     , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material     , false, 1},
			{"Solder Mask Thickness 1", TLINE_PARAMETERTYPE_REAL               , default_solder_mask_thickness_1, true , 0},
			{"Solder Mask Thickness 2", TLINE_PARAMETERTYPE_REAL               , default_solder_mask_thickness_2, true , 0},
			{"Solder Mask Material"   , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_solder_mask_material   , false, 0},
		},
		{"Single-ended"},
		&TLine_Microstrip_Single,
	});

	g_tline_types.push_back(TLineType{
		"Microstrip (differential)",
		"A differential pair above a ground plane. This is the simplest and most common differential PCB transmission line. "
		"Microstrips require very little space but are more susceptible to crosstalk than most other types of transmission lines.",
		{
			{"Track Width"            , TLINE_PARAMETERTYPE_REAL               , default_track_width            , true , 0},
			{"Track Spacing"          , TLINE_PARAMETERTYPE_REAL               , default_track_spacing          , true , 0},
			{"Track Thickness"        , TLINE_PARAMETERTYPE_REAL               , default_track_thickness        , true , 0},
			{"Track Material"         , TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR , default_track_material         , false, 1},
			{"Substrate Thickness"    , TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness    , true , 0},
			{"Substrate Material"     , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material     , false, 1},
			{"Solder Mask Thickness 1", TLINE_PARAMETERTYPE_REAL               , default_solder_mask_thickness_1, true , 0},
			{"Solder Mask Thickness 2", TLINE_PARAMETERTYPE_REAL               , default_solder_mask_thickness_2, true , 0},
			{"Solder Mask Material"   , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_solder_mask_material   , false, 0},
		},
		{"Differential", "Common-mode"},
		&TLine_Microstrip_Differential,
	});

}
