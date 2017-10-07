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

void TLine_CoplanarWaveguide_Single(TLineContext &context) {

	VDataPathDict root(context.m_parameters);

	real_t track_width = FloatUnscale(root.GetMember("track_width").AsFloat());
	real_t ground_spacing = FloatUnscale(root.GetMember("ground_spacing").AsFloat());
	real_t track_thickness = FloatUnscale(root.GetMember("track_thickness").AsFloat());
	const MaterialConductor *track_material = FindConductor(root, "track_material", context.m_material_database);
	real_t substrate_thickness = FloatUnscale(root.GetMember("substrate_thickness").AsFloat());
	const MaterialDielectric *substrate_material = FindDielectric(root, "substrate_material", context.m_material_database);

	real_t space_x = (track_width + track_thickness + substrate_thickness) * 15.0;
	real_t space_y = (track_width + track_thickness + substrate_thickness) * 25.0;
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

	real_t step0 = REAL_MAX, step1 = fmin(ground_spacing, substrate_thickness) * 0.02;

	std::unique_ptr<GridMesh2D> mesh(new GridMesh2D(world_box, world_focus, 0.15, substrate_thickness * 1.0e-6));

	size_t port_ground = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);
	size_t port_signal = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);

	mesh->AddConductor(ground_box, step0, track_material, port_ground);
	mesh->AddConductor(track_box, step1, track_material, port_signal);
	mesh->AddConductor(ground1_box, step0, step1, step1, step1, track_material, port_ground);
	mesh->AddConductor(ground2_box, step1, step0, step1, step1, track_material, port_ground);
	mesh->AddDielectric(substrate_box, step0, substrate_material);

	context.m_output_mesh = std::move(mesh);
	TLineSolveModes(context, {0.0, 1.0}, {1.0});

}

void TLine_CoplanarWaveguide_Differential(TLineContext &context) {

	VDataPathDict root(context.m_parameters);

	real_t track_width = FloatUnscale(root.GetMember("track_width").AsFloat());
	real_t track_spacing = FloatUnscale(root.GetMember("track_spacing").AsFloat());
	real_t ground_spacing = FloatUnscale(root.GetMember("ground_spacing").AsFloat());
	real_t track_thickness = FloatUnscale(root.GetMember("track_thickness").AsFloat());
	const MaterialConductor *track_material = FindConductor(root, "track_material", context.m_material_database);
	real_t substrate_thickness = FloatUnscale(root.GetMember("substrate_thickness").AsFloat());
	const MaterialDielectric *substrate_material = FindDielectric(root, "substrate_material", context.m_material_database);

	real_t space_x = (track_width + track_thickness + substrate_thickness) * 15.0;
	real_t space_y = (track_width + track_thickness + substrate_thickness) * 25.0;
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

	real_t step0 = REAL_MAX, step1 = fmin(fmin(track_spacing, ground_spacing), substrate_thickness) * 0.02;

	std::unique_ptr<GridMesh2D> mesh(new GridMesh2D(world_box, world_focus, 0.15, substrate_thickness * 1.0e-6));

	size_t port_ground = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);
	size_t port_signal1 = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);
	size_t port_signal2 = mesh->AddPort(GridMesh2D::PORTTYPE_FIXED);

	mesh->AddConductor(ground_box, step0, track_material, port_ground);
	mesh->AddConductor(track1_box, step1, track_material, port_signal1);
	mesh->AddConductor(track2_box, step1, track_material, port_signal2);
	mesh->AddConductor(ground1_box, step0, step1, step1, step1, track_material, port_ground);
	mesh->AddConductor(ground2_box, step1, step0, step1, step1, track_material, port_ground);
	mesh->AddDielectric(substrate_box, step0, substrate_material);

	context.m_output_mesh = std::move(mesh);
	TLineSolveModes(context, {0.0, 1.0, -1.0, 0.0, 1.0, 1.0}, {2.0, 1.0});

}

void RegisterTLine_CoplanarWaveguide() {

	VData default_track_width = Json::FromString("1.0");
	VData default_track_spacing = Json::FromString("1.0");
	VData default_ground_spacing = Json::FromString("1.0");
	VData default_track_thickness = Json::FromString("0.035");
	VData default_substrate_thickness = Json::FromString("1.6");
	VData default_track_material = "Copper";
	VData default_substrate_material = "Isola DE104";

	g_tline_types.push_back(TLineType{
		"Coplanar Waveguide (single)",
		"A single track above a ground plane, with ground planes on both sides. "
		"Coplanar waveguides are much less susceptible to crosstalk compared to microstrips. "
		"Ground vias should be placed on both sides of the track at regular intervals (less than 1/10th of the wavelength) to ensure correct behavior. "
		"Isolation can be improved by adding more vias.",
		{
			{"Track Width"        , TLINE_PARAMETERTYPE_REAL               , default_track_width        , true , 0},
			{"Ground Spacing"     , TLINE_PARAMETERTYPE_REAL               , default_ground_spacing     , true , 0},
			{"Track Thickness"    , TLINE_PARAMETERTYPE_REAL               , default_track_thickness    , true , 0},
			{"Track Material"     , TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR , default_track_material     , false, 1},
			{"Substrate Thickness", TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness, true , 0},
			{"Substrate Material" , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material , false, 0},
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
			{"Track Width"        , TLINE_PARAMETERTYPE_REAL               , default_track_width        , true , 0},
			{"Track Spacing"      , TLINE_PARAMETERTYPE_REAL               , default_track_spacing      , true , 0},
			{"Ground Spacing"     , TLINE_PARAMETERTYPE_REAL               , default_ground_spacing     , true , 0},
			{"Track Thickness"    , TLINE_PARAMETERTYPE_REAL               , default_track_thickness    , true , 0},
			{"Track Material"     , TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR , default_track_material     , false, 1},
			{"Substrate Thickness", TLINE_PARAMETERTYPE_REAL               , default_substrate_thickness, true , 0},
			{"Substrate Material" , TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC, default_substrate_material , false, 0},
		},
		{"Differential", "Common-mode"},
		&TLine_CoplanarWaveguide_Differential,
	});

}
