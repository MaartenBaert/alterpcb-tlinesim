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

#include "AboutDialog.h"
#include "Basics.h"
#include "FixedScrollArea.h"
#include "VData.h"
#include "Qt.h"

class MaterialDatabase;
class MeshViewer;
class TLineContext;

class MainWindow : public QMainWindow {
	Q_OBJECT

private:
	enum SimulationType {
		SIMULATION_SINGLE_FREQUENCY,
		SIMULATION_FREQUENCY_SWEEP,
		SIMULATION_PARAMETER_SWEEP,
		SIMULATION_PARAMETER_TUNE,
		SIMULATION_COUNT, // must be last
	};

	enum MeshDetail {
		MESHDETAIL_VERYLOW,
		MESHDETAIL_LOWER,
		MESHDETAIL_LOW,
		MESHDETAIL_MEDIUM,
		MESHDETAIL_HIGH,
		MESHDETAIL_HIGHER,
		MESHDETAIL_VERYHIGH,
		MESHDETAIL_COUNT, // must be last
	};

private:
	static const QString WINDOW_CAPTION;

private:
	std::unique_ptr<MaterialDatabase> m_material_database;
	size_t m_tline_type;

	QComboBox *m_combobox_tline_types;
	QPlainTextEdit *m_textedit_description;

	FixedScrollArea *m_scrollarea_parameters;
	std::vector<QWidget*> m_widget_parameters;

	QComboBox *m_combobox_simulation_type;
	QPushButton *m_pushbutton_simulate;

	QLabel *m_label_frequency[2];
	QLineEdit *m_lineedit_frequency;

	QLabel *m_label_frequency_sweep[4];
	QLineEdit *m_lineedit_frequency_sweep_min, *m_lineedit_frequency_sweep_max, *m_lineedit_frequency_sweep_step;
	QLabel *m_label_frequency_sweep_file;
	QLineEdit *m_lineedit_frequency_sweep_file;
	QPushButton *m_pushbutton_frequency_sweep_browse;

	QLabel *m_label_parameter_sweep[4];
	QComboBox *m_combobox_parameter_sweep_parameter;
	QLineEdit *m_lineedit_parameter_sweep_min, *m_lineedit_parameter_sweep_max, *m_lineedit_parameter_sweep_step;
	QLabel *m_label_parameter_sweep_file;
	QLineEdit *m_lineedit_parameter_sweep_file;
	QPushButton *m_pushbutton_parameter_sweep_browse;

	QLabel *m_label_parameter_tune[3];
	QComboBox *m_combobox_parameter_tune_parameter, *m_combobox_parameter_tune_target_result;
	QLineEdit *m_lineedit_parameter_tune_target_value;

	QComboBox *m_combobox_solver_type;
	QComboBox *m_combobox_element_type;
	QComboBox *m_combobox_mesh_detail;

	FixedScrollArea *m_scrollarea_results;
	std::vector<QLineEdit*> m_lineedit_results;

	MeshViewer *m_meshviewer;
	QSlider *m_slider_zoom;
	QComboBox *m_combobox_image_type;
	QCheckBox *m_checkbox_mesh_overlay;
	QComboBox *m_combobox_modes;

public:
	MainWindow();
	~MainWindow();

private:
	void LoadMaterials();
	void SimulationInit(TLineContext &context);
	void SimulationShowResult(TLineContext &context);
	void SimulateSingleFrequency();
	void SimulateFrequencySweep();
	void SimulateParameterSweep();
	void SimulateParameterTune();

	void ProcessSlowEvents(int msec = 10);

private slots:
	void OnUpdateTLineType();
	void OnUpdateSimulationType();
	void OnFrequencySweepBrowse();
	void OnParameterSweepBrowse();
	void OnSimulate();

	void OnZoomChange();
	void OnImageTypeChange();
	void OnMeshOverlayChange();
	void OnModeChange();

	void OnAbout();

private:
	inline SimulationType GetSimulationType() { return (SimulationType) clamp(m_combobox_simulation_type->currentIndex(), 0, SIMULATION_COUNT - 1); }

};
