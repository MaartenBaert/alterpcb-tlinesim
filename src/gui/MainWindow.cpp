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

#include "MainWindow.h"

#include "ApplicationDirs.h"
#include "FindRoot.h"
#include "Icons.h"
#include "Json.h"
#include "LayoutHelper.h"
#include "MaterialDatabase.h"
#include "MeshViewer.h"
#include "QLineEditSmall.h"
#include "TLineTypes.h"

#include <fstream>
#include <iostream>

const QString MainWindow::WINDOW_CAPTION = "AlterPCB Transmission Line Simulator";

static real_t FloatFromVData(const VData &data) {
	if(data.GetType() == VDATA_INT)
		return (real_t) data.AsInt(); // implicit int to float conversion
	if(data.GetType() == VDATA_FLOAT)
		return FloatUnscale(data.AsFloat());
	throw std::runtime_error(MakeString("Expected float, got ", EnumToString(data.GetType()), " instead."));
}

static real_t FloatFromString(const std::string &str) {
	return FloatFromVData(Json::FromString(str));
}

static void MakeSweep(std::vector<real_t> &results, real_t min, real_t max, real_t step) {
	index_t num = std::max<index_t>(lrint(fmax(0.0, max - min) / step + 0.5 + 1e-12), 1);
	results.clear();
	results.resize(num);
	for(index_t i = 0; i < num; ++i) {
		results[i] = min + step * (real_t) i;
	}
}

MainWindow::MainWindow() {

	m_material_database.reset(new MaterialDatabase());
	m_tline_type = 0;

	setWindowTitle(WINDOW_CAPTION);

	QWidget *centralwidget = new QWidget(this);
	setCentralWidget(centralwidget);

	QGroupBox *groupbox_type = new QGroupBox("Transmission line type", this);
	{
		m_combobox_tline_types = new QComboBox(groupbox_type);
		for(const TLineType &type : g_tline_types) {
			m_combobox_tline_types->addItem(QString::fromStdString(type.m_name));
		}
		m_textedit_description = new QPlainTextEdit(groupbox_type);
		m_textedit_description->setFixedHeight(80);
		m_textedit_description->setReadOnly(true);

		connect(m_combobox_tline_types, SIGNAL(currentIndexChanged(int)), this, SLOT(OnUpdateTLineType()));

		QVBoxLayout *layout = new QVBoxLayout(groupbox_type);
		layout->addWidget(m_combobox_tline_types);
		layout->addWidget(m_textedit_description);
	}
	QGroupBox *groupbox_parameters = new QGroupBox("Parameters", this);
	{
		m_scrollarea_parameters = new QScrollArea(groupbox_parameters);
		m_scrollarea_parameters->setWidgetResizable(true);
		m_scrollarea_parameters->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
		m_scrollarea_parameters->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

		QVBoxLayout *layout = new QVBoxLayout(groupbox_parameters);
		layout->addWidget(m_scrollarea_parameters);
	}
	QGroupBox *groupbox_simulation = new QGroupBox("Simulation", this);
	{
		QLabel *label_simulation_type = new QLabel("Simulation type:", groupbox_simulation);
		m_combobox_simulation_type = new QComboBox(groupbox_simulation);
		m_combobox_simulation_type->addItem("Single Frequency");
		m_combobox_simulation_type->addItem("Frequency Sweep");
		m_combobox_simulation_type->addItem("Parameter Sweep");
		m_combobox_simulation_type->addItem("Parameter Tune");
		m_pushbutton_simulate = new QPushButton(g_icon_simulation, "Simulate", this);
		m_pushbutton_simulate->setIconSize(QSize(16, 16));

		m_label_frequency[0] = new QLabel("Frequency:", groupbox_simulation);
		m_lineedit_frequency = new QLineEdit("1", groupbox_simulation);
		m_label_frequency[1] = new QLabel("GHz", groupbox_simulation);

		m_label_frequency_sweep[0] = new QLabel("Frequency:", groupbox_simulation);
		m_lineedit_frequency_sweep_min = new QLineEditSmall("0.0", groupbox_simulation, 40);
		m_label_frequency_sweep[1] = new QLabel("to", groupbox_simulation);
		m_lineedit_frequency_sweep_max = new QLineEditSmall("10.0", groupbox_simulation, 40);
		m_label_frequency_sweep[2] = new QLabel("in steps of", groupbox_simulation);
		m_lineedit_frequency_sweep_step = new QLineEditSmall("1.0", groupbox_simulation, 40);
		m_label_frequency_sweep[3] = new QLabel("GHz", groupbox_simulation);
		m_label_frequency_sweep_file = new QLabel("Output File:", groupbox_simulation);
		m_lineedit_frequency_sweep_file = new QLineEdit(QDir::homePath() + "/frequency_sweep.txt", groupbox_simulation);
		m_pushbutton_frequency_sweep_browse = new QPushButton("Browse...", groupbox_simulation);

		m_label_parameter_sweep[0] = new QLabel("Parameter:", groupbox_simulation);
		m_combobox_parameter_sweep_parameter = new QComboBox(groupbox_simulation);
		m_label_parameter_sweep[1] = new QLabel("Value:", groupbox_simulation);
		m_lineedit_parameter_sweep_min = new QLineEditSmall("1.0", groupbox_simulation, 40);
		m_label_parameter_sweep[2] = new QLabel("to", groupbox_simulation);
		m_lineedit_parameter_sweep_max = new QLineEditSmall("2.0", groupbox_simulation, 40);
		m_label_parameter_sweep[3] = new QLabel("in steps of", groupbox_simulation);
		m_lineedit_parameter_sweep_step = new QLineEditSmall("0.1", groupbox_simulation, 40);
		m_label_parameter_sweep_file = new QLabel("Output File:", groupbox_simulation);
		m_lineedit_parameter_sweep_file = new QLineEdit(QDir::homePath() + "/parameter_sweep.txt", groupbox_simulation);
		m_pushbutton_parameter_sweep_browse = new QPushButton("Browse...", groupbox_simulation);

		m_label_parameter_tune[0] = new QLabel("Parameter:", groupbox_simulation);
		m_combobox_parameter_tune_parameter = new QComboBox(groupbox_simulation);
		m_label_parameter_tune[1] = new QLabel("Target Result:", groupbox_simulation);
		m_combobox_parameter_tune_target_result = new QComboBox(groupbox_simulation);
		m_label_parameter_tune[2] = new QLabel("Target Value:", groupbox_simulation);
		m_lineedit_parameter_tune_target_value = new QLineEdit("50.0", groupbox_simulation);

		connect(m_combobox_simulation_type, SIGNAL(currentIndexChanged(int)), this, SLOT(OnUpdateSimulationType()));
		connect(m_pushbutton_simulate, SIGNAL(clicked(bool)), this, SLOT(OnSimulate()));
		connect(m_combobox_parameter_sweep_parameter, SIGNAL(currentIndexChanged(int)), this, SLOT(OnUpdateSimulationType()));
		connect(m_combobox_parameter_tune_parameter, SIGNAL(currentIndexChanged(int)), this, SLOT(OnUpdateSimulationType()));
		connect(m_pushbutton_frequency_sweep_browse, SIGNAL(clicked(bool)), this, SLOT(OnFrequencySweepBrowse()));
		connect(m_pushbutton_parameter_sweep_browse, SIGNAL(clicked(bool)), this, SLOT(OnParameterSweepBrowse()));

		QVBoxLayout *layout = new QVBoxLayout(groupbox_simulation);
		{
			QHBoxLayout *layout2 = new QHBoxLayout();
			layout->addLayout(layout2);
			layout2->addWidget(label_simulation_type);
			layout2->addWidget(m_combobox_simulation_type);
			layout2->addWidget(m_pushbutton_simulate);
		}
		{
			QGridLayout *layout2 = new QGridLayout();
			layout->addLayout(layout2);
			layout2->addWidget(m_label_frequency[0], 0, 0);
			{
				QHBoxLayout *layout3 = new QHBoxLayout();
				layout2->addLayout(layout3, 0, 1);
				layout3->addWidget(m_lineedit_frequency);
				layout3->addWidget(m_label_frequency[1]);
			}
			layout2->addWidget(m_label_frequency_sweep[0], 1, 0);
			{
				QHBoxLayout *layout3 = new QHBoxLayout();
				layout2->addLayout(layout3, 1, 1);
				layout3->addWidget(m_lineedit_frequency_sweep_min);
				layout3->addWidget(m_label_frequency_sweep[1]);
				layout3->addWidget(m_lineedit_frequency_sweep_max);
				layout3->addWidget(m_label_frequency_sweep[2]);
				layout3->addWidget(m_lineedit_frequency_sweep_step);
				layout3->addWidget(m_label_frequency_sweep[3]);
			}
			layout2->addWidget(m_label_frequency_sweep_file, 2, 0);
			{
				QHBoxLayout *layout3 = new QHBoxLayout();
				layout2->addLayout(layout3, 2, 1);
				layout3->addWidget(m_lineedit_frequency_sweep_file);
				layout3->addWidget(m_pushbutton_frequency_sweep_browse);
			}
			layout2->addWidget(m_label_parameter_sweep[0], 3, 0);
			layout2->addWidget(m_combobox_parameter_sweep_parameter, 3, 1);
			layout2->addWidget(m_label_parameter_sweep[1], 4, 0);
			{
				QHBoxLayout *layout3 = new QHBoxLayout();
				layout2->addLayout(layout3, 4, 1);
				layout3->addWidget(m_lineedit_parameter_sweep_min);
				layout3->addWidget(m_label_parameter_sweep[2]);
				layout3->addWidget(m_lineedit_parameter_sweep_max);
				layout3->addWidget(m_label_parameter_sweep[3]);
				layout3->addWidget(m_lineedit_parameter_sweep_step);
			}
			layout2->addWidget(m_label_parameter_sweep_file, 5, 0);
			{
				QHBoxLayout *layout3 = new QHBoxLayout();
				layout2->addLayout(layout3, 5, 1);
				layout3->addWidget(m_lineedit_parameter_sweep_file);
				layout3->addWidget(m_pushbutton_parameter_sweep_browse);
			}
			layout2->addWidget(m_label_parameter_tune[0], 6, 0);
			layout2->addWidget(m_combobox_parameter_tune_parameter, 6, 1);
			layout2->addWidget(m_label_parameter_tune[1], 7, 0);
			layout2->addWidget(m_combobox_parameter_tune_target_result, 7, 1);
			layout2->addWidget(m_label_parameter_tune[2], 8, 0);
			layout2->addWidget(m_lineedit_parameter_tune_target_value, 8, 1);
		}

	}
	QGroupBox *groupbox_results = new QGroupBox("Results", this);
	{
		m_scrollarea_results = new QScrollArea(groupbox_results);
		m_scrollarea_results->setWidgetResizable(true);
		m_scrollarea_results->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
		m_scrollarea_results->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

		QVBoxLayout *layout = new QVBoxLayout(groupbox_results);
		layout->addWidget(m_scrollarea_results);
	}
	QGroupBox *groupbox_viewer = new QGroupBox("Viewer", this);
	{
		m_meshviewer = new MeshViewer(groupbox_viewer);
		QLabel *label_zoom = new QLabel("Zoom:", groupbox_viewer);
		m_slider_zoom = new QSlider(Qt::Horizontal, groupbox_viewer);
		m_slider_zoom->setRange(0, 2000000);
		m_slider_zoom->setValue(1500000);
		m_slider_zoom->setSingleStep(10000);
		m_slider_zoom->setPageStep(100000);
		QLabel *label_imagetype = new QLabel("View type:", groupbox_viewer);
		m_buttongroup_image_type = new QButtonGroup(groupbox_viewer);
		QRadioButton *radio_imagetype_mesh = new QRadioButton("Mesh", groupbox_viewer);
		QRadioButton *radio_imagetype_field_e = new QRadioButton("E-field", groupbox_viewer);
		QRadioButton *radio_imagetype_field_h = new QRadioButton("H-field", groupbox_viewer);
		m_buttongroup_image_type->addButton(radio_imagetype_mesh, MESHIMAGETYPE_MESH);
		m_buttongroup_image_type->addButton(radio_imagetype_field_e, MESHIMAGETYPE_FIELD_E);
		m_buttongroup_image_type->addButton(radio_imagetype_field_h, MESHIMAGETYPE_FIELD_H);
		radio_imagetype_mesh->setChecked(true);
		m_checkbox_mesh_overlay = new QCheckBox("Mesh Overlay", groupbox_viewer);
		QLabel *label_mode = new QLabel("Mode:", groupbox_viewer);
		m_combobox_modes = new QComboBox(groupbox_viewer);
		m_combobox_modes->setMinimumWidth(120);

		connect(m_slider_zoom, SIGNAL(valueChanged(int)), this, SLOT(OnZoomChange()));
		connect(m_buttongroup_image_type, SIGNAL(buttonClicked(int)), this, SLOT(OnImageTypeChange()));
		connect(m_checkbox_mesh_overlay, SIGNAL(stateChanged(int)), this, SLOT(OnMeshOverlayChange()));
		connect(m_combobox_modes, SIGNAL(currentIndexChanged(int)), this, SLOT(OnModeChange()));

		QVBoxLayout *layout = new QVBoxLayout(groupbox_viewer);
		layout->addWidget(m_meshviewer);
		{
			QHBoxLayout *layout2 = new QHBoxLayout();
			layout->addLayout(layout2);
			layout2->addWidget(label_zoom);
			layout2->addWidget(m_slider_zoom);
		}
		{
			QHBoxLayout *layout2 = new QHBoxLayout();
			layout->addLayout(layout2);
			layout2->addWidget(label_imagetype);
			layout2->addWidget(radio_imagetype_mesh);
			layout2->addWidget(radio_imagetype_field_e);
			layout2->addWidget(radio_imagetype_field_h);
			layout2->addStretch();
			layout2->addWidget(m_checkbox_mesh_overlay);
			layout2->addWidget(label_mode);
			layout2->addWidget(m_combobox_modes);
		}
	}

	QHBoxLayout *layout = new QHBoxLayout(centralwidget);
	layout->addWidget(groupbox_viewer, 1);
	{
		QVBoxLayout *layout2 = new QVBoxLayout();
		layout->addLayout(layout2);
		layout2->addWidget(groupbox_type);
		layout2->addWidget(groupbox_parameters, 1);
		layout2->addWidget(groupbox_simulation);
		layout2->addWidget(groupbox_results, 1);
	}

	LoadMaterials();
	OnUpdateTLineType();
	OnUpdateSimulationType();

	OnZoomChange();
	OnImageTypeChange();
	OnMeshOverlayChange();
	OnModeChange();

	show();

}

MainWindow::~MainWindow() {
	// nothing
}

void MainWindow::LoadMaterials() {
	if(g_application_data_dir.isEmpty()) {
		std::cerr << "Error: Could not load materials, data directory is missing." << std::endl;
		return;
	}
	try {
		m_material_database->LoadFile(g_application_data_dir.toStdString() + "/materials.json");
		m_material_database->Finish();
	}
	catch(const std::runtime_error &e) {
		std::cerr << "Error: Could not load material database:" << std::endl;
		std::cerr << "    " << e.what() << std::endl;
	}
}

void MainWindow::GetParameterValues(VData::Dict &result) {
	const TLineType &tline_type = g_tline_types[m_tline_type];
	for(size_t i = 0; i < tline_type.m_parameters.size(); ++i) {
		const TLineParameter &parameter = tline_type.m_parameters[i];
		stringtag_t key = StringRegistry::NewTag(CanonicalName(parameter.m_name));
		VData value;
		switch(parameter.m_type) {
			case TLINE_PARAMETERTYPE_BOOL: {
				QCheckBox *checkbox = static_cast<QCheckBox*>(m_widget_parameters[i]);
				value = checkbox->isChecked();
				break;
			}
			case TLINE_PARAMETERTYPE_REAL: {
				QLineEdit *lineedit_value = static_cast<QLineEdit*>(m_widget_parameters[i]);
				Json::FromString(value, lineedit_value->text().toStdString());
				break;
			}
			case TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR:
			case TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC: {
				QComboBox *combobox_value = static_cast<QComboBox*>(m_widget_parameters[i]);
				value = combobox_value->currentText().toStdString();
				break;
			}
			case TLINE_PARAMETERTYPE_COUNT: {
				assert(false);
				break;
			}
		}
		result.EmplaceBack(VDataDictEntry(key, std::move(value)));
	}
}

void MainWindow::OnUpdateTLineType() {

	m_tline_type = clamp<index_t>(m_combobox_tline_types->currentIndex(), 0, g_tline_types.size() - 1);
	const TLineType &tline_type = g_tline_types[m_tline_type];
	const std::vector<MaterialConductor> &conductors = m_material_database->GetConductors();
	const std::vector<MaterialDielectric> &dielectrics = m_material_database->GetDielectrics();

	// update description
	m_textedit_description->setPlainText(QString::fromStdString(tline_type.m_description));

	/*for(ParameterWidgets &p : m_parameter_widgets) {
		delete p->m_label_name;
		delete p->m_lineedit;
		delete p->m_label_unit;
	}*/

	// update parameters
	m_widget_parameters.clear();
	m_combobox_parameter_sweep_parameter->clear();
	m_combobox_parameter_tune_parameter->clear();
	delete m_scrollarea_parameters->widget();
	{
		QWidget *widget = new QWidget(m_scrollarea_parameters), *widget_focus = m_scrollarea_parameters;
		QGridLayout *layout = new QGridLayout(widget);
		size_t row = 0;
		for(index_t i = 0; i < tline_type.m_parameters.size(); ++i) {
			const TLineParameter &parameter = tline_type.m_parameters[i];
			QLabel *label_name = new QLabel(QString::fromStdString(parameter.m_name) + ":", widget);
			QWidget *widget_value = NULL;
			QLabel *label_unit = NULL;
			switch(parameter.m_type) {
				case TLINE_PARAMETERTYPE_BOOL: {
					QCheckBox *checkbox_value = new QCheckBox(widget);
					checkbox_value->setChecked(parameter.m_default_value.AsBool());
					widget_value = checkbox_value;
					break;
				}
				case TLINE_PARAMETERTYPE_REAL: {
					QLineEdit *lineedit_value = new QLineEdit(widget);
					lineedit_value->setText(QString::fromStdString(Json::ToString(parameter.m_default_value)));
					widget_value = lineedit_value;
					label_unit = new QLabel((parameter.m_unit_mm)? "mm" : "", widget);
					break;
				}
				case TLINE_PARAMETERTYPE_MATERIAL_CONDUCTOR: {
					QComboBox *combobox_value = new QComboBox(widget);
					for(size_t j = 0; j < conductors.size(); ++j) {
						combobox_value->addItem(QString::fromStdString(conductors[j].m_name));
						if(parameter.m_default_value.AsString() == conductors[j].m_name) {
							combobox_value->setCurrentIndex(j);
						}
					}
					widget_value = combobox_value;
					break;
				}
				case TLINE_PARAMETERTYPE_MATERIAL_DIELECTRIC: {
					QComboBox *combobox_value = new QComboBox(widget);
					for(size_t j = 0; j < dielectrics.size(); ++j) {
						combobox_value->addItem(QString::fromStdString(dielectrics[j].m_name));
						if(parameter.m_default_value.AsString() == dielectrics[j].m_name) {
							combobox_value->setCurrentIndex(j);
						}
					}
					widget_value = combobox_value;
					break;
				}
				case TLINE_PARAMETERTYPE_COUNT: {
					assert(false);
					break;
				}
			}
			QWidget::setTabOrder(widget_focus, widget_value);
			widget_focus = widget_value;
			m_widget_parameters.push_back(widget_value);
			layout->addWidget(label_name, row, 0);
			if(label_unit == NULL) {
				layout->addWidget(widget_value, row, 1, 1, 2);
			} else {
				layout->addWidget(widget_value, row, 1);
				layout->addWidget(label_unit, row, 2);
			}
			++row;
			if(parameter.m_separator) {
				QFrame *line = new QFrame(widget);
				line->setFrameShape(QFrame::HLine);
				line->setFrameShadow(QFrame::Sunken);
				layout->addWidget(line, row, 0, 1, 3);
				++row;
			}
			if(parameter.m_type == TLINE_PARAMETERTYPE_REAL) {
				m_combobox_parameter_sweep_parameter->addItem(QString::fromStdString(parameter.m_name), i);
				m_combobox_parameter_tune_parameter->addItem(QString::fromStdString(parameter.m_name), i);
			}
		}
		layout->setRowStretch(row, 1);
		//layout->setColumnMinimumWidth(0, 200);
		//layout->setColumnMinimumWidth(1, 100);
		//layout->setColumnMinimumWidth(2, 60);
		m_scrollarea_parameters->setWidget(widget);
	}

	// update results
	m_lineedit_results.clear();
	m_combobox_parameter_tune_target_result->clear();
	delete m_scrollarea_results->widget();
	{
		QWidget *widget = new QWidget(m_scrollarea_results), *widget_focus = m_scrollarea_results;
		QGridLayout *layout = new QGridLayout(widget);
		size_t row = 0;
		for(size_t i = 0; i < tline_type.m_modes.size(); ++i) {
			for(size_t j = 0; j < TLINERESULT_COUNT; ++j) {
				QLabel *label_name = new QLabel(QString::fromStdString(tline_type.m_modes[i]) + " " + TLINERESULT_NAMES[j] + ":", widget);
				QLineEdit *lineedit_value = new QLineEdit("?", widget);
				lineedit_value->setReadOnly(true);
				QLabel *label_unit = new QLabel(TLINERESULT_UNITS[j], widget);
				QWidget::setTabOrder(widget_focus, lineedit_value);
				widget_focus = lineedit_value;
				m_lineedit_results.push_back(lineedit_value);
				layout->addWidget(label_name, row, 0);
				layout->addWidget(lineedit_value, row, 1);
				layout->addWidget(label_unit, row, 2);
				++row;
				m_combobox_parameter_tune_target_result->addItem(QString::fromStdString(tline_type.m_modes[i]) + " " + TLINERESULT_NAMES[j]);
			}
			if(i != tline_type.m_modes.size() - 1) {
				QFrame *line = new QFrame(widget);
				line->setFrameShape(QFrame::HLine);
				line->setFrameShadow(QFrame::Sunken);
				layout->addWidget(line, row, 0, 1, 3);
				++row;
			}
		}
		layout->setRowStretch(row, 1);
		//layout->setColumnMinimumWidth(0, 200);
		//layout->setColumnMinimumWidth(1, 100);
		//layout->setColumnMinimumWidth(2, 60);
		m_scrollarea_results->setWidget(widget);
	}

	// update modes
	//index_t current_index = m_combobox_modes->currentIndex();
	m_combobox_modes->clear();
	for(const std::string &mode : tline_type.m_modes) {
		m_combobox_modes->addItem(QString::fromStdString(mode));
	}
	/*if(current_index < mode_count) {
		m_combobox_modes->setCurrentIndex(current_index);
	}*/

	// clear mesh viewer
	m_meshviewer->SetMesh(NULL);

	// do suppressed update
	OnUpdateSimulationType();

}

void MainWindow::OnUpdateSimulationType() {
	SimulationType simulation_type = GetSimulationType();
	MultiGroupVisible({
		{{m_label_frequency[0], m_label_frequency[1], m_lineedit_frequency},
			simulation_type != SIMULATION_FREQUENCY_SWEEP},
		{{m_label_frequency_sweep[0], m_label_frequency_sweep[1], m_label_frequency_sweep[2], m_label_frequency_sweep[3],
		  m_lineedit_frequency_sweep_min, m_lineedit_frequency_sweep_max, m_lineedit_frequency_sweep_step,
		 m_label_frequency_sweep_file, m_lineedit_frequency_sweep_file, m_pushbutton_frequency_sweep_browse},
			simulation_type == SIMULATION_FREQUENCY_SWEEP},
		{{m_label_parameter_sweep[0], m_label_parameter_sweep[1], m_label_parameter_sweep[2], m_label_parameter_sweep[3],
		  m_combobox_parameter_sweep_parameter, m_lineedit_parameter_sweep_min, m_lineedit_parameter_sweep_max, m_lineedit_parameter_sweep_step,
		  m_label_parameter_sweep_file, m_lineedit_parameter_sweep_file, m_pushbutton_parameter_sweep_browse},
			simulation_type == SIMULATION_PARAMETER_SWEEP},
		{{m_label_parameter_tune[0], m_label_parameter_tune[1], m_label_parameter_tune[2],
		  m_combobox_parameter_tune_parameter, m_combobox_parameter_tune_target_result, m_lineedit_parameter_tune_target_value},
			simulation_type == SIMULATION_PARAMETER_TUNE},
	});
}

void MainWindow::OnFrequencySweepBrowse() {
	QString selected_file = QFileDialog::getSaveFileName(this, "Save results as", m_lineedit_frequency_sweep_file->text(), "Text Files (*.txt);;All Files (*)");
	QFileInfo fi(selected_file);
	if(fi.suffix().isEmpty()) {
		selected_file += ".txt";
	}
	m_lineedit_frequency_sweep_file->setText(selected_file);
}

void MainWindow::OnParameterSweepBrowse() {
	QString selected_file = QFileDialog::getSaveFileName(this, "Save results as", m_lineedit_parameter_sweep_file->text(), "Text Files (*.txt);;All Files (*)");
	QFileInfo fi(selected_file);
	if(fi.suffix().isEmpty()) {
		selected_file += ".txt";
	}
	m_lineedit_parameter_sweep_file->setText(selected_file);
}

void MainWindow::OnSimulate() {

	if(m_tline_type == INDEX_NONE)
		return;
	const TLineType &tline_type = g_tline_types[m_tline_type];

	try {

		// initialize context
		TLineContext context;
		context.m_material_database = m_material_database.get();
		GetParameterValues(context.m_parameters);

		// check simulation type
		switch(GetSimulationType()) {
			case SIMULATION_SINGLE_FREQUENCY: {

				// set frequency
				context.m_frequencies = {FloatFromString(m_lineedit_frequency->text().toStdString()) * 1e9};

				// simulate
				tline_type.m_simulate(context);

				break;
			}
			case SIMULATION_FREQUENCY_SWEEP: {

				// generate sweep
				real_t freq_min = FloatFromString(m_lineedit_frequency_sweep_min->text().toStdString()) * 1e9;
				real_t freq_max = FloatFromString(m_lineedit_frequency_sweep_max->text().toStdString()) * 1e9;
				real_t freq_step = FloatFromString(m_lineedit_frequency_sweep_step->text().toStdString()) * 1e9;
				MakeSweep(context.m_frequencies, freq_min, freq_max, freq_step);

				// simulate
				tline_type.m_simulate(context);

				// open output file
				std::string filename = m_lineedit_frequency_sweep_file->text().toLocal8Bit().constData();
				std::ofstream f;
				f.open(filename, std::ios_base::out | std::ios_base::trunc);
				if(f.fail())
					throw std::runtime_error("Could not open file '" + filename + "' for writing.");

				// write header
				f << "Frequency";
				for(size_t i = 0; i < tline_type.m_modes.size(); ++i) {
					for(size_t j = 0; j < TLINERESULT_COUNT; ++j) {
						f << '\t' << tline_type.m_modes[i] << ' ' << TLINERESULT_NAMES[j];
					}
				}
				f << std::endl;

				// write body
				for(index_t i = 0; i < context.m_frequencies.size(); ++i) {
					f << context.m_frequencies[i];
					real_t *results = context.m_results.data() + TLINERESULT_COUNT * tline_type.m_modes.size() * i;
					for(size_t j = 0; j < TLINERESULT_COUNT * tline_type.m_modes.size(); ++j) {
						f << '\t' << results[j];
					}
					f << std::endl;
				}

				break;
			}
			case SIMULATION_PARAMETER_SWEEP: {

				// set frequency
				context.m_frequencies = {FloatFromString(m_lineedit_frequency->text().toStdString()) * 1e9};

				// generate sweep
				real_t value_min = FloatFromString(m_lineedit_parameter_sweep_min->text().toStdString());
				real_t value_max = FloatFromString(m_lineedit_parameter_sweep_max->text().toStdString());
				real_t value_step = FloatFromString(m_lineedit_parameter_sweep_step->text().toStdString());
				std::vector<real_t> sweep_values;
				MakeSweep(sweep_values, value_min, value_max, value_step);

				// simulate
				index_t param_index = m_combobox_parameter_sweep_parameter->itemData(m_combobox_parameter_sweep_parameter->currentIndex()).toUInt();
				std::vector<real_t> combined_results;
				combined_results.resize(TLINERESULT_COUNT * tline_type.m_modes.size() * sweep_values.size());
				for(size_t i = 0; i < sweep_values.size(); ++i) {
					context.m_parameters[param_index].Value() = FloatScale(sweep_values[i]);
					tline_type.m_simulate(context);
					std::copy(context.m_results.begin(), context.m_results.end(), combined_results.begin() + TLINERESULT_COUNT * tline_type.m_modes.size() * i);
				}

				// open output file
				std::string filename = m_lineedit_parameter_sweep_file->text().toLocal8Bit().constData();
				std::ofstream f;
				f.open(filename, std::ios_base::out | std::ios_base::trunc);
				if(f.fail())
					throw std::runtime_error("Could not open file '" + filename + "' for writing.");

				// write header
				f << tline_type.m_parameters[param_index].m_name;
				for(size_t i = 0; i < tline_type.m_modes.size(); ++i) {
					for(size_t j = 0; j < TLINERESULT_COUNT; ++j) {
						f << '\t' << tline_type.m_modes[i] << ' ' << TLINERESULT_NAMES[j];
					}
				}
				f << std::endl;

				// write body
				for(size_t i = 0; i < sweep_values.size(); ++i) {
					f << sweep_values[i];
					real_t *results = combined_results.data() + TLINERESULT_COUNT * tline_type.m_modes.size() * i;
					for(size_t j = 0; j < TLINERESULT_COUNT * tline_type.m_modes.size(); ++j) {
						f << '\t' << results[j];
					}
					f << std::endl;
				}

				break;
			}
			case SIMULATION_PARAMETER_TUNE: {

				// set frequency
				context.m_frequencies = {FloatFromString(m_lineedit_frequency->text().toStdString()) * 1e9};

				// get parameter
				index_t param_index = m_combobox_parameter_tune_parameter->itemData(m_combobox_parameter_tune_parameter->currentIndex()).toUInt();
				const TLineParameter &parameter = tline_type.m_parameters[param_index];
				assert(parameter.m_type == TLINE_PARAMETERTYPE_REAL);

				// get target result
				index_t result_index = clamp<index_t>(m_combobox_parameter_tune_target_result->currentIndex(), 0, TLINERESULT_COUNT * tline_type.m_modes.size() - 1);
				real_t target_value = FloatFromString(m_lineedit_parameter_tune_target_value->text().toStdString());

				// simulate
				real_t initial_value = FloatFromVData(context.m_parameters[param_index].Value());
				if(!FinitePositive(initial_value)) {
					initial_value = FloatFromVData(tline_type.m_parameters[param_index].m_default_value);
				}
				real_t root_value = FindRootRelative([&context, &tline_type, param_index, result_index, target_value](real_t x) {
					context.m_parameters[param_index].Value() = FloatScale(x);
					tline_type.m_simulate(context);
					return context.m_results[result_index] - target_value;
				}, initial_value, 1e-8, target_value * 1e-8, 1e6);

				// write the result back
				QLineEdit *lineedit_value = static_cast<QLineEdit*>(m_widget_parameters[param_index]);
				lineedit_value->setText(QString::number(root_value));

				break;
			}
			case SIMULATION_COUNT: assert(false); break;
		}

		// show last result
		assert(context.m_frequencies.size() > 0);
		assert(context.m_results.size() == TLINERESULT_COUNT * tline_type.m_modes.size() * context.m_frequencies.size());
		real_t *last_results = context.m_results.data() + TLINERESULT_COUNT * tline_type.m_modes.size() * (context.m_frequencies.size() - 1);
		for(size_t i = 0; i < tline_type.m_modes.size(); ++i) {
			for(size_t j = 0; j < TLINERESULT_COUNT; ++j) {
				m_lineedit_results[TLINERESULT_COUNT * i + j]->setText(QString::number(last_results[TLINERESULT_COUNT * i + j]));
			}
		}

		// transfer mesh to viewer
		m_meshviewer->SetMesh(std::move(context.m_output_mesh));

	}
	catch(const std::runtime_error &e) {
		std::cerr << "Error: Could not simulate transmission line:" << std::endl;
		std::cerr << "    " << e.what() << std::endl;
	}

}

void MainWindow::OnZoomChange() {
	m_meshviewer->SetZoom((real_t) m_slider_zoom->value() * 0.000001);
}

void MainWindow::OnImageTypeChange() {
	m_meshviewer->SetImageType((MeshImageType) m_buttongroup_image_type->checkedId());
}

void MainWindow::OnMeshOverlayChange() {
	m_meshviewer->SetMeshOverlay(m_checkbox_mesh_overlay->isChecked());
}

void MainWindow::OnModeChange() {
	m_meshviewer->SetMode(clamp<index_t>(m_combobox_modes->currentIndex(), 0, g_tline_types[m_tline_type].m_modes.size() - 1));
}
