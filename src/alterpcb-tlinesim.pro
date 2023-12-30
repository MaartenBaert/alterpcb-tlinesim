QT += core gui

greaterThan(QT_MAJOR_VERSION, 4) {
	QT += widgets
}

TARGET = alterpcb-tlinesim
TEMPLATE = app

DEFINES += "ALTERPCB_VERSION=\\\"0.0.0\\\""

QMAKE_CXXFLAGS += -std=c++11 -Wconversion -Wsign-conversion -Wfloat-conversion
QMAKE_CXXFLAGS_RELEASE -= -O2 -g
QMAKE_CXXFLAGS_RELEASE += -O3 -DNDEBUG

INCLUDEPATH += common gui simulation
DEPENDPATH += common gui simulation

unix {
	CONFIG += link_pkgconfig
	PKGCONFIG += eigen3
}

########## Warning: Everything below this line is auto-generated and will be overwritten! ##########

HEADERS += \
	common/Basics.h \
	common/Color.h \
	common/ColorMap.h \
	common/Cow.h \
	common/Decimal.h \
	common/EnumTranslator.h \
	common/HashTable.h \
	common/Json.h \
	common/MiscMath.h \
	common/MurmurHash.h \
	common/NaturalSort.h \
	common/StringHelper.h \
	common/StringRegistry.h \
	common/VData.h \
	common/VDataReader.h \
	common/Vector.h \
	gui/AboutDialog.h \
	gui/ApplicationDirs.h \
	gui/FixedScrollArea.h \
	gui/GlobalDirs.h \
	gui/Icons.h \
	gui/LayoutHelper.h \
	gui/MainWindow.h \
	gui/MeshViewer.h \
	gui/QLineEditSmall.h \
	gui/QProgressDialogThreaded.h \
	gui/Qt.h \
	simulation/Eigen.h \
	simulation/EigenSparse.h \
	simulation/FindRoot.h \
	simulation/GenericMesh.h \
	simulation/GridMesh2D.h \
	simulation/MaterialDatabase.h \
	simulation/MatrixMarket.h \
	simulation/SparseMatrix.h \
	simulation/TLineTypes.h

SOURCES += \
	Main.cpp \
	common/Color.cpp \
	common/ColorMap.cpp \
	common/Decimal.cpp \
	common/Json.cpp \
	common/NaturalSort.cpp \
	common/StringRegistry.cpp \
	common/VData.cpp \
	gui/AboutDialog.cpp \
	gui/ApplicationDirs.cpp \
	gui/FixedScrollArea.cpp \
	gui/Icons.cpp \
	gui/MainWindow.cpp \
	gui/MeshViewer.cpp \
	gui/QLineEditSmall.cpp \
	gui/QProgressDialogThreaded.cpp \
	simulation/GenericMesh.cpp \
	simulation/GridMesh2D.cpp \
	simulation/MaterialDatabase.cpp \
	simulation/TLineTypes.cpp \
	simulation/TLine_CoplanarWaveguide.cpp \
	simulation/TLine_Microstrip.cpp \
	simulation/TLine_Stripline.cpp
