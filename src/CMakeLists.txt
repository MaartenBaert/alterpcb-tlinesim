set(CMAKE_AUTOMOC ON)

find_package(Qt5 5.7 COMPONENTS Core Gui Widgets REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(alterpcb-tlinesim
	common/Basics.h
	common/Color.cpp
	common/Color.h
	common/ColorMap.cpp
	common/ColorMap.h
	common/Cow.h
	common/Decimal.cpp
	common/Decimal.h
	common/EnumTranslator.h
	common/HashTable.h
	common/Json.cpp
	common/Json.h
	common/MiscMath.h
	common/MurmurHash.h
	common/NaturalSort.cpp
	common/NaturalSort.h
	common/StringHelper.h
	common/StringRegistry.cpp
	common/StringRegistry.h
	common/VData.cpp
	common/VData.h
	common/VDataReader.h
	common/Vector.h
	gui/AboutDialog.cpp
	gui/AboutDialog.h
	gui/ApplicationDirs.cpp
	gui/ApplicationDirs.h
	gui/FixedScrollArea.cpp
	gui/FixedScrollArea.h
	gui/GlobalDirs.h
	gui/Icons.cpp
	gui/Icons.h
	gui/LayoutHelper.h
	gui/MainWindow.cpp
	gui/MainWindow.h
	gui/MeshViewer.cpp
	gui/MeshViewer.h
	gui/QLineEditSmall.cpp
	gui/QLineEditSmall.h
	gui/QProgressDialogThreaded.cpp
	gui/QProgressDialogThreaded.h
	gui/Qt.h
	simulation/Eigen.h
	simulation/EigenSparse.h
	simulation/FemMatrix.cpp
	simulation/FemMatrix.h
	simulation/FindRoot.h
	simulation/GenericMesh.cpp
	simulation/GenericMesh.h
	simulation/GridMesh2D.cpp
	simulation/GridMesh2D.h
	simulation/MaterialDatabase.cpp
	simulation/MaterialDatabase.h
	simulation/Materials.cpp
	simulation/Materials.h
	simulation/MatrixMarket.h
	simulation/SparseMatrix.h
	simulation/TLine_CoplanarWaveguide.cpp
	simulation/TLine_Microstrip.cpp
	simulation/TLine_Stripline.cpp
	simulation/TLineTypes.cpp
	simulation/TLineTypes.h
	Main.cpp
)
target_include_directories(alterpcb-tlinesim PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/common
	${CMAKE_CURRENT_SOURCE_DIR}/gui
	${CMAKE_CURRENT_SOURCE_DIR}/simulation
	/home/maarten/nobackup/build-from-source/eigen
)
target_link_libraries(alterpcb-tlinesim
	Qt5::Core
	Qt5::Gui
	Qt5::Widgets
	# Eigen3::Eigen
)
target_compile_definitions(alterpcb-tlinesim PRIVATE
	-DALTERPCB_VERSION="${PROJECT_VERSION}"
)
