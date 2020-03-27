# Set up the environment
CMAKE_MINIMUM_REQUIRED(VERSION 3.8)
PROJECT(PVblockMeshReader)

FIND_PACKAGE(ParaView REQUIRED)

# Only ParaView < 5.7 uses include file
IF(${ParaView_VERSION} VERSION_LESS 5.7)
    INCLUDE(${PARAVIEW_USE_FILE})
ENDIF()

LINK_DIRECTORIES(
    $ENV{FOAM_LIBBIN}
    $ENV{FOAM_EXT_LIBBIN}
)

INCLUDE_DIRECTORIES(
    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude
    $ENV{WM_PROJECT_DIR}/src/OSspecific/$ENV{WM_OSTYPE}/lnInclude
    ${PROJECT_SOURCE_DIR}/../vtkPVblockMesh
    ${PROJECT_SOURCE_DIR}/vtk
)

ADD_DEFINITIONS(
    -std=c++11
    -DWM_$ENV{WM_PRECISION_OPTION}
    -DWM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}
)

# Set the output library destination to the plugin directory
SET(
    LIBRARY_OUTPUT_PATH $ENV{PV_PLUGIN_PATH}
    CACHE INTERNAL
    "Single output directory for building all libraries."
)

# Add the plugin. ParaView-4.0.? requires a GUI_RESOURCE_FILES XML file. As of
# version 4.1.? this is no longer needed.
IF(${ParaView_VERSION} VERSION_LESS 4.1)
    ADD_PARAVIEW_PLUGIN(
        PVblockMeshReader_SM
        "1.0"
        SERVER_MANAGER_XML PVblockMeshReader_SM.xml
        SERVER_MANAGER_SOURCES vtk/vtkPVblockMeshReader.cxx
        GUI_RESOURCE_FILES PVblockMeshReader.xml
        )
ELSEIF(${ParaView_VERSION} VERSION_LESS 5.7)
    ADD_PARAVIEW_PLUGIN(
        PVblockMeshReader_SM
        "1.0"
        SERVER_MANAGER_XML PVblockMeshReader_SM.xml
        SERVER_MANAGER_SOURCES vtk/vtkPVblockMeshReader.cxx
    )
ELSE()
    PARAVIEW_ADD_PLUGIN(
        PVblockMeshReader_SM
        VERSION "1.0"
        MODULES PVblockMeshReader_VTK
        MODULE_FILES "vtk/vtk.module"
    )
    SET_TARGET_PROPERTIES(
        PVblockMeshReader_SM
        PROPERTIES
        LIBRARY_OUTPUT_NAME      "libPVblockMeshReader_SM"
        LIBRARY_OUTPUT_DIRECTORY "$ENV{PV_PLUGIN_PATH}"
    )
ENDIF()

# Build the client-side plugin. Paraview-5.5.? needs QT-5 libraries listed in
# the link command in order for the headers to be available.
IF(${ParaView_VERSION} VERSION_LESS 5.5)
    TARGET_LINK_LIBRARIES(
        PVblockMeshReader_SM
        LINK_PUBLIC
        vtkPVblockMesh
        OpenFOAM
    )
ELSE()
    TARGET_LINK_LIBRARIES(
        PVblockMeshReader_SM
        LINK_PUBLIC
        vtkPVblockMesh
        OpenFOAM
        Qt5::Core
        Qt5::Gui
    )
ENDIF()
