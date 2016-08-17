# create a plugin that adds a reader to the ParaView GUI
# it is added in the file dialog when doing opens/saves.

# The qrc file is processed by Qt's resource compiler (rcc)
# the qrc file must have a resource prefix of "/ParaViewResources"
# and ParaView will read anything contained under that prefix
# the pqReader.xml file contains xml defining readers with their
# file extensions and descriptions.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

LINK_DIRECTORIES(
    $ENV{FOAM_LIBBIN}
    $ENV{FOAM_EXT_LIBBIN}
)

INCLUDE_DIRECTORIES(
    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude
    $ENV{WM_PROJECT_DIR}/src/OSspecific/$ENV{WM_OSTYPE}/lnInclude
    $ENV{WM_PROJECT_DIR}/src/meshing/blockMesh/lnInclude
    ${PROJECT_SOURCE_DIR}/../vtkPVblockMesh
)

ADD_DEFINITIONS(
    -std=c++0x
    -DWM_$ENV{WM_PRECISION_OPTION}
    -DWM_LABEL_SIZE=$ENV{WM_LABEL_SIZE}
)

# Set output library destination to plugin directory
SET(
    LIBRARY_OUTPUT_PATH $ENV{PV_PLUGIN_PATH}
    CACHE INTERNAL
    "Single output directory for building all libraries."
)

#
# Define combined plugin
#

# Extend the auto-generated panel
QT4_WRAP_CPP(MOC_SRCS pqPVblockMeshReaderPanel.h)

ADD_PARAVIEW_OBJECT_PANEL(IFACES IFACE_SRCS
  CLASS_NAME pqPVblockMeshReaderPanel
  XML_NAME  PVblockMeshReader # name of SourceProxy in *SM.xml
  XML_GROUP sources
)

IF("${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}" EQUAL 4.0)

    # Special build options, specifically for ParaView 4.0.*
    # Because as of 4.1.0, GUI_RESOURCE_FILES is no longer used.
    ADD_PARAVIEW_PLUGIN(
        PVblockMeshReader_SM "1.0"
        SERVER_MANAGER_XML PVblockMeshReader_SM.xml
        SERVER_MANAGER_SOURCES  vtkPVblockMeshReader.cxx
        GUI_INTERFACES ${IFACES}
        GUI_SOURCES pqPVblockMeshReaderPanel.cxx
        ${MOC_SRCS} ${UI_SRCS} ${IFACE_SRCS}
        GUI_RESOURCE_FILES PVblockMeshReader.xml
    )
ELSE("${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}" EQUAL 4.0)
    ADD_PARAVIEW_PLUGIN(
        PVblockMeshReader_SM "1.0"
        SERVER_MANAGER_XML PVblockMeshReader_SM.xml
        SERVER_MANAGER_SOURCES  vtkPVblockMeshReader.cxx
        GUI_INTERFACES ${IFACES}
        GUI_SOURCES pqPVblockMeshReaderPanel.cxx
        ${MOC_SRCS} ${UI_SRCS} ${IFACE_SRCS}
    )
ENDIF("${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}" EQUAL 4.0)

# Build the client-side plugin

TARGET_LINK_LIBRARIES(
    PVblockMeshReader_SM
    LINK_PUBLIC
    OpenFOAM
    blockMesh
    vtkPVblockMesh
)

#-----------------------------------------------------------------------------
