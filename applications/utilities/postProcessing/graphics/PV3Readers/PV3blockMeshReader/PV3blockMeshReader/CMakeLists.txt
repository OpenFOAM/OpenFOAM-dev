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
    ${PROJECT_SOURCE_DIR}/../vtkPV3blockMesh
)

ADD_DEFINITIONS(
    -DWM_$ENV{WM_PRECISION_OPTION}
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
QT4_WRAP_CPP(MOC_SRCS pqPV3blockMeshReaderPanel.h)

ADD_PARAVIEW_OBJECT_PANEL(IFACES IFACE_SRCS
  CLASS_NAME pqPV3blockMeshReaderPanel
  XML_NAME  PV3blockMeshReader # name of SourceProxy in *SM.xml
  XML_GROUP sources
)

ADD_PARAVIEW_PLUGIN(
    PV3blockMeshReader_SM "1.0"
    SERVER_MANAGER_XML PV3blockMeshReader_SM.xml
    SERVER_MANAGER_SOURCES  vtkPV3blockMeshReader.cxx
    GUI_INTERFACES ${IFACES}
    GUI_SOURCES pqPV3blockMeshReaderPanel.cxx
        ${MOC_SRCS} ${UI_SRCS} ${IFACE_SRCS}
    GUI_RESOURCE_FILES PV3blockMeshReader.xml
)


# #
# # Define the server-side portion of the reader plugin
# #
# ADD_PARAVIEW_PLUGIN(PV3blockMeshReader_SM "1.0"
#   SERVER_MANAGER_XML PV3blockMeshReader_SM.xml
#   SERVER_MANAGER_SOURCES vtkPV3blockMeshReader.cxx
# )
# #
# # Define the client-side portion of the reader plugin
# #
# ADD_PARAVIEW_PLUGIN(
#     PV3blockMeshReader "1.0"
#     GUI_RESOURCES PV3blockMeshReader.qrc
# )


# Build the client-side plugin

TARGET_LINK_LIBRARIES(
    PV3blockMeshReader_SM
    OpenFOAM
    blockMesh
    vtkPV3blockMesh
)
#-----------------------------------------------------------------------------
