# create a plugin that adds a reader to the ParaView GUI
# it is added in the file dialog when doing opens/saves.

# The qrc file is processed by Qt's resource compiler (rcc)
# the qrc file must have a resource prefix of "/ParaViewResources"
# and ParaView will read anything contained under that prefix
# the pqReader.xml file contains xml defining readers with their
# file extensions and descriptions.

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

LINK_DIRECTORIES(
    $ENV{FOAM_LIBBIN}
    $ENV{FOAM_EXT_LIBBIN}
)

INCLUDE_DIRECTORIES(
    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude
    $ENV{WM_PROJECT_DIR}/src/OSspecific/$ENV{WM_OSTYPE}/lnInclude
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude
    ${PROJECT_SOURCE_DIR}/../vtkPVFoam
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
# Defined combined plugin
#

# Extend the auto-generated panel
QT4_WRAP_CPP(MOC_SRCS pqPVFoamReaderPanel.h)

ADD_PARAVIEW_OBJECT_PANEL(IFACES IFACE_SRCS
  CLASS_NAME pqPVFoamReaderPanel
  XML_NAME  PVFoamReader  # name of SourceProxy in *SM.xml
  XML_GROUP sources
)

IF("${PARAVIEW_VERSION_MINOR}" EQUAL 0)
    ADD_PARAVIEW_PLUGIN(
        PVFoamReader_SM "1.0"
        SERVER_MANAGER_XML PVFoamReader_SM.xml
        SERVER_MANAGER_SOURCES vtkPVFoamReader.cxx
        GUI_INTERFACES ${IFACES}
        GUI_SOURCES pqPVFoamReaderPanel.cxx
            ${MOC_SRCS} ${UI_SRCS} ${IFACE_SRCS}
        GUI_RESOURCE_FILES PVFoamReader.xml
    )
ELSE("${PARAVIEW_VERSION_MINOR}" EQUAL 0)
    ADD_PARAVIEW_PLUGIN(
        PVFoamReader_SM "1.0"
        SERVER_MANAGER_XML PVFoamReader_SM.xml
        SERVER_MANAGER_SOURCES vtkPVFoamReader.cxx
        GUI_INTERFACES ${IFACES}
        GUI_SOURCES pqPVFoamReaderPanel.cxx
            ${MOC_SRCS} ${UI_SRCS} ${IFACE_SRCS}
    )
ENDIF("${PARAVIEW_VERSION_MINOR}" EQUAL 0)

TARGET_LINK_LIBRARIES(
    PVFoamReader_SM
    LINK_PUBLIC
    OpenFOAM
    finiteVolume
    vtkPVFoam
)

#-----------------------------------------------------------------------------
