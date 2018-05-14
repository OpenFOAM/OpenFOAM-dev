# Set up the environment
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
IF(${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR} VERSION_LESS 4.1)
    ADD_PARAVIEW_PLUGIN(
        PVFoamReader_SM
        "1.0"
        SERVER_MANAGER_XML PVFoamReader_SM.xml
        SERVER_MANAGER_SOURCES vtkPVFoamReader.cxx
        GUI_RESOURCE_FILES PVFoamReader.xml
    )
ELSE()
    ADD_PARAVIEW_PLUGIN(
        PVFoamReader_SM
        "1.0"
        SERVER_MANAGER_XML PVFoamReader_SM.xml
        SERVER_MANAGER_SOURCES vtkPVFoamReader.cxx
    )
ENDIF()

# Build the client-side plugin. Paraview-5.5.? needs QT-5 libraries listed in
# the link command in order for the headers to be available.
IF(${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR} VERSION_LESS 5.5)
    TARGET_LINK_LIBRARIES(
        PVFoamReader_SM
        LINK_PUBLIC
        vtkPVFoam
        finiteVolume
        OpenFOAM
    )
ELSE()
    TARGET_LINK_LIBRARIES(
        PVFoamReader_SM
        LINK_PUBLIC
        vtkPVFoam
        finiteVolume
        OpenFOAM
        Qt5::Core
        Qt5::Gui
    )
ENDIF()
