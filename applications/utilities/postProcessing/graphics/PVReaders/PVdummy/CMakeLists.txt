CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

PROJECT(PVdummy)

FIND_PACKAGE(ParaView REQUIRED)

SET(TARGET_NAME PVdummy)

ADD_EXECUTABLE(PVdummy Make/$ENV{WM_OPTIONS}/.cpp)

# Set dependencies
IF(${ParaView_VERSION} VERSION_LESS 5.5)
    TARGET_LINK_LIBRARIES(
        ${TARGET_NAME}
        LINK_PRIVATE
        pqApplicationComponents
    )
ELSEIF(${ParaView_VERSION} VERSION_LESS 5.7)
    # Assume QT version 5 if it hasn't been defined
    IF(NOT PARAVIEW_QT_MAJOR_VERSION)
        set(PARAVIEW_QT_MAJOR_VERSION 5)
    ENDIF()
    # Paraview-5.5.x and higher needs QT libraries listed in the link command
    # in order for the headers to be available
    TARGET_LINK_LIBRARIES(
        ${TARGET_NAME}
        LINK_PUBLIC
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Core
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Gui
        LINK_PRIVATE
        pqApplicationComponents
    )
ELSE()
    # Assume QT version 5 if it hasn't been defined
    IF(NOT PARAVIEW_QT_MAJOR_VERSION)
        set(PARAVIEW_QT_MAJOR_VERSION 5)
    ENDIF()
    # Paraview-5.5.x and higher needs QT libraries listed in the link command
    # in order for the headers to be available
    TARGET_LINK_LIBRARIES(
        ${TARGET_NAME}
        LINK_PUBLIC
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Core
        Qt${PARAVIEW_QT_MAJOR_VERSION}::Gui
        LINK_PRIVATE
        ParaView::pqApplicationComponents
    )
ENDIF()
