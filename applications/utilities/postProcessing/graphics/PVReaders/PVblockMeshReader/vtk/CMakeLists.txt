VTK_MODULE_ADD_MODULE(PVblockMeshReader_VTK CLASSES vtkPVblockMeshReader)

SET_TARGET_PROPERTIES(
    PVblockMeshReader_VTK
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "$ENV{PV_PLUGIN_PATH}"
)

TARGET_LINK_LIBRARIES(
    PVblockMeshReader_VTK
    LINK_PUBLIC
    vtkPVblockMesh
    OpenFOAM
    Qt5::Core
    Qt5::Gui
)
