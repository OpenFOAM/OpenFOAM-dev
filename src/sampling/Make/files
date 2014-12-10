probes/probes.C
probes/patchProbes.C
probes/probesGrouping.C
probes/probesFunctionObject/probesFunctionObject.C

sampledSet/circle/circleSet.C
sampledSet/cloud/cloudSet.C
sampledSet/patchCloud/patchCloudSet.C
sampledSet/polyLine/polyLineSet.C
sampledSet/face/faceOnlySet.C
sampledSet/midPoint/midPointSet.C
sampledSet/midPointAndFace/midPointAndFaceSet.C
sampledSet/patchSeed/patchSeedSet.C
sampledSet/sampledSet/sampledSet.C
sampledSet/sampledSets/sampledSets.C
sampledSet/sampledSets/sampledSetsGrouping.C
sampledSet/sampledSetsFunctionObject/sampledSetsFunctionObject.C
sampledSet/triSurfaceMeshPointSet/triSurfaceMeshPointSet.C
sampledSet/uniform/uniformSet.C
sampledSet/array/arraySet.C

cuttingPlane/cuttingPlane.C

sampledSurface/sampledPatch/sampledPatch.C
sampledSurface/sampledPatchInternalField/sampledPatchInternalField.C
sampledSurface/sampledPlane/sampledPlane.C
sampledSurface/isoSurface/isoSurface.C
sampledSurface/isoSurface/sampledIsoSurface.C
sampledSurface/isoSurface/isoSurfaceCell.C
sampledSurface/isoSurface/sampledIsoSurfaceCell.C
sampledSurface/distanceSurface/distanceSurface.C
sampledSurface/sampledCuttingPlane/sampledCuttingPlane.C
sampledSurface/sampledSurface/sampledSurface.C
sampledSurface/sampledSurfaces/sampledSurfaces.C
sampledSurface/sampledSurfaces/sampledSurfacesGrouping.C
sampledSurface/sampledSurfacesFunctionObject/sampledSurfacesFunctionObject.C
sampledSurface/sampledTriSurfaceMesh/sampledTriSurfaceMesh.C
sampledSurface/thresholdCellFaces/thresholdCellFaces.C
sampledSurface/thresholdCellFaces/sampledThresholdCellFaces.C

surfWriters = sampledSurface/writers

$(surfWriters)/surfaceWriter.C
$(surfWriters)/dx/dxSurfaceWriter.C
$(surfWriters)/ensight/ensightSurfaceWriter.C
$(surfWriters)/ensight/ensightPTraits.C
$(surfWriters)/foamFile/foamFileSurfaceWriter.C
$(surfWriters)/nastran/nastranSurfaceWriter.C
$(surfWriters)/proxy/proxySurfaceWriter.C
$(surfWriters)/raw/rawSurfaceWriter.C
$(surfWriters)/starcd/starcdSurfaceWriter.C
$(surfWriters)/vtk/vtkSurfaceWriter.C

graphField/writePatchGraph.C
graphField/writeCellGraph.C
graphField/makeGraph.C

meshToMesh = meshToMeshInterpolation/meshToMesh
$(meshToMesh)/meshToMesh.C
$(meshToMesh)/meshToMeshParallelOps.C
meshToMeshMethods = meshToMeshInterpolation/meshToMesh/calcMethod
$(meshToMeshMethods)/meshToMeshMethod/meshToMeshMethod.C
$(meshToMeshMethods)/meshToMeshMethod/meshToMeshMethodNew.C
$(meshToMeshMethods)/cellVolumeWeight/cellVolumeWeightMethod.C
$(meshToMeshMethods)/direct/directMethod.C
$(meshToMeshMethods)/mapNearest/mapNearestMethod.C


LIB = $(FOAM_LIBBIN)/libsampling
