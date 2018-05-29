2013-10-01 Meshing the angledDuct geometry.

constant/triSurface/angledDuct.stl
    outer geometry

constant/triSurface/boundaryAndFaceZones.stl
    boundary and faceZones to extract feature
    edges from.

    (in angledDuctImplicit:
        setSet:
            cellSet porosity new zoneToCell porosity
            cellSet other new cellToCell porosity
            cellSet other invert
            faceSet porosityFaces new cellToFace other all
            faceSet porosityFaces subset cellToFace porosity all
            faceZoneSet porosityFaces new setToFaceZone porosityFaces

        surfaceMeshTriangulate -faceZones '(porosityFaces)' boundaryAndFaceZones.stl
    )

surfaceFeatures



constant/triSurface/porosity_inflated.stl
    block around porosity
    (slightly inflated)

Done paraview:
- start off from porosity.stl
- rotate to align with x axis and translate so symmetric in y:
    Filter->transform
    - translate (0 -0.025 0)
    - rotate (0 0 -45)
- inflate y and z:
    Filter->transform
    - scale (1 1.5 1.5)
- translate back:
    Filter->transform
    - translate (0 0.025 0)
- rotate back:
    Filter->transform
    - rotate (0 0 45)
