Random rndGen(653213);

// Determine mesh bounding boxes:
List<treeBoundBox> meshBb
(
    1,
    treeBoundBox
    (
        boundBox(coarseMesh.points(), false)
    ).extend(rndGen, 1e-3)
);

// Dummy bounds dictionary
dictionary dict;
dict.add("bounds", meshBb);
dict.add
(
    "distributionType",
    distributedTriSurfaceMesh::distributionTypeNames_
    [
        distributedTriSurfaceMesh::FROZEN
    ]
);
dict.add("mergeDistance", SMALL);

labelHashSet includePatches;
forAll(patches, patchI)
{
    const polyPatch& pp = patches[patchI];
    if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp))
    {
        includePatches.insert(patchI);
    }
}

labelList triSurfaceToAgglom(5*nFineFaces);

const triSurface localSurface = triangulate
(
    patches,
    includePatches,
    finalAgglom,
    triSurfaceToAgglom,
    globalNumbering,
    coarsePatches
);


distributedTriSurfaceMesh surfacesMesh
(
    IOobject
    (
        "wallSurface.stl",
        runTime.constant(),     // directory
        "triSurface",           // instance
        runTime,                // registry
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    localSurface,
    dict
);


triSurfaceToAgglom.resize(surfacesMesh.size());

//surfacesMesh.searchableSurface::write();

surfacesMesh.setField(triSurfaceToAgglom);
