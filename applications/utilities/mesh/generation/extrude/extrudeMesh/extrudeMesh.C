/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    extrudeMesh

Description
    Extrude mesh from existing patch (by default outwards facing normals;
    optional flips faces) or from patch read from file.

    Note: Merges close points so be careful.

    Type of extrusion prescribed by run-time selectable model.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dimensionedTypes.H"
#include "IFstream.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "edgeCollapser.H"
#include "globalMeshData.H"
#include "perfectInterface.H"
#include "addPatchCellLayer.H"
#include "fvMesh.H"
#include "MeshedSurfaces.H"
#include "globalIndex.H"
#include "cellSet.H"

#include "extrudedMesh.H"
#include "extrudeModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

enum ExtrudeMode
{
    MESH,
    PATCH,
    SURFACE
};

namespace Foam
{
    template<>
    const char* NamedEnum<ExtrudeMode, 3>::names[] =
    {
        "mesh",
        "patch",
        "surface"
    };
}

static const NamedEnum<ExtrudeMode, 3> ExtrudeModeNames;


void createDummyFvMeshFiles(const polyMesh& mesh, const word& regionName)
{
    // Create dummy system/fv*
    {
        IOobject io
        (
            "fvSchemes",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        Info<< "Testing:" << io.objectPath() << endl;

        if (!io.headerOk())
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            dictionary divDict;
            dummyDict.add("divSchemes", divDict);
            dictionary gradDict;
            dummyDict.add("gradSchemes", gradDict);
            dictionary laplDict;
            dummyDict.add("laplacianSchemes", laplDict);

            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
    {
        IOobject io
        (
            "fvSolution",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.headerOk())
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
}


label findPatchID(const polyBoundaryMesh& patches, const word& name)
{
    const label patchID = patches.findPatchID(name);

    if (patchID == -1)
    {
        FatalErrorIn("findPatchID(const polyBoundaryMesh&, const word&)")
            << "Cannot find patch " << name
            << " in the source mesh.\n"
            << "Valid patch names are " << patches.names()
            << exit(FatalError);
    }
    return patchID;
}


labelList patchFaces(const polyBoundaryMesh& patches, const wordList& names)
{
    label n = 0;

    forAll(names, i)
    {
        const polyPatch& pp = patches[findPatchID(patches, names[i])];

        n += pp.size();
    }
    labelList faceLabels(n);
    n = 0;
    forAll(names, i)
    {
        const polyPatch& pp = patches[findPatchID(patches, names[i])];

        forAll(pp, j)
        {
            faceLabels[n++] = pp.start()+j;
        }
    }

    return faceLabels;
}


void updateFaceLabels(const mapPolyMesh& map, labelList& faceLabels)
{
    const labelList& reverseMap = map.reverseFaceMap();

    labelList newFaceLabels(faceLabels.size());
    label newI = 0;

    forAll(faceLabels, i)
    {
        label oldFaceI = faceLabels[i];

        if (reverseMap[oldFaceI] >= 0)
        {
            newFaceLabels[newI++] = reverseMap[oldFaceI];
        }
    }
    newFaceLabels.setSize(newI);
    faceLabels.transfer(newFaceLabels);
}


void updateCellSet(const mapPolyMesh& map, labelHashSet& cellLabels)
{
    const labelList& reverseMap = map.reverseCellMap();

    labelHashSet newCellLabels(2*cellLabels.size());

    forAll(cellLabels, i)
    {
        label oldCellI = cellLabels[i];

        if (reverseMap[oldCellI] >= 0)
        {
            newCellLabels.insert(reverseMap[oldCellI]);
        }
    }
    cellLabels.transfer(newCellLabels);
}


int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTimeExtruded.H"

    // Get optional regionName
    word regionName;
    word regionDir;
    if (args.optionReadIfPresent("region", regionName))
    {
        regionDir = regionName;
        Info<< "Create mesh " << regionName << " for time = "
            << runTimeExtruded.timeName() << nl << endl;
    }
    else
    {
        regionName = fvMesh::defaultRegion;
        Info<< "Create mesh for time = "
            << runTimeExtruded.timeName() << nl << endl;
    }


    IOdictionary dict
    (
        IOobject
        (
            "extrudeMeshDict",
            runTimeExtruded.system(),
            runTimeExtruded,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );

    // Point generator
    autoPtr<extrudeModel> model(extrudeModel::New(dict));

    // Whether to flip normals
    const Switch flipNormals(dict.lookup("flipNormals"));

    // What to extrude
    const ExtrudeMode mode = ExtrudeModeNames.read
    (
        dict.lookup("constructFrom")
    );

    // Any merging of small edges
    const scalar mergeTol(dict.lookupOrDefault<scalar>("mergeTol", 1e-4));

    Info<< "Extruding from " << ExtrudeModeNames[mode]
        << " using model " << model().type() << endl;
    if (flipNormals)
    {
        Info<< "Flipping normals before extruding" << endl;
    }
    if (mergeTol > 0)
    {
        Info<< "Collapsing edges < " << mergeTol << " of bounding box" << endl;
    }
    else
    {
        Info<< "Not collapsing any edges after extrusion" << endl;
    }
    Info<< endl;


    // Generated mesh (one of either)
    autoPtr<fvMesh> meshFromMesh;
    autoPtr<polyMesh> meshFromSurface;

    // Faces on front and back for stitching (in case of mergeFaces)
    word frontPatchName;
    labelList frontPatchFaces;
    word backPatchName;
    labelList backPatchFaces;

    // Optional added cells (get written to cellSet)
    labelHashSet addedCellsSet;

    if (mode == PATCH || mode == MESH)
    {
        if (flipNormals && mode == MESH)
        {
            FatalErrorIn(args.executable())
                << "Flipping normals not supported for extrusions from mesh."
                << exit(FatalError);
        }

        fileName sourceCasePath(dict.lookup("sourceCase"));
        sourceCasePath.expand();
        fileName sourceRootDir = sourceCasePath.path();
        fileName sourceCaseDir = sourceCasePath.name();
        if (Pstream::parRun())
        {
            sourceCaseDir =
                sourceCaseDir
               /"processor" + Foam::name(Pstream::myProcNo());
        }
        wordList sourcePatches;
        dict.lookup("sourcePatches") >> sourcePatches;

        if (sourcePatches.size() == 1)
        {
            frontPatchName = sourcePatches[0];
        }

        Info<< "Extruding patches " << sourcePatches
            << " on mesh " << sourceCasePath << nl
            << endl;

        Time runTime
        (
            Time::controlDictName,
            sourceRootDir,
            sourceCaseDir
        );

        #include "createMesh.H"

        const polyBoundaryMesh& patches = mesh.boundaryMesh();


        // Extrusion engine. Either adding to existing mesh or
        // creating separate mesh.
        addPatchCellLayer layerExtrude(mesh, (mode == MESH));

        const labelList meshFaces(patchFaces(patches, sourcePatches));

        if (mode == PATCH && flipNormals)
        {
            // Cheat. Flip patch faces in mesh. This invalidates the
            // mesh (open cells) but does produce the correct extrusion.
            polyTopoChange meshMod(mesh);
            forAll(meshFaces, i)
            {
                label meshFaceI = meshFaces[i];

                label patchI = patches.whichPatch(meshFaceI);
                label own = mesh.faceOwner()[meshFaceI];
                label nei = -1;
                if (patchI == -1)
                {
                    nei = mesh.faceNeighbour()[meshFaceI];
                }

                label zoneI = mesh.faceZones().whichZone(meshFaceI);
                bool zoneFlip = false;
                if (zoneI != -1)
                {
                    label index = mesh.faceZones()[zoneI].whichFace(meshFaceI);
                    zoneFlip = mesh.faceZones()[zoneI].flipMap()[index];
                }

                meshMod.modifyFace
                (
                    mesh.faces()[meshFaceI].reverseFace(),  // modified face
                    meshFaceI,                      // label of face
                    own,                            // owner
                    nei,                            // neighbour
                    true,                           // face flip
                    patchI,                         // patch for face
                    zoneI,                          // zone for face
                    zoneFlip                        // face flip in zone
                );
            }

            // Change the mesh. No inflation.
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

            // Update fields
            mesh.updateMesh(map);

            // Move mesh (since morphing does not do this)
            if (map().hasMotionPoints())
            {
                mesh.movePoints(map().preMotionPoints());
            }
        }



        indirectPrimitivePatch extrudePatch
        (
            IndirectList<face>
            (
                mesh.faces(),
                meshFaces
            ),
            mesh.points()
        );

        // Determine extrudePatch normal
        pointField extrudePatchPointNormals
        (
            PatchTools::pointNormals(mesh, extrudePatch)
        );


        // Precalculate mesh edges for pp.edges.
        const labelList meshEdges
        (
            extrudePatch.meshEdges
            (
                mesh.edges(),
                mesh.pointEdges()
            )
        );

        // Global face indices engine
        const globalIndex globalFaces(mesh.nFaces());

        // Determine extrudePatch.edgeFaces in global numbering (so across
        // coupled patches)
        labelListList edgeGlobalFaces
        (
            addPatchCellLayer::globalEdgeFaces
            (
                mesh,
                globalFaces,
                extrudePatch
            )
        );


        // Determine what patches boundary edges need to get extruded into.
        // This might actually cause edge-connected processors to become
        // face-connected so might need to introduce new processor boundaries.
        // Calculates:
        //  - per pp.edge the patch to extrude into
        //  - any additional processor boundaries (patchToNbrProc = map from
        //    new patchID to neighbour processor)
        //  - number of new patches (nPatches)

        labelList sidePatchID;
        label nPatches;
        Map<label> nbrProcToPatch;
        Map<label> patchToNbrProc;
        addPatchCellLayer::calcSidePatch
        (
            mesh,
            globalFaces,
            edgeGlobalFaces,
            extrudePatch,

            sidePatchID,
            nPatches,
            nbrProcToPatch,
            patchToNbrProc
        );


        // Add any patches.

        label nAdded = nPatches - mesh.boundaryMesh().size();
        reduce(nAdded, sumOp<label>());

        Info<< "Adding overall " << nAdded << " processor patches." << endl;

        if (nAdded > 0)
        {
            DynamicList<polyPatch*> newPatches(nPatches);
            forAll(mesh.boundaryMesh(), patchI)
            {
                newPatches.append
                (
                    mesh.boundaryMesh()[patchI].clone
                    (
                        mesh.boundaryMesh()
                    ).ptr()
                );
            }
            for
            (
                label patchI = mesh.boundaryMesh().size();
                patchI < nPatches;
                patchI++
            )
            {
                label nbrProcI = patchToNbrProc[patchI];

                word name =
                        "procBoundary"
                      + Foam::name(Pstream::myProcNo())
                      + "to"
                      + Foam::name(nbrProcI);

                Pout<< "Adding patch " << patchI
                    << " name:" << name
                    << " between " << Pstream::myProcNo()
                    << " and " << nbrProcI
                    << endl;


                newPatches.append
                (
                    new processorPolyPatch
                    (
                        name,
                        0,                  // size
                        mesh.nFaces(),      // start
                        patchI,             // index
                        mesh.boundaryMesh(),// polyBoundaryMesh
                        Pstream::myProcNo(),// myProcNo
                        nbrProcI            // neighbProcNo
                    )
                );
            }

            // Add patches. Do no parallel updates.
            mesh.removeFvBoundary();
            mesh.addFvPatches(newPatches, true);
        }



        // Only used for addPatchCellLayer into new mesh
        labelList exposedPatchID;
        if (mode == PATCH)
        {
            dict.lookup("exposedPatchName") >> backPatchName;
            exposedPatchID.setSize
            (
                extrudePatch.size(),
                findPatchID(patches, backPatchName)
            );
        }

        // Determine points and extrusion
        pointField layer0Points(extrudePatch.nPoints());
        pointField displacement(extrudePatch.nPoints());
        forAll(displacement, pointI)
        {
            const vector& patchNormal = extrudePatchPointNormals[pointI];

            // layer0 point
            layer0Points[pointI] = model()
            (
                extrudePatch.localPoints()[pointI],
                patchNormal,
                0
            );
            // layerN point
            point extrudePt = model()
            (
                extrudePatch.localPoints()[pointI],
                patchNormal,
                model().nLayers()
            );
            displacement[pointI] = extrudePt - layer0Points[pointI];
        }


        // Check if wedge (has layer0 different from original patch points)
        // If so move the mesh to starting position.
        if (gMax(mag(layer0Points-extrudePatch.localPoints())) > SMALL)
        {
            Info<< "Moving mesh to layer0 points since differ from original"
                << " points - this can happen for wedge extrusions." << nl
                << endl;

            pointField newPoints(mesh.points());
            forAll(extrudePatch.meshPoints(), i)
            {
                newPoints[extrudePatch.meshPoints()[i]] = layer0Points[i];
            }
            mesh.movePoints(newPoints);
        }


        // Layers per face
        labelList nFaceLayers(extrudePatch.size(), model().nLayers());
        // Layers per point
        labelList nPointLayers(extrudePatch.nPoints(), model().nLayers());
        // Displacement for first layer
        vectorField firstLayerDisp(displacement*model().sumThickness(1));

        // Expansion ratio not used.
        scalarField ratio(extrudePatch.nPoints(), 1.0);

        // Topo change container. Either copy an existing mesh or start
        // with empty storage (number of patches only needed for checking)
        autoPtr<polyTopoChange> meshMod
        (
            (
                mode == MESH
              ? new polyTopoChange(mesh)
              : new polyTopoChange(patches.size())
            )
        );

        layerExtrude.setRefinement
        (
            globalFaces,
            edgeGlobalFaces,

            ratio,              // expansion ratio
            extrudePatch,       // patch faces to extrude
            sidePatchID,        // if boundary edge: patch to use
            exposedPatchID,     // if new mesh: patches for exposed faces
            nFaceLayers,
            nPointLayers,
            firstLayerDisp,
            meshMod()
        );

        // Reset points according to extrusion model
        forAll(layerExtrude.addedPoints(), pointI)
        {
            const labelList& pPoints = layerExtrude.addedPoints()[pointI];
            forAll(pPoints, pPointI)
            {
                label meshPointI = pPoints[pPointI];

                point modelPt
                (
                    model()
                    (
                        extrudePatch.localPoints()[pointI],
                        extrudePatchPointNormals[pointI],
                        pPointI+1       // layer
                    )
                );

                const_cast<DynamicList<point>&>
                (
                    meshMod().points()
                )[meshPointI] = modelPt;
            }
        }

        // Store faces on front and exposed patch (if mode=patch there are
        // only added faces so cannot used map to old faces)
        const labelListList& layerFaces = layerExtrude.layerFaces();
        backPatchFaces.setSize(layerFaces.size());
        frontPatchFaces.setSize(layerFaces.size());
        forAll(backPatchFaces, patchFaceI)
        {
            backPatchFaces[patchFaceI]  = layerFaces[patchFaceI].first();
            frontPatchFaces[patchFaceI] = layerFaces[patchFaceI].last();
        }


        // Create dummy fvSchemes, fvSolution
        createDummyFvMeshFiles(mesh, regionDir);

        // Create actual mesh from polyTopoChange container
        autoPtr<mapPolyMesh> map = meshMod().makeMesh
        (
            meshFromMesh,
            IOobject
            (
                regionName,
                runTimeExtruded.constant(),
                runTimeExtruded,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh
        );

        layerExtrude.updateMesh
        (
            map(),
            identity(extrudePatch.size()),
            identity(extrudePatch.nPoints())
        );

        // Calculate face labels for front and back.
        frontPatchFaces = renumber
        (
            map().reverseFaceMap(),
            frontPatchFaces
        );
        backPatchFaces = renumber
        (
            map().reverseFaceMap(),
            backPatchFaces
        );

        // Store added cells
        if (mode == MESH)
        {
            const labelListList addedCells
            (
                layerExtrude.addedCells
                (
                    meshFromMesh,
                    layerExtrude.layerFaces()
                )
            );
            forAll(addedCells, faceI)
            {
                const labelList& aCells = addedCells[faceI];
                forAll(aCells, i)
                {
                    addedCellsSet.insert(aCells[i]);
                }
            }
        }
    }
    else
    {
        // Read from surface
        fileName surfName(dict.lookup("surface"));
        surfName.expand();

        Info<< "Extruding surfaceMesh read from file " << surfName << nl
            << endl;

        MeshedSurface<face> fMesh(surfName);

        if (flipNormals)
        {
            Info<< "Flipping faces." << nl << endl;
            faceList& faces = const_cast<faceList&>(fMesh.faces());
            forAll(faces, i)
            {
                faces[i] = fMesh[i].reverseFace();
            }
        }

        Info<< "Extruding surface with :" << nl
                << "    points     : " << fMesh.points().size() << nl
                << "    faces      : " << fMesh.size() << nl
                << "    normals[0] : " << fMesh.faceNormals()[0]
                << nl
                << endl;

        meshFromSurface.reset
        (
            new extrudedMesh
            (
                IOobject
                (
                    extrudedMesh::defaultRegion,
                    runTimeExtruded.constant(),
                    runTimeExtruded
                ),
                fMesh,
                model()
            )
        );


        // Get the faces on front and back
        frontPatchName = "originalPatch";
        frontPatchFaces = patchFaces
        (
            meshFromSurface().boundaryMesh(),
            wordList(1, frontPatchName)
        );
        backPatchName = "otherSide";
        backPatchFaces = patchFaces
        (
            meshFromSurface().boundaryMesh(),
            wordList(1, backPatchName)
        );
    }


    polyMesh& mesh =
    (
        meshFromMesh.valid()
      ? meshFromMesh()
      : meshFromSurface()
    );


    const boundBox& bb = mesh.bounds();
    const vector span = bb.span();
    const scalar mergeDim = mergeTol * bb.minDim();

    Info<< "Mesh bounding box : " << bb << nl
        << "        with span : " << span << nl
        << "Merge distance    : " << mergeDim << nl
        << endl;


    // Collapse edges
    // ~~~~~~~~~~~~~~

    if (mergeDim > 0)
    {
        Info<< "Collapsing edges < " << mergeDim << " ..." << nl << endl;

        // Edge collapsing engine
        edgeCollapser collapser(mesh);

        const edgeList& edges = mesh.edges();
        const pointField& points = mesh.points();

        PackedBoolList collapseEdge(mesh.nEdges());
        Map<point> collapsePointToLocation(mesh.nPoints());

        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            scalar d = e.mag(points);

            if (d < mergeDim)
            {
                Info<< "Merging edge " << e << " since length " << d
                    << " << " << mergeDim << nl;

                collapseEdge[edgeI] = true;
                collapsePointToLocation.set(e[1], points[e[0]]);
            }
        }

        List<pointEdgeCollapse> allPointInfo;
        const globalIndex globalPoints(mesh.nPoints());
        labelList pointPriority(mesh.nPoints(), 0);

        collapser.consistentCollapse
        (
            globalPoints,
            pointPriority,
            collapsePointToLocation,
            collapseEdge,
            allPointInfo
        );

        // Topo change container
        polyTopoChange meshMod(mesh);

        // Put all modifications into meshMod
        bool anyChange = collapser.setRefinement(allPointInfo, meshMod);

        if (anyChange)
        {
            // Construct new mesh from polyTopoChange.
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

            // Update fields
            mesh.updateMesh(map);

            // Update stored data
            updateFaceLabels(map(), frontPatchFaces);
            updateFaceLabels(map(), backPatchFaces);
            updateCellSet(map(), addedCellsSet);

            // Move mesh (if inflation used)
            if (map().hasMotionPoints())
            {
                mesh.movePoints(map().preMotionPoints());
            }
        }
    }


    // Merging front and back patch faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Switch mergeFaces(dict.lookup("mergeFaces"));
    if (mergeFaces)
    {
        if (mode == MESH)
        {
            FatalErrorIn(args.executable())
                << "Cannot stitch front and back of extrusion since"
                << " in 'mesh' mode (extrusion appended to mesh)."
                << exit(FatalError);
        }

        Info<< "Assuming full 360 degree axisymmetric case;"
            << " stitching faces on patches "
            << frontPatchName << " and "
            << backPatchName << " together ..." << nl << endl;

        if (frontPatchFaces.size() != backPatchFaces.size())
        {
            FatalErrorIn(args.executable())
                << "Differing number of faces on front ("
                << frontPatchFaces.size() << ") and back ("
                << backPatchFaces.size() << ")"
                << exit(FatalError);
        }



        polyTopoChanger stitcher(mesh);
        stitcher.setSize(1);

        const word cutZoneName("originalCutFaceZone");

        List<faceZone*> fz
        (
            1,
            new faceZone
            (
                cutZoneName,
                frontPatchFaces,
                boolList(frontPatchFaces.size(), false),
                0,
                mesh.faceZones()
            )
        );

        mesh.addZones(List<pointZone*>(0), fz, List<cellZone*>(0));

        // Add the perfect interface mesh modifier
        perfectInterface perfectStitcher
        (
            "couple",
            0,
            stitcher,
            cutZoneName,
            word::null,         // dummy patch name
            word::null          // dummy patch name
        );

        // Topo change container
        polyTopoChange meshMod(mesh);

        perfectStitcher.setRefinement
        (
            indirectPrimitivePatch
            (
                IndirectList<face>
                (
                    mesh.faces(),
                    frontPatchFaces
                ),
                mesh.points()
            ),
            indirectPrimitivePatch
            (
                IndirectList<face>
                (
                    mesh.faces(),
                    backPatchFaces
                ),
                mesh.points()
            ),
            meshMod
        );

        // Construct new mesh from polyTopoChange.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(map);

        // Update local data
        updateCellSet(map(), addedCellsSet);

        // Move mesh (if inflation used)
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
    }

    mesh.setInstance(runTimeExtruded.constant());
    Info<< "Writing mesh to " << mesh.objectPath() << nl << endl;

    if (!mesh.write())
    {
        FatalErrorIn(args.executable()) << "Failed writing mesh"
            << exit(FatalError);
    }

    // Need writing cellSet
    label nAdded = returnReduce(addedCellsSet.size(), sumOp<label>());
    if (nAdded > 0)
    {
        cellSet addedCells(mesh, "addedCells", addedCellsSet);
        Info<< "Writing added cells to cellSet " << addedCells.name()
            << nl << endl;
        if (!addedCells.write())
        {
            FatalErrorIn(args.executable()) << "Failed writing cellSet"
                << addedCells.name()
                << exit(FatalError);
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
