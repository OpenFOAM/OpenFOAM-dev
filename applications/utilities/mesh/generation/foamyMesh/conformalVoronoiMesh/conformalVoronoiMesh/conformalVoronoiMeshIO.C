/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "conformalVoronoiMesh.H"
#include "IOstreams.H"
#include "OFstream.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "ListOps.H"
#include "polyMeshFilter.H"
#include "polyTopoChange.H"
#include "PrintTable.H"
#include "pointMesh.H"
#include "indexedVertexOps.H"
#include "DelaunayMeshTools.H"
#include "syncTools.H"
#include "faceSet.H"
#include "OBJstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::timeCheck
(
    const string& description
) const
{
    timeCheck(time(), description, foamyHexMeshControls().timeChecks());
}


void Foam::conformalVoronoiMesh::timeCheck
(
    const Time& runTime,
    const string& description,
    const bool check
)
{
    if (check)
    {
        Info<< nl << "--- [ cpuTime "
            << runTime.elapsedCpuTime() << " s, "
            << "delta " << runTime.cpuTimeIncrement()<< " s";

        if (description != word::null)
        {
            Info<< ", " << description << " ";
        }
        else
        {
            Info<< " ";
        }

        Info<< "] --- " << endl;

        memInfo m;

        if (m.valid())
        {
            PrintTable<word, label> memoryTable
            (
                "Memory Usage (kB): "
              + description
            );

            memoryTable.add("mSize", m.size());
            memoryTable.add("mPeak", m.peak());
            memoryTable.add("mRss", m.rss());

            Info<< incrIndent;
            memoryTable.print(Info, true, true);
            Info<< decrIndent;
        }
    }
}


void Foam::conformalVoronoiMesh::writeMesh(const fileName& instance)
{
    DelaunayMeshTools::writeInternalDelaunayVertices(instance, *this);

    // Per cell the Delaunay vertex
    labelList cellToDelaunayVertex;
    // Per patch, per face the Delaunay vertex
    labelListList patchToDelaunayVertex;

    labelList dualPatchStarts;

    {
        pointField points;
        labelList boundaryPts;
        faceList faces;
        labelList owner;
        labelList neighbour;
        wordList patchNames;
        PtrList<dictionary> patchDicts;
        pointField cellCentres;
        PackedBoolList boundaryFacesToRemove;

        calcDualMesh
        (
            points,
            boundaryPts,
            faces,
            owner,
            neighbour,
            patchNames,
            patchDicts,
            cellCentres,
            cellToDelaunayVertex,
            patchToDelaunayVertex,
            boundaryFacesToRemove
        );

        Info<< nl << "Writing polyMesh to " << instance << endl;

        writeMesh
        (
            Foam::polyMesh::defaultRegion,
            instance,
            points,
            boundaryPts,
            faces,
            owner,
            neighbour,
            patchNames,
            patchDicts,
            cellCentres,
            boundaryFacesToRemove
        );

        dualPatchStarts.setSize(patchDicts.size());

        forAll(dualPatchStarts, patchi)
        {
            dualPatchStarts[patchi] =
                readLabel(patchDicts[patchi].lookup("startFace"));
        }
    }

    if (foamyHexMeshControls().writeCellShapeControlMesh())
    {
        cellShapeControls().shapeControlMesh().write();
    }

    if (foamyHexMeshControls().writeBackgroundMeshDecomposition())
    {
        Info<< nl << "Writing " << "backgroundMeshDecomposition" << endl;

        // Have to explicitly update the mesh instance.
        const_cast<fvMesh&>(decomposition_().mesh()).setInstance
        (
            time().timeName()
        );

        decomposition_().mesh().write();
    }

    if (foamyHexMeshControls().writeTetDualMesh())
    {
        label celli = 0;
        for
        (
            Finite_cells_iterator cit = finite_cells_begin();
            cit != finite_cells_end();
            ++cit
        )
        {
            if
            (
                !cit->hasFarPoint()
             && !is_infinite(cit)
            )
            {
                cit->cellIndex() = celli++;
            }
        }

        Info<< nl << "Writing " << "tetDualMesh" << endl;

        DistributedDelaunayMesh<Delaunay>::labelTolabelPairHashTable vertexMap;
        labelList cellMap;
        autoPtr<polyMesh> tetMesh =
            createMesh("tetDualMesh", vertexMap, cellMap);

        tetMesh().write();

//        // Determine map from Delaunay vertex to Dual mesh
//        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//        // From all Delaunay vertices to cell (positive index)
//        // or patch face (negative index)
//        labelList vertexToDualAddressing(number_of_vertices(), 0);
//
//        forAll(cellToDelaunayVertex, celli)
//        {
//            label vertI = cellToDelaunayVertex[celli];
//
//            if (vertexToDualAddressing[vertI] != 0)
//            {
//                FatalErrorInFunction
//                    << "Delaunay vertex " << vertI
//                    << " from cell " << celli
//                    << " is already mapped to "
//                    << vertexToDualAddressing[vertI]
//                    << exit(FatalError);
//            }
//            vertexToDualAddressing[vertI] = celli+1;
//        }
//
//        forAll(patchToDelaunayVertex, patchi)
//        {
//            const labelList& patchVertices = patchToDelaunayVertex[patchi];
//
//            forAll(patchVertices, i)
//            {
//                label vertI = patchVertices[i];
//
//                if (vertexToDualAddressing[vertI] > 0)
//                {
//                    FatalErrorInFunction
//                        << "Delaunay vertex " << vertI
//                        << " from patch " << patchi
//                        << " local index " << i
//                        << " is already mapped to cell "
//                        << vertexToDualAddressing[vertI]-1
//                        << exit(FatalError);
//                }
//
//                // Vertex might be used by multiple faces. Which one to
//                // use? For now last one wins.
//                label dualFacei = dualPatchStarts[patchi]+i;
//                vertexToDualAddressing[vertI] = -dualFacei-1;
//            }
//        }
//
//
//        // Calculate tet mesh addressing
//        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//        pointField points;
//        labelList boundaryPts(number_of_finite_cells(), -1);
//        // From tet point back to Delaunay vertex index
//        labelList pointToDelaunayVertex;
//        faceList faces;
//        labelList owner;
//        labelList neighbour;
//        wordList patchTypes;
//        wordList patchNames;
//        PtrList<dictionary> patchDicts;
//        pointField cellCentres;
//
//        calcTetMesh
//        (
//            points,
//            pointToDelaunayVertex,
//            faces,
//            owner,
//            neighbour,
//            patchTypes,
//            patchNames,
//            patchDicts
//        );
//
//
//
//        // Calculate map from tet points to dual mesh cells/patch faces
//        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//        labelIOList pointDualAddressing
//        (
//            IOobject
//            (
//                "pointDualAddressing",
//                instance,
//                "tetDualMesh"/polyMesh::meshSubDir,
//                runTime_,
//                IOobject::NO_READ,
//                IOobject::AUTO_WRITE,
//                false
//            ),
//            UIndirectList<label>
//            (
//                vertexToDualAddressing,
//                pointToDelaunayVertex
//            )()
//        );
//
//        label pointi = findIndex(pointDualAddressing, -1);
//        if (pointi != -1)
//        {
//            WarningInFunction
//                << "Delaunay vertex " << pointi
//                << " does not have a corresponding dual cell." << endl;
//        }
//
//        Info<< "Writing map from tetDualMesh points to Voronoi mesh to "
//            << pointDualAddressing.objectPath() << endl;
//        pointDualAddressing.write();
//
//
//
//        // Write tet points corresponding to the Voronoi cell/face centre
//        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//        {
//            // Read Voronoi mesh
//            fvMesh mesh
//            (
//                IOobject
//                (
//                    Foam::polyMesh::defaultRegion,
//                    instance,
//                    runTime_,
//                    IOobject::MUST_READ
//                )
//            );
//            pointIOField dualPoints
//            (
//                IOobject
//                (
//                    "dualPoints",
//                    instance,
//                    "tetDualMesh"/polyMesh::meshSubDir,
//                    runTime_,
//                    IOobject::NO_READ,
//                    IOobject::AUTO_WRITE,
//                    false
//                ),
//                points
//            );
//
//            forAll(pointDualAddressing, pointi)
//            {
//                label index = pointDualAddressing[pointi];
//
//                if (index > 0)
//                {
//                    label celli = index-1;
//                    dualPoints[pointi] = mesh.cellCentres()[celli];
//                }
//                else if (index < 0)
//                {
//                    label facei = -index-1;
//                    if (facei >= mesh.nInternalFaces())
//                    {
//                        dualPoints[pointi] = mesh.faceCentres()[facei];
//                    }
//                }
//            }
//
//            Info<< "Writing tetDualMesh points mapped onto Voronoi mesh to "
//                << dualPoints.objectPath() << endl
//                << "Replace the polyMesh/points with these." << endl;
//            dualPoints.write();
//        }
//
//
//        Info<< nl << "Writing tetDualMesh to " << instance << endl;
//
//        PackedBoolList boundaryFacesToRemove;
//        writeMesh
//        (
//            "tetDualMesh",
//            instance,
//            points,
//            boundaryPts,
//            faces,
//            owner,
//            neighbour,
//            patchTypes,
//            patchNames,
//            patchDicts,
//            cellCentres,
//            boundaryFacesToRemove
//        );
    }
}


Foam::autoPtr<Foam::fvMesh> Foam::conformalVoronoiMesh::createDummyMesh
(
    const IOobject& io,
    const wordList& patchNames,
    const PtrList<dictionary>& patchDicts
) const
{
    autoPtr<fvMesh> meshPtr
    (
        new fvMesh
        (
            io,
            xferCopy(pointField()),
            xferCopy(faceList()),
            xferCopy(cellList())
        )
    );
    fvMesh& mesh = meshPtr();

    List<polyPatch*> patches(patchDicts.size());

    forAll(patches, patchi)
    {
        if
        (
            patchDicts.set(patchi)
         && (
                word(patchDicts[patchi].lookup("type"))
             == processorPolyPatch::typeName
            )
        )
        {
            patches[patchi] = new processorPolyPatch
            (
                0,          // patchSizes[p],
                0,          // patchStarts[p],
                patchi,
                mesh.boundaryMesh(),
                readLabel(patchDicts[patchi].lookup("myProcNo")),
                readLabel(patchDicts[patchi].lookup("neighbProcNo")),
                coupledPolyPatch::COINCIDENTFULLMATCH
            );
        }
        else
        {
            patches[patchi] = polyPatch::New
            (
                patchDicts[patchi].lookup("type"),
                patchNames[patchi],
                0,          // patchSizes[p],
                0,          // patchStarts[p],
                patchi,
                mesh.boundaryMesh()
            ).ptr();
        }
    }

    mesh.addFvPatches(patches);

    return meshPtr;
}


void Foam::conformalVoronoiMesh::checkProcessorPatchesMatch
(
    const PtrList<dictionary>& patchDicts
) const
{
    // Check patch sizes
    labelListList procPatchSizes
    (
        Pstream::nProcs(),
        labelList(Pstream::nProcs(), -1)
    );

    forAll(patchDicts, patchi)
    {
        if
        (
            patchDicts.set(patchi)
         && (
                word(patchDicts[patchi].lookup("type"))
             == processorPolyPatch::typeName
            )
        )
        {
            const label procNeighb =
                readLabel(patchDicts[patchi].lookup("neighbProcNo"));

            procPatchSizes[Pstream::myProcNo()][procNeighb]
                = readLabel(patchDicts[patchi].lookup("nFaces"));
        }
    }

    Pstream::gatherList(procPatchSizes);

    if (Pstream::master())
    {
        bool allMatch = true;

        forAll(procPatchSizes, proci)
        {
            const labelList& patchSizes = procPatchSizes[proci];

            forAll(patchSizes, patchi)
            {
                if (patchSizes[patchi] != procPatchSizes[patchi][proci])
                {
                    allMatch = false;

                    Info<< indent << "Patches " << proci << " and " << patchi
                        << " have different sizes: " << patchSizes[patchi]
                        << " and " << procPatchSizes[patchi][proci] << endl;
                }
            }
        }

        if (allMatch)
        {
            Info<< indent << "All processor patches have matching numbers of "
                << "faces" << endl;
        }
    }
}


void Foam::conformalVoronoiMesh::reorderPoints
(
    pointField& points,
    labelList& boundaryPts,
    faceList& faces,
    const label nInternalFaces
) const
{
    Info<< incrIndent << indent << "Reordering points into internal/external"
        << endl;

    labelList oldToNew(points.size(), label(0));

    // Find points that are internal
    for (label fI = nInternalFaces; fI < faces.size(); ++fI)
    {
        const face& f = faces[fI];

        forAll(f, fpI)
        {
            oldToNew[f[fpI]] = 1;
        }
    }

    const label nInternalPoints = points.size() - sum(oldToNew);

    label countInternal = 0;
    label countExternal = nInternalPoints;

    forAll(points, pI)
    {
        if (oldToNew[pI] == 0)
        {
            oldToNew[pI] = countInternal++;
        }
        else
        {
            oldToNew[pI] = countExternal++;
        }
    }

    Info<< indent
        << "Number of internal points: " << countInternal << nl
        << indent << "Number of external points: " << countExternal
        << decrIndent << endl;

    inplaceReorder(oldToNew, points);
    inplaceReorder(oldToNew, boundaryPts);

    forAll(faces, fI)
    {
        face& f = faces[fI];

        forAll(f, fpI)
        {
            f[fpI] = oldToNew[f[fpI]];
        }
    }
}


void Foam::conformalVoronoiMesh::reorderProcessorPatches
(
    const word& meshName,
    const fileName& instance,
    const pointField& points,
    faceList& faces,
    const wordList& patchNames,
    const PtrList<dictionary>& patchDicts
) const
{
    Info<< incrIndent << indent << "Reordering processor patches" << endl;

    Info<< incrIndent;
    checkProcessorPatchesMatch(patchDicts);

    // Create dummy mesh with correct proc boundaries to do sorting
    autoPtr<fvMesh> sortMeshPtr
    (
        createDummyMesh
        (
            IOobject
            (
                meshName,
                instance,
                runTime_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            patchNames,
            patchDicts
        )
    );
    const fvMesh& sortMesh = sortMeshPtr();

    // Rotation on new faces.
    labelList rotation(faces.size(), label(0));
    labelList faceMap(faces.size(), label(-1));

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    // Send ordering
    forAll(sortMesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = sortMesh.boundaryMesh()[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            refCast<const processorPolyPatch>(pp).initOrder
            (
                pBufs,
                primitivePatch
                (
                    SubList<face>
                    (
                        faces,
                        readLabel(patchDicts[patchi].lookup("nFaces")),
                        readLabel(patchDicts[patchi].lookup("startFace"))
                    ),
                    points
                )
            );
        }
    }

    pBufs.finishedSends();

    Info<< incrIndent << indent << "Face ordering initialised..." << endl;

    // Receive and calculate ordering
    bool anyChanged = false;

    forAll(sortMesh.boundaryMesh(), patchi)
    {
        const polyPatch& pp = sortMesh.boundaryMesh()[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            const label nPatchFaces =
                readLabel(patchDicts[patchi].lookup("nFaces"));
            const label patchStartFace =
                readLabel(patchDicts[patchi].lookup("startFace"));

            labelList patchFaceMap(nPatchFaces, label(-1));
            labelList patchFaceRotation(nPatchFaces, label(0));

            bool changed = refCast<const processorPolyPatch>(pp).order
            (
                pBufs,
                primitivePatch
                (
                    SubList<face>
                    (
                        faces,
                        nPatchFaces,
                        patchStartFace
                    ),
                    points
                ),
                patchFaceMap,
                patchFaceRotation
            );

            if (changed)
            {
                // Merge patch face reordering into mesh face reordering table
                forAll(patchFaceRotation, patchFacei)
                {
                    rotation[patchFacei + patchStartFace]
                        = patchFaceRotation[patchFacei];
                }

                forAll(patchFaceMap, patchFacei)
                {
                    if (patchFaceMap[patchFacei] != patchFacei)
                    {
                        faceMap[patchFacei + patchStartFace]
                            = patchFaceMap[patchFacei] + patchStartFace;
                    }
                }

                anyChanged = true;
            }
        }
    }

    Info<< incrIndent << indent << "Faces matched." << endl;

    reduce(anyChanged, orOp<bool>());

    if (anyChanged)
    {
        label nReorderedFaces = 0;

        forAll(faceMap, facei)
        {
           if (faceMap[facei] != -1)
           {
               nReorderedFaces++;
           }
        }

        if (nReorderedFaces > 0)
        {
            inplaceReorder(faceMap, faces);
        }

        // Rotate faces (rotation is already in new face indices).
        label nRotated = 0;

        forAll(rotation, facei)
        {
            if (rotation[facei] != 0)
            {
                faces[facei] = rotateList(faces[facei], rotation[facei]);
                nRotated++;
            }
        }

        Info<< indent << returnReduce(nReorderedFaces, sumOp<label>())
            << " faces have been reordered" << nl
            << indent << returnReduce(nRotated, sumOp<label>())
            << " faces have been rotated"
            << decrIndent << decrIndent
            << decrIndent << decrIndent << endl;
    }
}


void Foam::conformalVoronoiMesh::writeMesh
(
    const word& meshName,
    const fileName& instance,
    pointField& points,
    labelList& boundaryPts,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    const wordList& patchNames,
    const PtrList<dictionary>& patchDicts,
    const pointField& cellCentres,
    PackedBoolList& boundaryFacesToRemove
) const
{
    if (foamyHexMeshControls().objOutput())
    {
        DelaunayMeshTools::writeObjMesh
        (
            time().path()/word(meshName + ".obj"),
            points,
            faces
        );
    }

    const label nInternalFaces = readLabel(patchDicts[0].lookup("startFace"));

    reorderPoints(points, boundaryPts, faces, nInternalFaces);

    if (Pstream::parRun())
    {
        reorderProcessorPatches
        (
            meshName,
            instance,
            points,
            faces,
            patchNames,
            patchDicts
        );
    }

    Info<< incrIndent;
    Info<< indent << "Constructing mesh" << endl;

    timeCheck("Before fvMesh construction");

    fvMesh mesh
    (
        IOobject
        (
            meshName,
            instance,
            runTime_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        xferMove(points),
        xferMove(faces),
        xferMove(owner),
        xferMove(neighbour)
    );

    Info<< indent << "Adding patches to mesh" << endl;

    List<polyPatch*> patches(patchNames.size());

    label nValidPatches = 0;

    forAll(patches, p)
    {
        label totalPatchSize = readLabel(patchDicts[p].lookup("nFaces"));

        if
        (
            patchDicts.set(p)
         && (
                word(patchDicts[p].lookup("type"))
             == processorPolyPatch::typeName
            )
        )
        {
            const_cast<dictionary&>(patchDicts[p]).set
            (
                "transform",
                "coincidentFullMatch"
            );

            // Do not create empty processor patches
            if (totalPatchSize > 0)
            {
                patches[nValidPatches] = new processorPolyPatch
                (
                    patchNames[p],
                    patchDicts[p],
                    nValidPatches,
                    mesh.boundaryMesh(),
                    processorPolyPatch::typeName
                );

                nValidPatches++;
            }
        }
        else
        {
            // Check that the patch is not empty on every processor
            reduce(totalPatchSize, sumOp<label>());

            if (totalPatchSize > 0)
            {
                patches[nValidPatches] = polyPatch::New
                (
                    patchNames[p],
                    patchDicts[p],
                    nValidPatches,
                    mesh.boundaryMesh()
                ).ptr();

                nValidPatches++;
            }
        }
    }

    patches.setSize(nValidPatches);

    mesh.addFvPatches(patches);



    // Add zones to the mesh
    addZones(mesh, cellCentres);



    Info<< indent << "Add pointZones" << endl;
    {
        label sz = mesh.pointZones().size();

        DynamicList<label> bPts(boundaryPts.size());

        forAll(dualMeshPointTypeNames_, typeI)
        {
            forAll(boundaryPts, ptI)
            {
                const label& bPtType = boundaryPts[ptI];

                if (bPtType == typeI)
                {
                    bPts.append(ptI);
                }
            }

//            syncTools::syncPointList(mesh, bPts, maxEqOp<label>(), -1);

            Info<< incrIndent << indent
                << "Adding " << bPts.size()
                << " points of type " << dualMeshPointTypeNames_.words()[typeI]
                << decrIndent << endl;

            mesh.pointZones().append
            (
                new pointZone
                (
                    dualMeshPointTypeNames_.words()[typeI],
                    bPts,
                    sz + typeI,
                    mesh.pointZones()
                )
            );

            bPts.clear();
        }
    }



    // Add indirectPatchFaces to a face zone
    Info<< indent << "Adding indirect patch faces set" << endl;

    syncTools::syncFaceList
    (
        mesh,
        boundaryFacesToRemove,
        orEqOp<unsigned int>()
    );

    labelList addr(boundaryFacesToRemove.count());
    label count = 0;

    forAll(boundaryFacesToRemove, facei)
    {
        if (boundaryFacesToRemove[facei])
        {
            addr[count++] = facei;
        }
    }

    addr.setSize(count);

    faceSet indirectPatchFaces
    (
        mesh,
        "indirectPatchFaces",
        addr,
        IOobject::AUTO_WRITE
    );

    indirectPatchFaces.sync(mesh);


    Info<< decrIndent;

    timeCheck("Before fvMesh filtering");

    autoPtr<polyMeshFilter> meshFilter;

    label nInitialBadFaces = 0;

    if (foamyHexMeshControls().filterEdges())
    {
        Info<< nl << "Filtering edges on polyMesh" << nl << endl;

        meshFilter.reset(new polyMeshFilter(mesh, boundaryPts));

        // Filter small edges only. This reduces the number of faces so that
        // the face filtering is sped up.
        nInitialBadFaces = meshFilter().filterEdges(0);
        {
            const autoPtr<fvMesh>& newMesh = meshFilter().filteredMesh();

            polyTopoChange meshMod(newMesh());

            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

            polyMeshFilter::copySets(newMesh(), mesh);
        }
    }

    if (foamyHexMeshControls().filterFaces())
    {
        labelIOList boundaryPtsIO
        (
            IOobject
            (
                "pointPriority",
                instance,
                time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            labelList(mesh.nPoints(), labelMin)
        );

        forAll(mesh.points(), ptI)
        {
            boundaryPtsIO[ptI] = mesh.pointZones().whichZone(ptI);
        }


        Info<< nl << "Filtering faces on polyMesh" << nl << endl;

        meshFilter.reset(new polyMeshFilter(mesh, boundaryPtsIO));

        meshFilter().filter(nInitialBadFaces);
        {
            const autoPtr<fvMesh>& newMesh = meshFilter().filteredMesh();

            polyTopoChange meshMod(newMesh());

            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

            polyMeshFilter::copySets(newMesh(), mesh);
        }
    }

    timeCheck("After fvMesh filtering");

    mesh.setInstance(instance);

    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }
    else
    {
        Info<< nl << "Written filtered mesh to "
            << mesh.polyMesh::instance() << nl
            << endl;
    }

    {
        pointScalarField boundaryPtsScalarField
        (
            IOobject
            (
                "boundaryPoints_collapsed",
                instance,
                time(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pointMesh::New(mesh),
            scalar(labelMin)
        );

        labelIOList boundaryPtsIO
        (
            IOobject
            (
                "pointPriority",
                instance,
                time(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            labelList(mesh.nPoints(), labelMin)
        );

        forAll(mesh.points(), ptI)
        {
            boundaryPtsScalarField[ptI] = mesh.pointZones().whichZone(ptI);
            boundaryPtsIO[ptI] = mesh.pointZones().whichZone(ptI);
        }

        boundaryPtsScalarField.write();
        boundaryPtsIO.write();
    }

//    writeCellSizes(mesh);

//    writeCellAlignments(mesh);

//    writeCellCentres(mesh);

    findRemainingProtrusionSet(mesh);
}


void Foam::conformalVoronoiMesh::writeCellSizes
(
    const fvMesh& mesh
) const
{
    {
        timeCheck("Start writeCellSizes");

        Info<< nl << "Create targetCellSize volScalarField" << endl;

        volScalarField targetCellSize
        (
            IOobject
            (
                "targetCellSize",
                mesh.polyMesh::instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("cellSize", dimLength, 0),
            zeroGradientFvPatchScalarField::typeName
        );

        scalarField& cellSize = targetCellSize.primitiveFieldRef();

        const vectorField& C = mesh.cellCentres();

        forAll(cellSize, i)
        {
            cellSize[i] = cellShapeControls().cellSize(C[i]);
        }

        // Info<< nl << "Create targetCellVolume volScalarField" << endl;

        // volScalarField targetCellVolume
        // (
        //     IOobject
        //     (
        //         "targetCellVolume",
        //         mesh.polyMesh::instance(),
        //         mesh,
        //         IOobject::NO_READ,
        //         IOobject::AUTO_WRITE
        //     ),
        //     mesh,
        //     dimensionedScalar("cellVolume", dimLength, 0),
        //     zeroGradientFvPatchScalarField::typeName
        // );

        // targetCellVolume.primitiveFieldRef() = pow3(cellSize);

        // Info<< nl << "Create actualCellVolume volScalarField" << endl;

        // volScalarField actualCellVolume
        // (
        //     IOobject
        //     (
        //         "actualCellVolume",
        //         mesh.polyMesh::instance(),
        //         mesh,
        //         IOobject::NO_READ,
        //         IOobject::AUTO_WRITE
        //     ),
        //     mesh,
        //     dimensionedScalar("cellVolume", dimVolume, 0),
        //     zeroGradientFvPatchScalarField::typeName
        // );

        // actualCellVolume.primitiveFieldRef() = mesh.cellVolumes();

        // Info<< nl << "Create equivalentCellSize volScalarField" << endl;

        // volScalarField equivalentCellSize
        // (
        //     IOobject
        //     (
        //         "equivalentCellSize",
        //         mesh.polyMesh::instance(),
        //         mesh,
        //         IOobject::NO_READ,
        //         IOobject::AUTO_WRITE
        //     ),
        //     mesh,
        //     dimensionedScalar("cellSize", dimLength, 0),
        //     zeroGradientFvPatchScalarField::typeName
        // );

        // equivalentCellSize.primitiveFieldRef() = pow
        // (
        //     actualCellVolume.primitiveField(),
        //     1.0/3.0
        // );

        targetCellSize.correctBoundaryConditions();
        // targetCellVolume.correctBoundaryConditions();
        // actualCellVolume.correctBoundaryConditions();
        // equivalentCellSize.correctBoundaryConditions();

        targetCellSize.write();
        // targetCellVolume.write();
        // actualCellVolume.write();
        // equivalentCellSize.write();
    }

    // {
    //     polyMesh tetMesh
    //     (
    //         IOobject
    //         (
    //             "tetDualMesh",
    //             runTime_.constant(),
    //             runTime_,
    //             IOobject::MUST_READ
    //         )
    //     );

    //     pointMesh ptMesh(tetMesh);

    //     pointScalarField ptTargetCellSize
    //     (
    //         IOobject
    //         (
    //             "ptTargetCellSize",
    //             runTime_.timeName(),
    //             tetMesh,
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         ptMesh,
    //         dimensionedScalar("ptTargetCellSize", dimLength, 0),
    //         pointPatchVectorField::calculatedType()
    //     );

    //     scalarField& cellSize = ptTargetCellSize.primitiveFieldRef();

    //     const vectorField& P = tetMesh.points();

    //     forAll(cellSize, i)
    //     {
    //         cellSize[i] = cellShapeControls().cellSize(P[i]);
    //     }

    //     ptTargetCellSize.write();
    // }
}


void Foam::conformalVoronoiMesh::writeCellAlignments
(
    const fvMesh& mesh
) const
{
//    Info<< nl << "Create cellAlignments volTensorField" << endl;
//
//    volTensorField cellAlignments
//    (
//        IOobject
//        (
//            "cellAlignments",
//            mesh.polyMesh::instance(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        tensor::I,
//        zeroGradientFvPatchTensorField::typeName
//    );
//
//    tensorField& cellAlignment = cellAlignments.primitiveFieldRef();
//
//    const vectorField& C = mesh.cellCentres();
//
//    vectorField xDir(cellAlignment.size());
//    vectorField yDir(cellAlignment.size());
//    vectorField zDir(cellAlignment.size());
//
//    forAll(cellAlignment, i)
//    {
//        cellAlignment[i] = cellShapeControls().cellAlignment(C[i]);
//        xDir[i] = cellAlignment[i] & vector(1, 0, 0);
//        yDir[i] = cellAlignment[i] & vector(0, 1, 0);
//        zDir[i] = cellAlignment[i] & vector(0, 0, 1);
//    }
//
//    OFstream xStr("xDir.obj");
//    OFstream yStr("yDir.obj");
//    OFstream zStr("zDir.obj");
//
//    forAll(xDir, i)
//    {
//        meshTools::writeOBJ(xStr, C[i], C[i] + xDir[i]);
//        meshTools::writeOBJ(yStr, C[i], C[i] + yDir[i]);
//        meshTools::writeOBJ(zStr, C[i], C[i] + zDir[i]);
//    }
//
//    cellAlignments.correctBoundaryConditions();
//
//    cellAlignments.write();
}


void Foam::conformalVoronoiMesh::writeCellCentres
(
    const fvMesh& mesh
) const
{
    Info<< "Writing components of cellCentre positions to volScalarFields"
        << " ccx, ccy, ccz in " <<  runTime_.timeName() << endl;

    for (direction i=0; i<vector::nComponents; i++)
    {
        volScalarField cci
        (
            IOobject
            (
                "cc" + word(vector::componentNames[i]),
                runTime_.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh.C().component(i)
        );

        cci.write();
    }
}


Foam::labelHashSet Foam::conformalVoronoiMesh::findRemainingProtrusionSet
(
    const polyMesh& mesh
) const
{
    timeCheck("Start findRemainingProtrusionSet");

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelHashSet protrudingBoundaryPoints;

    forAll(patches, patchi)
    {
        const polyPatch& patch = patches[patchi];

        forAll(patch.localPoints(), pLPI)
        {
            label meshPtI = patch.meshPoints()[pLPI];

            const Foam::point& pt = patch.localPoints()[pLPI];

            if
            (
                geometryToConformTo_.wellOutside
                (
                    pt,
                    sqr(targetCellSize(pt))
                )
            )
            {
                protrudingBoundaryPoints.insert(meshPtI);
            }
        }
    }

    cellSet protrudingCells
    (
        mesh,
        "foamyHexMesh_remainingProtrusions",
        mesh.nCells()/1000
    );

    forAllConstIter(labelHashSet, protrudingBoundaryPoints, iter)
    {
        const label pointi = iter.key();
        const labelList& pCells = mesh.pointCells()[pointi];

        forAll(pCells, pCI)
        {
            protrudingCells.insert(pCells[pCI]);
        }
    }

    label protrudingCellsSize = protrudingCells.size();

    reduce(protrudingCellsSize, sumOp<label>());

    if (foamyHexMeshControls().objOutput() && protrudingCellsSize > 0)
    {
        Info<< nl << "Found " << protrudingCellsSize
            << " cells protruding from the surface, writing cellSet "
            << protrudingCells.name()
            << endl;

        protrudingCells.write();
    }

    return protrudingCells;
}


void Foam::conformalVoronoiMesh::writePointPairs
(
    const fileName& fName
) const
{
    OBJstream os(fName);

    for
    (
        Delaunay::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        ++eit
    )
    {
        Cell_handle c = eit->first;
        Vertex_handle vA = c->vertex(eit->second);
        Vertex_handle vB = c->vertex(eit->third);

        if (ptPairs_.isPointPair(vA, vB))
        {
            os.write
            (
                linePointRef(topoint(vA->point()), topoint(vB->point()))
            );
        }
    }
}


// ************************************************************************* //
