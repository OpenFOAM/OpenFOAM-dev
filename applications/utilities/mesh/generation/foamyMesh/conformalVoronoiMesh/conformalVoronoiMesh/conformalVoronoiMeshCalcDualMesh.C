/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
#include "motionSmoother.H"
#include "backgroundMeshDecomposition.H"
#include "polyMeshGeometry.H"
#include "indexedCellChecks.H"
#include "OBJstream.H"
#include "indexedCellOps.H"
#include "ListOps.H"
#include "DelaunayMeshTools.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::calcDualMesh
(
    pointField& points,
    labelList& boundaryPts,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchNames,
    PtrList<dictionary>& patchDicts,
    pointField& cellCentres,
    labelList& cellToDelaunayVertex,
    labelListList& patchToDelaunayVertex,
    PackedBoolList& boundaryFacesToRemove
)
{
    timeCheck("Start calcDualMesh");

    setVertexSizeAndAlignment();

    timeCheck("After setVertexSizeAndAlignment");

    indexDualVertices(points, boundaryPts);

    {
        Info<< nl << "Merging identical points" << endl;

        // There is no guarantee that a merge of close points is no-risk
        mergeIdenticalDualVertices(points, boundaryPts);
    }

    // Final dual face and owner neighbour construction

    timeCheck("Before createFacesOwnerNeighbourAndPatches");

    createFacesOwnerNeighbourAndPatches
    (
        points,
        faces,
        owner,
        neighbour,
        patchNames,
        patchDicts,
        patchToDelaunayVertex,  // from patch face to Delaunay vertex (slavePp)
        boundaryFacesToRemove,
        false
    );

    // deferredCollapseFaceSet(owner, neighbour, deferredCollapseFaces);

    cellCentres = DelaunayMeshTools::allPoints(*this);

    cellToDelaunayVertex = removeUnusedCells(owner, neighbour);

    cellCentres = pointField(cellCentres, cellToDelaunayVertex);

    removeUnusedPoints(faces, points, boundaryPts);

    timeCheck("End of calcDualMesh");
}


void Foam::conformalVoronoiMesh::calcTetMesh
(
    pointField& points,
    labelList& pointToDelaunayVertex,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchNames,
    PtrList<dictionary>& patchDicts
)
{
    labelList vertexMap(number_of_vertices());

    label vertI = 0;

    points.setSize(number_of_vertices());
    pointToDelaunayVertex.setSize(number_of_vertices());

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint() || vit->boundaryPoint())
        {
            vertexMap[vit->index()] = vertI;
            points[vertI] = topoint(vit->point());
            pointToDelaunayVertex[vertI] = vit->index();
            vertI++;
        }
    }

    points.setSize(vertI);
    pointToDelaunayVertex.setSize(vertI);

    label celli = 0;

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if (cit->internalOrBoundaryDualVertex())
        {
             cit->cellIndex() = celli++;
        }
        else
        {
            cit->cellIndex() = Cb::ctFar;
        }
    }

    patchNames = geometryToConformTo_.patchNames();

    patchNames.setSize(patchNames.size() + 1);

    patchNames[patchNames.size() - 1] = "foamyHexMesh_defaultPatch";

    label nPatches = patchNames.size();

    List<DynamicList<face>> patchFaces(nPatches, DynamicList<face>(0));

    List<DynamicList<label>> patchOwners(nPatches, DynamicList<label>(0));

    faces.setSize(number_of_finite_facets());

    owner.setSize(number_of_finite_facets());

    neighbour.setSize(number_of_finite_facets());

    label facei = 0;

    labelList verticesOnTriFace(3, label(-1));

    face newFace(verticesOnTriFace);

    for
    (
        Delaunay::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle c1(fit->first);
        const label oppositeVertex = fit->second;
        const Cell_handle c2(c1->neighbor(oppositeVertex));

        if (c1->hasFarPoint() && c2->hasFarPoint())
        {
            // Both tets are outside, skip
            continue;
        }

        label c1I = c1->cellIndex();
        label c2I = c2->cellIndex();

        label ownerCell = -1;
        label neighbourCell = -1;

        for (label i = 0; i < 3; i++)
        {
            verticesOnTriFace[i] = vertexMap
            [
                c1->vertex(vertex_triple_index(oppositeVertex, i))->index()
            ];
        }

        newFace = face(verticesOnTriFace);

        if (c1->hasFarPoint() || c2->hasFarPoint())
        {
            // Boundary face...
            if (c1->hasFarPoint())
            {
                //... with c1 outside
                ownerCell = c2I;
            }
            else
            {
                // ... with c2 outside
                ownerCell = c1I;

                reverse(newFace);
            }

            label patchIndex = geometryToConformTo_.findPatch
            (
                newFace.centre(points)
            );

            if (patchIndex == -1)
            {
                patchIndex = patchNames.size() - 1;

                WarningInFunction
                    << newFace.centre(points) << nl
                    << "did not find a surface patch. Adding to "
                    << patchNames[patchIndex]
                    << endl;
            }

            patchFaces[patchIndex].append(newFace);
            patchOwners[patchIndex].append(ownerCell);
        }
        else
        {
            // Internal face...
            if (c1I < c2I)
            {
                // ...with c1 as the ownerCell
                ownerCell = c1I;
                neighbourCell = c2I;

                reverse(newFace);
            }
            else
            {
                // ...with c2 as the ownerCell
                ownerCell = c2I;
                neighbourCell = c1I;
            }

            faces[facei] = newFace;
            owner[facei] = ownerCell;
            neighbour[facei] = neighbourCell;
            facei++;
        }
    }

    label nInternalFaces = facei;

    faces.setSize(nInternalFaces);
    owner.setSize(nInternalFaces);
    neighbour.setSize(nInternalFaces);

    sortFaces(faces, owner, neighbour);

//    PackedBoolList boundaryFacesToRemove;
//    List<DynamicList<bool>> indirectPatchFace;
//
//    addPatches
//    (
//        nInternalFaces,
//        faces,
//        owner,
//        patchDicts,
//        boundaryFacesToRemove,
//        patchFaces,
//        patchOwners,
//        indirectPatchFace
//    );
}


void Foam::conformalVoronoiMesh::mergeIdenticalDualVertices
(
    const pointField& pts,
    labelList& boundaryPts
)
{
    // Assess close points to be merged

    label nPtsMerged = 0;
    label nPtsMergedSum = 0;

    do
    {
        Map<label> dualPtIndexMap;

        nPtsMerged = mergeIdenticalDualVertices
        (
            pts,
            dualPtIndexMap
        );

        reindexDualVertices(dualPtIndexMap, boundaryPts);

        reduce(nPtsMerged, sumOp<label>());

        nPtsMergedSum += nPtsMerged;

    } while (nPtsMerged > 0);

    if (nPtsMergedSum > 0)
    {
        Info<< "    Merged " << nPtsMergedSum << " points " << endl;
    }
}


Foam::label Foam::conformalVoronoiMesh::mergeIdenticalDualVertices
(
    const pointField& pts,
    Map<label>& dualPtIndexMap
) const
{
    label nPtsMerged = 0;

    for
    (
        Delaunay::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        const Cell_handle c1(fit->first);
        const label oppositeVertex = fit->second;
        const Cell_handle c2(c1->neighbor(oppositeVertex));

        if (is_infinite(c1) || is_infinite(c2))
        {
            continue;
        }

        label& c1I = c1->cellIndex();
        label& c2I = c2->cellIndex();

        if ((c1I != c2I) && !c1->hasFarPoint() && !c2->hasFarPoint())
        {
            const Foam::point& p1 = pts[c1I];
            const Foam::point& p2 = pts[c2I];

            if (p1 == p2)
            {
//                if (c1->parallelDualVertex() || c2->parallelDualVertex())
//                {
//                    if (c1->vertexLowestProc() < c2->vertexLowestProc())
//                    {
//                        dualPtIndexMap.insert(c1I, c1I);
//                        dualPtIndexMap.insert(c2I, c1I);
//                    }
//                    else
//                    {
//                        dualPtIndexMap.insert(c1I, c2I);
//                        dualPtIndexMap.insert(c2I, c2I);
//                    }
//                }
                if (c1I < c2I)
                {
                    dualPtIndexMap.insert(c1I, c1I);
                    dualPtIndexMap.insert(c2I, c1I);
                }
                else
                {
                    dualPtIndexMap.insert(c1I, c2I);
                    dualPtIndexMap.insert(c2I, c2I);
                }

                nPtsMerged++;
            }
        }
    }

    if (debug)
    {
        Info<< "mergeIdenticalDualVertices:" << endl
            << "    zero-length edges     : "
            << returnReduce(nPtsMerged, sumOp<label>()) << endl
            << endl;
    }

    return nPtsMerged;
}


//void Foam::conformalVoronoiMesh::smoothSurface
//(
//    pointField& pts,
//    const labelList& boundaryPts
//)
//{
//    label nCollapsedFaces = 0;
//
//    label iterI = 0;
//
//    do
//    {
//        Map<label> dualPtIndexMap;
//
//        nCollapsedFaces = smoothSurfaceDualFaces
//        (
//            pts,
//            boundaryPts,
//            dualPtIndexMap
//        );
//
//        reduce(nCollapsedFaces, sumOp<label>());
//
//        reindexDualVertices(dualPtIndexMap);
//
//        mergeIdenticalDualVertices(pts, boundaryPts);
//
//        if (nCollapsedFaces > 0)
//        {
//            Info<< "    Collapsed " << nCollapsedFaces << " boundary faces"
//                << endl;
//        }
//
//        if (++iterI > foamyHexMeshControls().maxCollapseIterations())
//        {
//            Info<< "    maxCollapseIterations reached, stopping collapse"
//                << endl;
//
//            break;
//        }
//
//    } while (nCollapsedFaces > 0);
//
//    // Force all points of boundary faces to be on the surface
////    for
////    (
////        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
////        cit != finite_cells_end();
////        ++cit
////    )
////    {
////        label ptI = cit->cellIndex();
////
////        label fC = cit->filterCount();
////
////        if (fC > foamyHexMeshControls().filterCountSkipThreshold())
////        {
////            // This vertex has been limited too many times, skip
////            continue;
////        }
////
////        // Only cells with indices > -1 are valid
////        if (ptI > -1)
////        {
////            if (boundaryPts[ptI] != -1)
////            {
////                Foam::point& pt = pts[ptI];
////
////                pointIndexHit surfHit;
////                label hitSurface;
////
////                geometryToConformTo_.findSurfaceNearest
////                (
////                    pt,
////                    sqr(great),
////                    surfHit,
////                    hitSurface
////                );
////
////                if (surfHit.hit())
////                {
////                    pt += (surfHit.hitPoint() - pt)
////                         *pow
////                          (
////                              foamyHexMeshControls()
////                                .filterErrorReductionCoeff(),
////                              fC
////                          );
////                }
////            }
////        }
////    }
////
////    mergeCloseDualVertices(pts, boundaryPts);
//}
//
//
//Foam::label Foam::conformalVoronoiMesh::smoothSurfaceDualFaces
//(
//    pointField& pts,
//    const labelList& boundaryPts,
//    Map<label>& dualPtIndexMap
//) const
//{
//    label nCollapsedFaces = 0;
//
//    const scalar cosPerpendicularToleranceAngle = cos
//    (
//        degToRad(foamyHexMeshControls().surfaceStepFaceAngle())
//    );
//
//    for
//    (
//        Delaunay::Finite_edges_iterator eit = finite_edges_begin();
//        eit != finite_edges_end();
//        ++eit
//    )
//    {
//        Cell_circulator ccStart = incident_cells(*eit);
//        Cell_circulator cc = ccStart;
//
//        bool skipFace = false;
//
//        do
//        {
//            if (dualPtIndexMap.found(cc->cellIndex()))
//            {
//                // One of the points of this face has already been
//                // collapsed this sweep, leave for next sweep
//
//                skipFace = true;
//
//                break;
//            }
//
//        } while (++cc != ccStart);
//
//        if (skipFace)
//        {
//            continue;
//        }
//
//        if (isBoundaryDualFace(eit))
//        {
//            face dualFace = buildDualFace(eit);
//
//            if (dualFace.size() < 3)
//            {
//                // This face has been collapsed already
//                continue;
//            }
//
//            label maxFC = maxFilterCount(eit);
//
//            if (maxFC > foamyHexMeshControls().filterCountSkipThreshold())
//            {
//                // A vertex on this face has been limited too many
//                // times, skip
//                continue;
//            }
//
//            Cell_handle c = eit->first;
//            Vertex_handle vA = c->vertex(eit->second);
//            Vertex_handle vB = c->vertex(eit->third);
//
//            if
//            (
//                vA->internalBoundaryPoint() && vA->surfacePoint()
//             && vB->externalBoundaryPoint() && vB->surfacePoint()
//            )
//            {
//                if (vA->index() == vB->index() - 1)
//                {
//                    continue;
//                }
//            }
//            else if
//            (
//                vA->externalBoundaryPoint() && vA->surfacePoint()
//             && vB->internalBoundaryPoint() && vB->surfacePoint()
//            )
//            {
//                if (vA->index() == vB->index() + 1)
//                {
//                    continue;
//                }
//            }
////            else if
////            (
////                vA->internalBoundaryPoint() && vA->featureEdgePoint()
////             && vB->externalBoundaryPoint() && vB->featureEdgePoint()
////            )
////            {
////                if (vA->index() == vB->index() - 1)
////                {
////                    continue;
////                }
////            }
////            else if
////            (
////                vA->externalBoundaryPoint() && vA->featureEdgePoint()
////             && vB->internalBoundaryPoint() && vB->featureEdgePoint()
////            )
////            {
////                if (vA->index() == vB->index() + 1)
////                {
////                    continue;
////                }
////            }
////            else if
////            (
////                vA->internalBoundaryPoint() && vA->featurePoint()
////             && vB->externalBoundaryPoint() && vB->featurePoint()
////            )
////            {
////                if (vA->index() == vB->index() - 1)
////                {
////                    continue;
////                }
////            }
////            else if
////            (
////                vA->externalBoundaryPoint() && vA->featurePoint()
////             && vB->internalBoundaryPoint() && vB->featurePoint()
////            )
////            {
////                if (vA->index() == vB->index() + 1)
////                {
////                    continue;
////                }
////            }
//
//
//            if ((faceNormal & surfaceNormal) < cosPerpendicularToleranceAngle)
//            {
//                scalar targetFaceSize = averageAnyCellSize(vA, vB);
//
//                // Selecting faces to collapse based on angle to
//                // surface, so set collapseSizeLimitCoeff to great to
//                // allow collapse of all faces
//
//                faceCollapseMode mode = collapseFace
//                (
//                    dualFace,
//                    pts,
//                    boundaryPts,
//                    dualPtIndexMap,
//                    targetFaceSize,
//                    great,
//                    maxFC
//                );
//
//                if (mode == fcmPoint || mode == fcmEdge)
//                {
//                    nCollapsedFaces++;
//                }
//            }
//        }
//    }
//
//    return nCollapsedFaces;
//}


void Foam::conformalVoronoiMesh::deferredCollapseFaceSet
(
    labelList& owner,
    labelList& neighbour,
    const HashSet<labelPair, labelPair::Hash<>>& deferredCollapseFaces
) const
{
    DynamicList<label> faceLabels;

    forAll(neighbour, nI)
    {
        if (deferredCollapseFaces.found(Pair<label>(owner[nI], neighbour[nI])))
        {
            faceLabels.append(nI);
        }
    }

    Pout<< "facesToCollapse" << nl << faceLabels << endl;
}


Foam::autoPtr<Foam::polyMesh>
Foam::conformalVoronoiMesh::createPolyMeshFromPoints
(
    const pointField& pts
) const
{
    faceList faces;
    labelList owner;
    labelList neighbour;
    wordList patchNames;
    PtrList<dictionary> patchDicts;
    pointField cellCentres;
    labelListList patchToDelaunayVertex;
    PackedBoolList boundaryFacesToRemove;

    timeCheck("Start of checkPolyMeshQuality");

    Info<< nl << "Creating polyMesh to assess quality" << endl;

    createFacesOwnerNeighbourAndPatches
    (
        pts,
        faces,
        owner,
        neighbour,
        patchNames,
        patchDicts,
        patchToDelaunayVertex,
        boundaryFacesToRemove,
        false
    );

    cellCentres = DelaunayMeshTools::allPoints(*this);

    labelList cellToDelaunayVertex(removeUnusedCells(owner, neighbour));
    cellCentres = pointField(cellCentres, cellToDelaunayVertex);

    autoPtr<polyMesh> meshPtr
    (
        new polyMesh
        (
            IOobject
            (
                "foamyHexMesh_temporary",
                runTime_.timeName(),
                runTime_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            xferCopy(pts),
            xferMove(faces),
            xferMove(owner),
            xferMove(neighbour)
        )
    );

    polyMesh& pMesh = meshPtr();

    List<polyPatch*> patches(patchNames.size());

    label nValidPatches = 0;

    forAll(patches, p)
    {
        label totalPatchSize = readLabel(patchDicts[p].lookup("nFaces"));

        if
        (
            patchDicts.set(p)
         && word(patchDicts[p].lookup("type")) == processorPolyPatch::typeName
        )
        {
            // Do not create empty processor patches
            if (totalPatchSize > 0)
            {
                patchDicts[p].set("transform", "coincidentFullMatch");

                patches[nValidPatches] = new processorPolyPatch
                (
                    patchNames[p],
                    patchDicts[p],
                    nValidPatches,
                    pMesh.boundaryMesh(),
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
                    pMesh.boundaryMesh()
                ).ptr();

                nValidPatches++;
            }
        }
    }

    patches.setSize(nValidPatches);

    pMesh.addPatches(patches);

    return meshPtr;
}


void Foam::conformalVoronoiMesh::checkCellSizing()
{
    Info<< "Checking cell sizes..."<< endl;

    timeCheck("Start of Cell Sizing");

    labelList boundaryPts(number_of_finite_cells(), internal);
    pointField ptsField;

    indexDualVertices(ptsField, boundaryPts);

    // Merge close dual vertices.
    mergeIdenticalDualVertices(ptsField, boundaryPts);

    autoPtr<polyMesh> meshPtr = createPolyMeshFromPoints(ptsField);
    const polyMesh& pMesh = meshPtr();

    // pMesh.write();

    // Find cells with poor quality
    DynamicList<label> checkFaces(identity(pMesh.nFaces()));
    labelHashSet wrongFaces(pMesh.nFaces()/100);

    Info<< "Running checkMesh on mesh with " << pMesh.nCells()
        << " cells "<< endl;

    const dictionary& dict
        = foamyHexMeshControls().foamyHexMeshDict();

    const dictionary& meshQualityDict
        = dict.subDict("meshQualityControls");

    const scalar maxNonOrtho =
        readScalar(meshQualityDict.lookup("maxNonOrtho", true));

    label nWrongFaces = 0;

    if (maxNonOrtho < 180.0 - small)
    {
        polyMeshGeometry::checkFaceDotProduct
        (
            false,
            maxNonOrtho,
            pMesh,
            pMesh.cellCentres(),
            pMesh.faceAreas(),
            checkFaces,
            List<labelPair>(),
            &wrongFaces
        );

        label nNonOrthogonal = returnReduce(wrongFaces.size(), sumOp<label>());

        Info<< "    non-orthogonality > " << maxNonOrtho
            << " degrees : " << nNonOrthogonal << endl;

        nWrongFaces += nNonOrthogonal;
    }

    labelHashSet protrudingCells = findOffsetPatchFaces(pMesh, 0.25);

    label nProtrudingCells = protrudingCells.size();

    Info<< "    protruding/intruding cells : " << nProtrudingCells << endl;

    nWrongFaces += nProtrudingCells;

//    motionSmoother::checkMesh
//    (
//        false,
//        pMesh,
//        meshQualityDict,
//        checkFaces,
//        wrongFaces
//    );

    Info<< "    Found total of " << nWrongFaces << " bad faces" << endl;

    {
        labelHashSet cellsToResizeMap(pMesh.nFaces()/100);

        // Find cells that are attached to the faces in wrongFaces.
        forAllConstIter(labelHashSet, wrongFaces, iter)
        {
            const label faceOwner = pMesh.faceOwner()[iter.key()];
            const label faceNeighbour = pMesh.faceNeighbour()[iter.key()];

            if (!cellsToResizeMap.found(faceOwner))
            {
                cellsToResizeMap.insert(faceOwner);
            }

            if (!cellsToResizeMap.found(faceNeighbour))
            {
                cellsToResizeMap.insert(faceNeighbour);
            }
        }

        cellsToResizeMap += protrudingCells;

        pointField cellsToResize(cellsToResizeMap.size());

        label count = 0;
        for (label celli = 0; celli < pMesh.nCells(); ++celli)
        {
            if (cellsToResizeMap.found(celli))
            {
                cellsToResize[count++] = pMesh.cellCentres()[celli];
            }
        }

        Info<< "    DISABLED: Automatically re-sizing " << cellsToResize.size()
            << " cells that are attached to the bad faces: " << endl;

        // cellSizeControl_.setCellSizes(cellsToResize);
    }

    timeCheck("End of Cell Sizing");

    Info<< "Finished checking cell sizes"<< endl;
}


Foam::labelHashSet Foam::conformalVoronoiMesh::findOffsetPatchFaces
(
    const polyMesh& mesh,
    const scalar allowedOffset
) const
{
    timeCheck("Start findRemainingProtrusionSet");

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    cellSet offsetBoundaryCells
    (
        mesh,
        "foamyHexMesh_protrudingCells",
        mesh.nCells()/1000
    );

    forAll(patches, patchi)
    {
        const polyPatch& patch = patches[patchi];

        const faceList& localFaces = patch.localFaces();
        const pointField& localPoints = patch.localPoints();

        const labelList& fCell = patch.faceCells();

        forAll(localFaces, pLFI)
        {
            const face& f = localFaces[pLFI];

            const Foam::point& faceCentre = f.centre(localPoints);

            const scalar targetSize = targetCellSize(faceCentre);

            pointIndexHit pHit;
            label surfHit = -1;

            geometryToConformTo_.findSurfaceNearest
            (
                faceCentre,
                sqr(targetSize),
                pHit,
                surfHit
            );

            if
            (
                pHit.hit()
             && (mag(pHit.hitPoint() - faceCentre) > allowedOffset*targetSize)
            )
            {
                offsetBoundaryCells.insert(fCell[pLFI]);
            }
        }
    }

    if (foamyHexMeshControls().objOutput())
    {
        offsetBoundaryCells.write();
    }

    return offsetBoundaryCells;
}


Foam::labelHashSet Foam::conformalVoronoiMesh::checkPolyMeshQuality
(
    const pointField& pts
) const
{
    autoPtr<polyMesh> meshPtr = createPolyMeshFromPoints(pts);
    polyMesh& pMesh = meshPtr();

    timeCheck("polyMesh created, checking quality");

    labelHashSet wrongFaces(pMesh.nFaces()/100);

    DynamicList<label> checkFaces(pMesh.nFaces());

    const vectorField& fAreas = pMesh.faceAreas();

    scalar faceAreaLimit = small;

    forAll(fAreas, fI)
    {
        if (mag(fAreas[fI]) > faceAreaLimit)
        {
            checkFaces.append(fI);
        }
    }

    Info<< nl << "Excluding "
        << returnReduce(fAreas.size() - checkFaces.size(), sumOp<label>())
        << " faces from check, < " << faceAreaLimit << " area" << endl;

    const dictionary& dict
        = foamyHexMeshControls().foamyHexMeshDict();

    const dictionary& meshQualityDict
        = dict.subDict("meshQualityControls");

    motionSmoother::checkMesh
    (
        false,
        pMesh,
        meshQualityDict,
        checkFaces,
        wrongFaces
    );

    {
        // Check for cells with more than 1 but fewer than 4 faces
        label nInvalidPolyhedra = 0;

        const cellList& cells = pMesh.cells();

        forAll(cells, cI)
        {
            if (cells[cI].size() < 4 && cells[cI].size() > 0)
            {
                // Pout<< "cell " << cI << " " << cells[cI]
                //     << " has " << cells[cI].size() << " faces."
                //     << endl;

                nInvalidPolyhedra++;

                forAll(cells[cI], cFI)
                {
                    wrongFaces.insert(cells[cI][cFI]);
                }
            }
        }

        Info<< "    cells with more than 1 but fewer than 4 faces          : "
            << returnReduce(nInvalidPolyhedra, sumOp<label>())
            << endl;

        // Check for cells with one internal face only

        labelList nInternalFaces(pMesh.nCells(), label(0));

        for (label fI = 0; fI < pMesh.nInternalFaces(); fI++)
        {
            nInternalFaces[pMesh.faceOwner()[fI]]++;
            nInternalFaces[pMesh.faceNeighbour()[fI]]++;
        }

        const polyBoundaryMesh& patches = pMesh.boundaryMesh();

        forAll(patches, patchi)
        {
            if (patches[patchi].coupled())
            {
                const labelUList& owners = patches[patchi].faceCells();

                forAll(owners, i)
                {
                    nInternalFaces[owners[i]]++;
                }
            }
        }

        label oneInternalFaceCells = 0;

        forAll(nInternalFaces, cI)
        {
            if (nInternalFaces[cI] <= 1)
            {
                oneInternalFaceCells++;

                forAll(cells[cI], cFI)
                {
                    wrongFaces.insert(cells[cI][cFI]);
                }
            }
        }

        Info<< "    cells with with zero or one non-boundary face          : "
            << returnReduce(oneInternalFaceCells, sumOp<label>())
            << endl;
    }


    PackedBoolList ptToBeLimited(pts.size(), false);

    forAllConstIter(labelHashSet, wrongFaces, iter)
    {
        const face f = pMesh.faces()[iter.key()];

        forAll(f, fPtI)
        {
            ptToBeLimited[f[fPtI]] = true;
        }
    }

    // // Limit connected cells

    // labelHashSet limitCells(pMesh.nCells()/100);

    // const labelListList& ptCells = pMesh.pointCells();

    // forAllConstIter(labelHashSet, wrongFaces, iter)
    // {
    //     const face f = pMesh.faces()[iter.key()];

    //     forAll(f, fPtI)
    //     {
    //         label ptI = f[fPtI];

    //         const labelList& pC = ptCells[ptI];

    //         forAll(pC, pCI)
    //         {
    //             limitCells.insert(pC[pCI]);
    //         }
    //     }
    // }

    // const labelListList& cellPts = pMesh.cellPoints();

    // forAllConstIter(labelHashSet, limitCells, iter)
    // {
    //     label celli = iter.key();

    //     const labelList& cP = cellPts[celli];

    //     forAll(cP, cPI)
    //     {
    //         ptToBeLimited[cP[cPI]] = true;
    //     }
    // }


    // Apply Delaunay cell filterCounts and determine the maximum
    // overall filterCount

    label maxFilterCount = 0;

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        label cI = cit->cellIndex();

        if (cI >= 0)
        {
            if (ptToBeLimited[cI] == true)
            {
                cit->filterCount()++;
            }

            if (cit->filterCount() > maxFilterCount)
            {
                maxFilterCount = cit->filterCount();
            }
        }
    }

    Info<< nl << "Maximum number of filter limits applied: "
        << returnReduce(maxFilterCount, maxOp<label>()) << endl;

    return wrongFaces;
}


Foam::label Foam::conformalVoronoiMesh::classifyBoundaryPoint
(
    Cell_handle cit
) const
{
    if (cit->boundaryDualVertex())
    {
        if (cit->featurePointDualVertex())
        {
            return featurePoint;
        }
        else if (cit->featureEdgeDualVertex())
        {
            return featureEdge;
        }
        else
        {
            return surface;
        }
    }
    else if (cit->baffleSurfaceDualVertex())
    {
        return surface;
    }
    else if (cit->baffleEdgeDualVertex())
    {
        return featureEdge;
    }
    else
    {
        return internal;
    }
}


void Foam::conformalVoronoiMesh::indexDualVertices
(
    pointField& pts,
    labelList& boundaryPts
)
{
    // Indexing Delaunay cells, which are the dual vertices

    this->resetCellCount();

    label nConstrainedVertices = 0;
    if (foamyHexMeshControls().guardFeaturePoints())
    {
        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            ++vit
        )
        {
            if (vit->constrained())
            {
                vit->index() = number_of_finite_cells() + nConstrainedVertices;
                nConstrainedVertices++;
            }
        }
    }

    pts.setSize(number_of_finite_cells() + nConstrainedVertices);
    boundaryPts.setSize
    (
        number_of_finite_cells() + nConstrainedVertices,
        internal
    );

    if (foamyHexMeshControls().guardFeaturePoints())
    {
        nConstrainedVertices = 0;
        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            ++vit
        )
        {
            if (vit->constrained())
            {
                pts[number_of_finite_cells() + nConstrainedVertices] =
                    topoint(vit->point());

                boundaryPts[number_of_finite_cells() + nConstrainedVertices] =
                    constrained;

                nConstrainedVertices++;
            }
        }
    }

    // OBJstream snapping1("snapToSurface1.obj");
    // OBJstream snapping2("snapToSurface2.obj");
    // OFstream tetToSnapTo("tetsToSnapTo.obj");

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
//        if (tetrahedron(cit).volume() == 0)
//        {
//            Pout<< "ZERO VOLUME TET" << endl;
//            Pout<< cit->info();
//            Pout<< "Dual = " << cit->dual();
//        }

        if (!cit->hasFarPoint())
        {
            cit->cellIndex() = getNewCellIndex();

            // For nearly coplanar Delaunay cells that are present on different
            // processors the result of the circumcentre calculation depends on
            // the ordering of the vertices, so synchronise it across processors

            if (Pstream::parRun() && cit->parallelDualVertex())
            {
                typedef CGAL::Exact_predicates_exact_constructions_kernel Exact;
                typedef CGAL::Point_3<Exact> ExactPoint;

                List<labelPair> cellVerticesPair(4);
                List<ExactPoint> cellVertices(4);

                for (label vI = 0; vI < 4; ++vI)
                {
                    cellVerticesPair[vI] = labelPair
                    (
                        cit->vertex(vI)->procIndex(),
                        cit->vertex(vI)->index()
                    );

                    cellVertices[vI] = ExactPoint
                    (
                        cit->vertex(vI)->point().x(),
                        cit->vertex(vI)->point().y(),
                        cit->vertex(vI)->point().z()
                    );
                }

                // Sort the vertices so that they will be in the same order on
                // each processor
                labelList oldToNew;
                sortedOrder(cellVerticesPair, oldToNew);
                oldToNew = invert(oldToNew.size(), oldToNew);
                inplaceReorder(oldToNew, cellVertices);

                ExactPoint synchronisedDual = CGAL::circumcenter
                (
                    cellVertices[0],
                    cellVertices[1],
                    cellVertices[2],
                    cellVertices[3]
                );

                pts[cit->cellIndex()] = Foam::point
                (
                    CGAL::to_double(synchronisedDual.x()),
                    CGAL::to_double(synchronisedDual.y()),
                    CGAL::to_double(synchronisedDual.z())
                );
            }
            else
            {
                pts[cit->cellIndex()] = cit->dual();
            }

            // Feature point snapping
            if (foamyHexMeshControls().snapFeaturePoints())
            {
                if (cit->featurePointDualVertex())
                {
                    pointFromPoint dual = cit->dual();

                    pointIndexHit fpHit;
                    label featureHit;

                    // Find nearest feature point and compare
                    geometryToConformTo_.findFeaturePointNearest
                    (
                        dual,
                        sqr(targetCellSize(dual)),
                        fpHit,
                        featureHit
                    );

                    if (fpHit.hit())
                    {
                        if (debug)
                        {
                            Info<< "Dual        = " << dual << nl
                                << "    Nearest = " << fpHit.hitPoint() << endl;
                        }

                        pts[cit->cellIndex()] = fpHit.hitPoint();
                    }
                }
            }

//            {
//                // Snapping points far outside
//                if (cit->boundaryDualVertex() && !cit->parallelDualVertex())
//                {
//                    pointFromPoint dual = cit->dual();
//
//                    pointIndexHit hitInfo;
//                    label surfHit;
//
//                    // Find nearest surface point
//                    geometryToConformTo_.findSurfaceNearest
//                    (
//                        dual,
//                        sqr(targetCellSize(dual)),
//                        hitInfo,
//                        surfHit
//                    );
//
//                    if (!hitInfo.hit())
//                    {
//                        // Project dual to nearest point on tet
//
//                        tetPointRef tet
//                        (
//                            topoint(cit->vertex(0)->point()),
//                            topoint(cit->vertex(1)->point()),
//                            topoint(cit->vertex(2)->point()),
//                            topoint(cit->vertex(3)->point())
//                        );
//
//                        pointFromPoint nearestPointOnTet =
//                            tet.nearestPoint(dual).rawPoint();
//
//                        // Get nearest point on surface from tet.
//                        geometryToConformTo_.findSurfaceNearest
//                        (
//                            nearestPointOnTet,
//                            sqr(targetCellSize(nearestPointOnTet)),
//                            hitInfo,
//                            surfHit
//                        );
//
//                        vector snapDir = nearestPointOnTet - dual;
//                        snapDir /= mag(snapDir) + small;
//
//                        drawDelaunayCell(tetToSnapTo, cit, offset);
//                        offset += 1;
//
//                        vectorField norm(1);
//                        allGeometry_[surfHit].getNormal
//                        (
//                            List<pointIndexHit>(1, hitInfo),
//                            norm
//                        );
//                        norm[0] /= mag(norm[0]) + small;
//
//                        if
//                        (
//                            hitInfo.hit()
//                         && (mag(snapDir & norm[0]) > 0.5)
//                        )
//                        {
//                            snapping1.write
//                            (
//                                linePointRef(dual, nearestPointOnTet)
//                            );
//
//                            snapping2.write
//                            (
//                                linePointRef
//                                (
//                                    nearestPointOnTet,
//                                    hitInfo.hitPoint()
//                                )
//                            );
//
//                            pts[cit->cellIndex()] = hitInfo.hitPoint();
//                        }
//                    }
//                }
//            }

            boundaryPts[cit->cellIndex()] = classifyBoundaryPoint(cit);
        }
        else
        {
            cit->cellIndex() = Cb::ctFar;
        }
    }

    // pts.setSize(this->cellCount());

    // boundaryPts.setSize(this->cellCount());
}


void Foam::conformalVoronoiMesh::reindexDualVertices
(
    const Map<label>& dualPtIndexMap,
    labelList& boundaryPts
)
{
    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        if (dualPtIndexMap.found(cit->cellIndex()))
        {
            cit->cellIndex() = dualPtIndexMap[cit->cellIndex()];
            boundaryPts[cit->cellIndex()] =
                max
                (
                    boundaryPts[cit->cellIndex()],
                    boundaryPts[dualPtIndexMap[cit->cellIndex()]]
                );
        }
    }
}


Foam::label Foam::conformalVoronoiMesh::createPatchInfo
(
    wordList& patchNames,
    PtrList<dictionary>& patchDicts
) const
{
    patchNames = geometryToConformTo_.patchNames();

    patchDicts.setSize(patchNames.size() + 1);

    const PtrList<dictionary>& patchInfo = geometryToConformTo_.patchInfo();

    forAll(patchNames, patchi)
    {
        if (patchInfo.set(patchi))
        {
            patchDicts.set(patchi, new dictionary(patchInfo[patchi]));
        }
        else
        {
            patchDicts.set(patchi, new dictionary());
            patchDicts[patchi].set
            (
                "type",
                wallPolyPatch::typeName
            );
        }
    }

    patchNames.setSize(patchNames.size() + 1);
    label defaultPatchIndex = patchNames.size() - 1;
    patchNames[defaultPatchIndex] = "foamyHexMesh_defaultPatch";
    patchDicts.set(defaultPatchIndex, new dictionary());
    patchDicts[defaultPatchIndex].set
    (
        "type",
        wallPolyPatch::typeName
    );

    label nProcPatches = 0;

    if (Pstream::parRun())
    {
        List<boolList> procUsedList
        (
            Pstream::nProcs(),
            boolList(Pstream::nProcs(), false)
        );

        boolList& procUsed = procUsedList[Pstream::myProcNo()];

        // Determine which processor patches are required
        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            vit++
        )
        {
            // This test is not sufficient if one of the processors does
            // not receive a referred vertex from another processor, but does
            // send one to the other processor.
            if (vit->referred())
            {
                procUsed[vit->procIndex()] = true;
            }
        }

        // Because the previous test was insufficient, combine the lists.
        Pstream::gatherList(procUsedList);
        Pstream::scatterList(procUsedList);

        forAll(procUsedList, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                if (procUsedList[proci][Pstream::myProcNo()])
                {
                    procUsed[proci] = true;
                }
            }
        }

        forAll(procUsed, pUI)
        {
            if (procUsed[pUI])
            {
                nProcPatches++;
            }
        }

        label nNonProcPatches = patchNames.size();
        label nTotalPatches = nNonProcPatches + nProcPatches;

        patchNames.setSize(nTotalPatches);
        patchDicts.setSize(nTotalPatches);
        for (label pI = nNonProcPatches; pI < nTotalPatches; ++pI)
        {
            patchDicts.set(pI, new dictionary());
        }

        label procAddI = 0;

        forAll(procUsed, pUI)
        {
            if (procUsed[pUI])
            {
                patchNames[nNonProcPatches + procAddI] =
                    processorPolyPatch::newName(Pstream::myProcNo(), pUI);

                patchDicts[nNonProcPatches + procAddI].set
                (
                    "type",
                    processorPolyPatch::typeName
                );

                patchDicts[nNonProcPatches + procAddI].set
                (
                    "myProcNo",
                    Pstream::myProcNo()
                );

                patchDicts[nNonProcPatches + procAddI].set("neighbProcNo", pUI);

                procAddI++;
            }
        }
    }

    return defaultPatchIndex;
}


Foam::vector Foam::conformalVoronoiMesh::calcSharedPatchNormal
(
    Cell_handle c1,
    Cell_handle c2
) const
{
    List<Foam::point> patchEdge(2, point::max);

    // Get shared Facet
    for (label cI = 0; cI < 4; ++cI)
    {
        if (c1->neighbor(cI) != c2 && !c1->vertex(cI)->constrained())
        {
            if (c1->vertex(cI)->internalBoundaryPoint())
            {
                patchEdge[0] = topoint(c1->vertex(cI)->point());
            }
            else
            {
                patchEdge[1] = topoint(c1->vertex(cI)->point());
            }
        }
    }

    Info<< "    " << patchEdge << endl;

    return vector(patchEdge[1] - patchEdge[0]);
}


bool Foam::conformalVoronoiMesh::boundaryDualFace
(
    Cell_handle c1,
    Cell_handle c2
) const
{
    label nInternal = 0;
    label nExternal = 0;

    for (label cI = 0; cI < 4; ++cI)
    {
        if (c1->neighbor(cI) != c2 && !c1->vertex(cI)->constrained())
        {
            if (c1->vertex(cI)->internalBoundaryPoint())
            {
                nInternal++;
            }
            else if (c1->vertex(cI)->externalBoundaryPoint())
            {
                nExternal++;
            }
        }
    }

    Info<< "in = " << nInternal << " out = " << nExternal << endl;

    return (nInternal == 1 && nExternal == 1);
}


void Foam::conformalVoronoiMesh::createFacesOwnerNeighbourAndPatches
(
    const pointField& pts,
    faceList& faces,
    labelList& owner,
    labelList& neighbour,
    wordList& patchNames,
    PtrList<dictionary>& patchDicts,
    labelListList& patchPointPairSlaves,
    PackedBoolList& boundaryFacesToRemove,
    bool includeEmptyPatches
) const
{
    const label defaultPatchIndex = createPatchInfo(patchNames, patchDicts);

    const label nPatches = patchNames.size();

    labelList procNeighbours(nPatches, label(-1));
    forAll(procNeighbours, patchi)
    {
        if (patchDicts[patchi].found("neighbProcNo"))
        {
            procNeighbours[patchi] =
            (
                patchDicts[patchi].found("neighbProcNo")
              ? readLabel(patchDicts[patchi].lookup("neighbProcNo"))
              : -1
            );
        }
    }

    List<DynamicList<face>> patchFaces(nPatches, DynamicList<face>(0));
    List<DynamicList<label>> patchOwners(nPatches, DynamicList<label>(0));
    // Per patch face the index of the slave node of the point pair
    List<DynamicList<label>> patchPPSlaves(nPatches, DynamicList<label>(0));

    List<DynamicList<bool>> indirectPatchFace(nPatches, DynamicList<bool>(0));


    faces.setSize(number_of_finite_edges());
    owner.setSize(number_of_finite_edges());
    neighbour.setSize(number_of_finite_edges());
    boundaryFacesToRemove.setSize(number_of_finite_edges(), false);

    labelPairPairDynListList procPatchSortingIndex(nPatches);

    label dualFacei = 0;

    if (foamyHexMeshControls().guardFeaturePoints())
    {
        OBJstream startCellStr("startingCell.obj");
        OBJstream featurePointFacesStr("ftPtFaces.obj");
        OBJstream featurePointDualsStr("ftPtDuals.obj");
        OFstream cellStr("vertexCells.obj");

        label vcount = 1;

        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            ++vit
        )
        {
            if (vit->constrained())
            {
                // Find a starting cell
                std::list<Cell_handle> vertexCells;
                finite_incident_cells(vit, std::back_inserter(vertexCells));

                Cell_handle startCell;

                for
                (
                    std::list<Cell_handle>::iterator vcit = vertexCells.begin();
                    vcit != vertexCells.end();
                    ++vcit
                )
                {
                    if ((*vcit)->featurePointExternalCell())
                    {
                        startCell = *vcit;
                    }

                    if ((*vcit)->real())
                    {
                        featurePointDualsStr.write
                        (
                            linePointRef(topoint(vit->point()), (*vcit)->dual())
                        );
                    }
                }

                // Error if startCell is null
                if (startCell == nullptr)
                {
                    Pout<< "Start cell is null!" << endl;
                }

                // Need to pick a direction to walk in
                Cell_handle vc1 = startCell;
                Cell_handle vc2;

                Info<< "c1 index = " << vc1->cellIndex() << " "
                    << vc1->dual() << endl;

                for (label cI = 0; cI < 4; ++cI)
                {
                    Info<< "c1 = " << cI << " "
                        << vc1->neighbor(cI)->cellIndex() << " v = "
                        << vc1->neighbor(cI)->dual() << endl;

                    Info<< vc1->vertex(cI)->info();
                }

                Cell_handle nextCell;

                for (label cI = 0; cI < 4; ++cI)
                {
                    if (vc1->vertex(cI)->externalBoundaryPoint())
                    {
                        vc2 = vc1->neighbor(cI);

                        Info<< "    c2 is neighbor "
                            << vc2->cellIndex()
                            << " of c1" << endl;

                        for (label cI = 0; cI < 4; ++cI)
                        {
                            Info<< "    c2 = " << cI << " "
                                << vc2->neighbor(cI)->cellIndex() << " v = "
                                << vc2->vertex(cI)->index() << endl;
                        }

                        face f(3);
                        f[0] = vit->index();
                        f[1] = vc1->cellIndex();
                        f[2] = vc2->cellIndex();

                        Info<< "f " << f << endl;
                        forAll(f, pI)
                        {
                            Info<< "    " << pts[f[pI]] << endl;
                        }

                        vector correctNormal = calcSharedPatchNormal(vc1, vc2);
                        correctNormal /= mag(correctNormal);

                        Info<< "    cN " << correctNormal << endl;

                        vector fN = f.area(pts);

                        if (mag(fN) < small)
                        {
                            nextCell = vc2;
                            continue;
                        }

                        fN /= mag(fN);
                        Info<< "    fN " << fN << endl;

                        if ((fN & correctNormal) > 0)
                        {
                            nextCell = vc2;
                            break;
                        }
                    }
                }

                vc2 = nextCell;

                label own = vit->index();
                face f(3);
                f[0] = own;

                Info<< "Start walk from " << vc1->cellIndex()
                    << " to " << vc2->cellIndex() << endl;

                // Walk while not at start cell

                label iter = 0;
                do
                {
                    Info<< "     Walk from " << vc1->cellIndex()
                        << " " << vc1->dual()
                        << " to " << vc2->cellIndex()
                        << " " << vc2->dual()
                        << endl;

                    startCellStr.write(linePointRef(vc1->dual(), vc2->dual()));

                    // Get patch by getting face between cells and the two
                    // points on the face that are not the feature vertex
                    label patchIndex =
                        geometryToConformTo_.findPatch
                        (
                            topoint(vit->point())
                        );

                    f[1] = vc1->cellIndex();
                    f[2] = vc2->cellIndex();

                    patchFaces[patchIndex].append(f);
                    patchOwners[patchIndex].append(own);
                    patchPPSlaves[patchIndex].append(own);

                    // Find next cell
                    Cell_handle nextCell;

                    Info<< "    c1 vertices " << vc2->dual() << endl;
                    for (label cI = 0; cI < 4; ++cI)
                    {
                        Info<< "        " << vc2->vertex(cI)->info();
                    }
                    Info<< "    c1 neighbour vertices " << endl;
                    for (label cI = 0; cI < 4; ++cI)
                    {
                        if
                        (
                            !vc2->vertex(cI)->constrained()
                         && vc2->neighbor(cI) != vc1
                         && !is_infinite(vc2->neighbor(cI))
                         &&
                            (
                                vc2->neighbor(cI)->featurePointExternalCell()
                             || vc2->neighbor(cI)->featurePointInternalCell()
                            )
                         && vc2->neighbor(cI)->hasConstrainedPoint()
                        )
                        {
                            DelaunayMeshTools::drawDelaunayCell
                            (
                                cellStr,
                                vc2->neighbor(cI),
                                vcount++
                            );

                            Info<< "        neighbour " << cI << " "
                                << vc2->neighbor(cI)->dual() << endl;
                            for (label I = 0; I < 4; ++I)
                            {
                                Info<< "            "
                                    << vc2->neighbor(cI)->vertex(I)->info();
                            }
                        }
                    }

                    for (label cI = 0; cI < 4; ++cI)
                    {
                        if
                        (
                            !vc2->vertex(cI)->constrained()
                         && vc2->neighbor(cI) != vc1
                         && !is_infinite(vc2->neighbor(cI))
                         &&
                            (
                                vc2->neighbor(cI)->featurePointExternalCell()
                             || vc2->neighbor(cI)->featurePointInternalCell()
                            )
                         && vc2->neighbor(cI)->hasConstrainedPoint()
                        )
                        {
                            // check if shared edge is internal/internal
                            if (boundaryDualFace(vc2, vc2->neighbor(cI)))
                            {
                                nextCell = vc2->neighbor(cI);
                                break;
                            }
                        }
                    }

                    vc1 = vc2;
                    vc2 = nextCell;

                    iter++;
                } while (vc1 != startCell && iter < 100);
            }
        }
    }

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

        if (vA->constrained() && vB->constrained())
        {
            continue;
        }

        if
        (
            (vA->constrained() && vB->internalOrBoundaryPoint())
         || (vB->constrained() && vA->internalOrBoundaryPoint())
        )
        {
            face newDualFace = buildDualFace(eit);

            label own = -1;
            label nei = -1;

            if (ownerAndNeighbour(vA, vB, own, nei))
            {
                reverse(newDualFace);
            }

            // internal face
            faces[dualFacei] = newDualFace;
            owner[dualFacei] = own;
            neighbour[dualFacei] = nei;

            dualFacei++;
        }
        else if
        (
            (vA->internalOrBoundaryPoint() && !vA->referred())
         || (vB->internalOrBoundaryPoint() && !vB->referred())
        )
        {
            if
            (
                (vA->internalPoint() && vB->externalBoundaryPoint())
             || (vB->internalPoint() && vA->externalBoundaryPoint())
            )
            {
                Cell_circulator ccStart = incident_cells(*eit);
                Cell_circulator cc1 = ccStart;
                Cell_circulator cc2 = cc1;

                cc2++;

                bool skipEdge = false;

                do
                {
                    if
                    (
                        cc1->hasFarPoint() || cc2->hasFarPoint()
                     || is_infinite(cc1) || is_infinite(cc2)
                    )
                    {
                        Pout<< "Ignoring edge between internal and external: "
                            << vA->info()
                            << vB->info();

                        skipEdge = true;
                        break;
                    }

                    cc1++;
                    cc2++;

                } while (cc1 != ccStart);


                // Do not create faces if the internal point is outside!
                // This occurs because the internal point is not determined to
                // be outside in the inside/outside test. This is most likely
                // due to the triangle.nearestPointClassify test not returning
                // edge/point as the nearest type.

                if (skipEdge)
                {
                    continue;
                }
            }

            face newDualFace = buildDualFace(eit);

            if (newDualFace.size() >= 3)
            {
                label own = -1;
                label nei = -1;

                if (ownerAndNeighbour(vA, vB, own, nei))
                {
                    reverse(newDualFace);
                }

                label patchIndex = -1;

                pointFromPoint ptA = topoint(vA->point());
                pointFromPoint ptB = topoint(vB->point());

                if (nei == -1)
                {
                    // boundary face

                    if (isProcBoundaryEdge(eit))
                    {
                        // One (and only one) of the points is an internal
                        // point from another processor

                        label procIndex = max(vA->procIndex(), vB->procIndex());

                        patchIndex = max
                        (
                            findIndex(procNeighbours, vA->procIndex()),
                            findIndex(procNeighbours, vB->procIndex())
                        );

                        // The lower processor index is the owner of the
                        // two for the purpose of sorting the patch faces.

                        if (Pstream::myProcNo() < procIndex)
                        {
                            // Use this processor's vertex index as the master
                            // for sorting

                            DynamicList<Pair<labelPair>>& sortingIndex =
                                procPatchSortingIndex[patchIndex];

                            if (vB->internalOrBoundaryPoint() && vB->referred())
                            {
                                sortingIndex.append
                                (
                                    Pair<labelPair>
                                    (
                                        labelPair(vA->index(), vA->procIndex()),
                                        labelPair(vB->index(), vB->procIndex())
                                    )
                                );
                            }
                            else
                            {
                                sortingIndex.append
                                (
                                    Pair<labelPair>
                                    (
                                        labelPair(vB->index(), vB->procIndex()),
                                        labelPair(vA->index(), vA->procIndex())
                                    )
                                );
                            }
                        }
                        else
                        {
                            // Use the other processor's vertex index as the
                            // master for sorting

                            DynamicList<Pair<labelPair>>& sortingIndex =
                                procPatchSortingIndex[patchIndex];

                            if (vA->internalOrBoundaryPoint() && vA->referred())
                            {
                                sortingIndex.append
                                (
                                    Pair<labelPair>
                                    (
                                        labelPair(vA->index(), vA->procIndex()),
                                        labelPair(vB->index(), vB->procIndex())
                                    )
                                );
                            }
                            else
                            {
                                sortingIndex.append
                                (
                                    Pair<labelPair>
                                    (
                                        labelPair(vB->index(), vB->procIndex()),
                                        labelPair(vA->index(), vA->procIndex())
                                    )
                                );
                            }
                        }

//                        Pout<< ptA << " " << ptB
//                            << " proc indices "
//                            << vA->procIndex() << " " << vB->procIndex()
//                            << " indices " << vA->index()
//                            << " " << vB->index()
//                            << " my proc " << Pstream::myProcNo()
//                            << " addedIndex "
//                            << procPatchSortingIndex[patchIndex].last()
//                            << endl;
                    }
                    else
                    {
                        patchIndex = geometryToConformTo_.findPatch(ptA, ptB);
                    }

                    if (patchIndex == -1)
                    {
                        // Did not find a surface patch between
                        // between Dv pair, finding nearest patch

//                         Pout<< "Did not find a surface patch between "
//                             << "for face, finding nearest patch to"
//                             << 0.5*(ptA + ptB) << endl;

                        patchIndex = geometryToConformTo_.findPatch
                        (
                            0.5*(ptA + ptB)
                        );
                    }

                    patchFaces[patchIndex].append(newDualFace);
                    patchOwners[patchIndex].append(own);

                    // If the two vertices are a pair, then the patch face is
                    // a desired one.
                    if
                    (
                        vA->boundaryPoint() && vB->boundaryPoint()
                     && !ptPairs_.isPointPair(vA, vB)
                     && !ftPtConformer_.featurePointPairs().isPointPair(vA, vB)
                    )
                    {
                        indirectPatchFace[patchIndex].append(true);
                    }
                    else
                    {
                        indirectPatchFace[patchIndex].append(false);
                    }

                    // Store the non-internal or boundary point
                    if (vA->internalOrBoundaryPoint())
                    {
                        patchPPSlaves[patchIndex].append(vB->index());
                    }
                    else
                    {
                        patchPPSlaves[patchIndex].append(vA->index());
                    }
                }
                else
                {
                    if
                    (
                        !vA->boundaryPoint()
                     || !vB->boundaryPoint()
                     || ptPairs_.isPointPair(vA, vB)
                     || ftPtConformer_.featurePointPairs().isPointPair(vA, vB)
                    )
                    {
                        patchIndex = geometryToConformTo_.findPatch(ptA, ptB);
                    }

                    if
                    (
                        patchIndex != -1
                     && geometryToConformTo_.patchInfo().set(patchIndex)
                    )
                    {
                        // baffle faces

                        patchFaces[patchIndex].append(newDualFace);
                        patchOwners[patchIndex].append(own);
                        indirectPatchFace[patchIndex].append(false);

                        reverse(newDualFace);

                        patchFaces[patchIndex].append(newDualFace);
                        patchOwners[patchIndex].append(nei);
                        indirectPatchFace[patchIndex].append(false);

                        if
                        (
                            labelPair(vB->index(), vB->procIndex())
                          < labelPair(vA->index(), vA->procIndex())
                        )
                        {
                            patchPPSlaves[patchIndex].append(vB->index());
                            patchPPSlaves[patchIndex].append(vB->index());
                        }
                        else
                        {
                            patchPPSlaves[patchIndex].append(vA->index());
                            patchPPSlaves[patchIndex].append(vA->index());
                        }

                    }
                    else
                    {
                        // internal face
                        faces[dualFacei] = newDualFace;
                        owner[dualFacei] = own;
                        neighbour[dualFacei] = nei;

                        dualFacei++;
                    }
                }
            }
        }
    }

    if (!patchFaces[defaultPatchIndex].empty())
    {
        Pout<< nl << patchFaces[defaultPatchIndex].size()
            << " faces were not able to have their patch determined from "
            << "the surface. "
            << nl <<  "Adding to patch " << patchNames[defaultPatchIndex]
            << endl;
    }

    label nInternalFaces = dualFacei;

    faces.setSize(nInternalFaces);
    owner.setSize(nInternalFaces);
    neighbour.setSize(nInternalFaces);

    timeCheck("polyMesh quality checked");

    sortFaces(faces, owner, neighbour);

    sortProcPatches
    (
        patchFaces,
        patchOwners,
        patchPPSlaves,
        procPatchSortingIndex
    );

    timeCheck("faces, owner, neighbour sorted");

    addPatches
    (
        nInternalFaces,
        faces,
        owner,
        patchDicts,
        boundaryFacesToRemove,
        patchFaces,
        patchOwners,
        indirectPatchFace
    );

    // Return     patchPointPairSlaves.setSize(nPatches);
    patchPointPairSlaves.setSize(nPatches);
    forAll(patchPPSlaves, patchi)
    {
        patchPointPairSlaves[patchi].transfer(patchPPSlaves[patchi]);
    }

    if (foamyHexMeshControls().objOutput())
    {
        Info<< "Writing processor interfaces" << endl;

        forAll(patchDicts, nbI)
        {
            if (patchFaces[nbI].size() > 0)
            {
                const label neighbour =
                (
                    patchDicts[nbI].found("neighbProcNo")
                  ? readLabel(patchDicts[nbI].lookup("neighbProcNo"))
                  : -1
                );

                faceList procPatchFaces = patchFaces[nbI];

                // Reverse faces as it makes it easier to analyse the output
                // using a diff
                if (neighbour < Pstream::myProcNo())
                {
                    forAll(procPatchFaces, fI)
                    {
                        procPatchFaces[fI] = procPatchFaces[fI].reverseFace();
                    }
                }

                if (neighbour != -1)
                {
                    word fName =
                        "processor_"
                      + name(Pstream::myProcNo())
                      + "_to_"
                      + name(neighbour)
                      + "_interface.obj";

                    DelaunayMeshTools::writeProcessorInterface
                    (
                        time().path()/fName,
                        *this,
                        procPatchFaces
                    );
                }
            }
        }
    }
}


void Foam::conformalVoronoiMesh::sortFaces
(
    faceList& faces,
    labelList& owner,
    labelList& neighbour
) const
{
    // Upper triangular order:
    // + owner is sorted in ascending cell order
    // + within each block of equal value for owner, neighbour is sorted in
    //   ascending cell order.
    // + faces sorted to correspond
    // e.g.
    // owner | neighbour
    // 0     | 2
    // 0     | 23
    // 0     | 71
    // 1     | 23
    // 1     | 24
    // 1     | 91

    List<labelPair> ownerNeighbourPair(owner.size());

    forAll(ownerNeighbourPair, oNI)
    {
        ownerNeighbourPair[oNI] = labelPair(owner[oNI], neighbour[oNI]);
    }

    Info<< nl
        << "Sorting faces, owner and neighbour into upper triangular order"
        << endl;

    labelList oldToNew;

    sortedOrder(ownerNeighbourPair, oldToNew);

    oldToNew = invert(oldToNew.size(), oldToNew);

    inplaceReorder(oldToNew, faces);
    inplaceReorder(oldToNew, owner);
    inplaceReorder(oldToNew, neighbour);
}


void Foam::conformalVoronoiMesh::sortProcPatches
(
    List<DynamicList<face>>& patchFaces,
    List<DynamicList<label>>& patchOwners,
    List<DynamicList<label>>& patchPointPairSlaves,
    labelPairPairDynListList& patchSortingIndices
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    forAll(patchSortingIndices, patchi)
    {
        faceList& faces = patchFaces[patchi];
        labelList& owner = patchOwners[patchi];
        DynamicList<label>& slaves = patchPointPairSlaves[patchi];
        DynamicList<Pair<labelPair>>& sortingIndices
            = patchSortingIndices[patchi];

        if (!sortingIndices.empty())
        {
            if
            (
                faces.size() != sortingIndices.size()
             || owner.size() != sortingIndices.size()
             || slaves.size() != sortingIndices.size()
            )
            {
                FatalErrorInFunction
                    << "patch size and size of sorting indices is inconsistent "
                    << " for patch " << patchi << nl
                    << " faces.size() " << faces.size() << nl
                    << " owner.size() " << owner.size() << nl
                    << " slaves.size() " << slaves.size() << nl
                    << " sortingIndices.size() "
                    << sortingIndices.size()
                    << exit(FatalError) << endl;
            }

            labelList oldToNew;

            sortedOrder(sortingIndices, oldToNew);

            oldToNew = invert(oldToNew.size(), oldToNew);

            inplaceReorder(oldToNew, sortingIndices);
            inplaceReorder(oldToNew, faces);
            inplaceReorder(oldToNew, owner);
            inplaceReorder(oldToNew, slaves);
        }
    }
}


void Foam::conformalVoronoiMesh::addPatches
(
    const label nInternalFaces,
    faceList& faces,
    labelList& owner,
    PtrList<dictionary>& patchDicts,
    PackedBoolList& boundaryFacesToRemove,
    const List<DynamicList<face>>& patchFaces,
    const List<DynamicList<label>>& patchOwners,
    const List<DynamicList<bool>>& indirectPatchFace
) const
{
    label nBoundaryFaces = 0;

    forAll(patchFaces, p)
    {
        patchDicts[p].set("nFaces", patchFaces[p].size());
        patchDicts[p].set("startFace", nInternalFaces + nBoundaryFaces);

        nBoundaryFaces += patchFaces[p].size();
    }

    faces.setSize(nInternalFaces + nBoundaryFaces);
    owner.setSize(nInternalFaces + nBoundaryFaces);
    boundaryFacesToRemove.setSize(nInternalFaces + nBoundaryFaces);

    label facei = nInternalFaces;

    forAll(patchFaces, p)
    {
        forAll(patchFaces[p], f)
        {
            faces[facei] = patchFaces[p][f];
            owner[facei] = patchOwners[p][f];
            boundaryFacesToRemove[facei] = indirectPatchFace[p][f];

            facei++;
        }
    }
}


void Foam::conformalVoronoiMesh::removeUnusedPoints
(
    faceList& faces,
    pointField& pts,
    labelList& boundaryPts
) const
{
    Info<< nl << "Removing unused points" << endl;

    PackedBoolList ptUsed(pts.size(), false);

    // Scan all faces to find all of the points that are used

    forAll(faces, fI)
    {
        const face& f = faces[fI];

        forAll(f, fPtI)
        {
            ptUsed[f[fPtI]] = true;
        }
    }

    label pointi = 0;

    labelList oldToNew(pts.size(), label(-1));

    // Move all of the used points to the start of the pointField and
    // truncate it

    forAll(ptUsed, ptUI)
    {
        if (ptUsed[ptUI] == true)
        {
            oldToNew[ptUI] = pointi++;
        }
    }

    inplaceReorder(oldToNew, pts);
    inplaceReorder(oldToNew, boundaryPts);

    Info<< "    Removing "
        << returnReduce(pts.size() - pointi, sumOp<label>())
        << " unused points"
        << endl;

    pts.setSize(pointi);
    boundaryPts.setSize(pointi);

    // Renumber the faces to use the new point numbers

    forAll(faces, fI)
    {
        inplaceRenumber(oldToNew, faces[fI]);
    }
}


Foam::labelList Foam::conformalVoronoiMesh::removeUnusedCells
(
    labelList& owner,
    labelList& neighbour
) const
{
    Info<< nl << "Removing unused cells" << endl;

    PackedBoolList cellUsed(vertexCount(), false);

    // Scan all faces to find all of the cells that are used

    forAll(owner, oI)
    {
        cellUsed[owner[oI]] = true;
    }

    forAll(neighbour, nI)
    {
        cellUsed[neighbour[nI]] = true;
    }

    label celli = 0;

    labelList oldToNew(cellUsed.size(), label(-1));

    // Move all of the used cellCentres to the start of the pointField and
    // truncate it

    forAll(cellUsed, cellUI)
    {
        if (cellUsed[cellUI] == true)
        {
            oldToNew[cellUI] = celli++;
        }
    }

    labelList newToOld(invert(celli, oldToNew));

    // Find all of the unused cells, create a list of them, then
    // subtract one from each owner and neighbour entry for each of
    // the unused cell indices that it is above.

    DynamicList<label> unusedCells;

    forAll(cellUsed, cUI)
    {
        if (cellUsed[cUI] == false)
        {
            unusedCells.append(cUI);
        }
    }

    if (unusedCells.size() > 0)
    {
        Info<< "    Removing "
            << returnReduce(unusedCells.size(), sumOp<label>())
            <<  " unused cell labels" << endl;

        forAll(owner, oI)
        {
            label& o = owner[oI];

            o -= findLower(unusedCells, o) + 1;
        }

        forAll(neighbour, nI)
        {
            label& n = neighbour[nI];

            n -= findLower(unusedCells, n) + 1;
        }
    }

    return newToOld;
}


// ************************************************************************* //
