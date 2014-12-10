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

\*---------------------------------------------------------------------------*/

#include "edgeCollapser.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "globalMeshData.H"
#include "syncTools.H"
#include "PointEdgeWave.H"
#include "globalIndex.H"
#include "removePoints.H"
#include "motionSmoother.H"
#include "OFstream.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(edgeCollapser, 0);
}


Foam::HashSet<Foam::label> Foam::edgeCollapser::checkBadFaces
(
    const polyMesh& mesh,
    const dictionary& meshQualityDict
)
{
    labelHashSet badFaces(mesh.nFaces()/100);
    DynamicList<label> checkFaces(mesh.nFaces());

    const vectorField& fAreas = mesh.faceAreas();

    scalar faceAreaLimit = SMALL;

    forAll(fAreas, fI)
    {
        if (mag(fAreas[fI]) > faceAreaLimit)
        {
            checkFaces.append(fI);
        }
    }

    Info<< endl;

    motionSmoother::checkMesh
    (
        false,
        mesh,
        meshQualityDict,
        checkFaces,
        badFaces
    );

    return badFaces;
}


Foam::label Foam::edgeCollapser::checkMeshQuality
(
    const polyMesh& mesh,
    const dictionary& meshQualityDict,
    PackedBoolList& isErrorPoint
)
{
    labelHashSet badFaces = edgeCollapser::checkBadFaces
    (
        mesh,
        meshQualityDict
    );

    label nBadFaces = returnReduce(badFaces.size(), sumOp<label>());

    forAllConstIter(labelHashSet, badFaces, iter)
    {
        const face& f = mesh.faces()[iter.key()];

        forAll(f, pI)
        {
            isErrorPoint[f[pI]] = true;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        isErrorPoint,
        orEqOp<unsigned int>(),
        0
    );

    return nBadFaces;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::edgeCollapser::edgesFromPoints
(
    const label& faceI,
    const labelList& pointLabels
) const
{
    labelList edgeLabels(pointLabels.size() - 1, -1);

    const labelList& faceEdges = mesh_.faceEdges()[faceI];
    const edgeList& edges = mesh_.edges();

    label count = 0;

    forAll(faceEdges, eI)
    {
        const label edgeI = faceEdges[eI];
        const edge& e = edges[edgeI];

        label match = 0;

        forAll(pointLabels, pI)
        {
            if (e[0] == pointLabels[pI])
            {
                match++;
            }

            if (e[1] == pointLabels[pI])
            {
                match++;
            }
        }

        if (match == 2)
        {
            // Edge found
            edgeLabels[count++] = edgeI;
        }
    }

    if (count != edgeLabels.size())
    {
        edgeLabels.setSize(count);
    }

    return edgeLabels;
}


void Foam::edgeCollapser::collapseToEdge
(
    const label faceI,
    const pointField& pts,
    const labelList& pointPriority,
    const vector& collapseAxis,
    const point& fC,
    const labelList& facePtsNeg,
    const labelList& facePtsPos,
    const scalarList& dNeg,
    const scalarList& dPos,
    const scalar dShift,
    PackedBoolList& collapseEdge,
    Map<point>& collapsePointToLocation
) const
{
    // Negative half

    Foam::point collapseToPtA(GREAT, GREAT, GREAT);
        //collapseAxis*(sum(dNeg)/dNeg.size() - dShift) + fC;

    label maxPriority = labelMin;
    DynamicList<label> maxPriorityPts(max(dNeg.size(), dPos.size()));

    forAll(facePtsNeg, fPtI)
    {
        const label facePointI = facePtsNeg[fPtI];
        const label facePtPriority = pointPriority[facePointI];

        if (facePtPriority > maxPriority)
        {
            maxPriority = facePtPriority;
            maxPriorityPts.clear();
            maxPriorityPts.append(facePointI);
        }
        else if (facePtPriority == maxPriority)
        {
            maxPriorityPts.append(facePointI);
        }
    }

    if (!maxPriorityPts.empty())
    {
        Foam::point averagePt(vector::zero);

        forAll(maxPriorityPts, ptI)
        {
            averagePt += pts[maxPriorityPts[ptI]];
        }

        collapseToPtA = averagePt/maxPriorityPts.size();
//        collapseToPtA = pts[maxPriorityPts.first()];
    }

    maxPriority = labelMin;
    maxPriorityPts.clear();

    labelList faceEdgesNeg = edgesFromPoints(faceI, facePtsNeg);

    forAll(faceEdgesNeg, edgeI)
    {
        collapseEdge[faceEdgesNeg[edgeI]] = true;
    }

    forAll(facePtsNeg, pI)
    {
        collapsePointToLocation.set(facePtsNeg[pI], collapseToPtA);
    }


    // Positive half
    Foam::point collapseToPtB(GREAT, GREAT, GREAT);
//        = collapseAxis*(sum(dPos)/dPos.size() - dShift) + fC;

    forAll(facePtsPos, fPtI)
    {
        const label facePointI = facePtsPos[fPtI];
        const label facePtPriority = pointPriority[facePointI];

        if (facePtPriority > maxPriority)
        {
            maxPriority = facePtPriority;
            maxPriorityPts.clear();
            maxPriorityPts.append(facePointI);
        }
        else if (facePtPriority == maxPriority)
        {
            maxPriorityPts.append(facePointI);
        }
    }

    if (!maxPriorityPts.empty())
    {
        Foam::point averagePt(vector::zero);

        forAll(maxPriorityPts, ptI)
        {
            averagePt += pts[maxPriorityPts[ptI]];
        }

        collapseToPtB = averagePt/maxPriorityPts.size();
//        collapseToPtB = pts[maxPriorityPts.first()];
    }

    labelList faceEdgesPos = edgesFromPoints(faceI, facePtsPos);

    forAll(faceEdgesPos, edgeI)
    {
        collapseEdge[faceEdgesPos[edgeI]] = true;
    }

    forAll(facePtsPos, pI)
    {
        collapsePointToLocation.set(facePtsPos[pI], collapseToPtB);
    }
}


void Foam::edgeCollapser::collapseToPoint
(
    const label& faceI,
    const pointField& pts,
    const labelList& pointPriority,
    const point& fC,
    const labelList& facePts,
    PackedBoolList& collapseEdge,
    Map<point>& collapsePointToLocation
) const
{
    const face& f = mesh_.faces()[faceI];

    Foam::point collapseToPt = fC;

    label maxPriority = labelMin;
    DynamicList<label> maxPriorityPts(f.size());

    forAll(facePts, fPtI)
    {
        const label facePointI = facePts[fPtI];
        const label facePtPriority = pointPriority[facePointI];

        if (facePtPriority > maxPriority)
        {
            maxPriority = facePtPriority;
            maxPriorityPts.clear();
            maxPriorityPts.append(facePointI);
        }
        else if (facePtPriority == maxPriority)
        {
            maxPriorityPts.append(facePointI);
        }
    }

    if (!maxPriorityPts.empty())
    {
        Foam::point averagePt(vector::zero);

        forAll(maxPriorityPts, ptI)
        {
            averagePt += pts[maxPriorityPts[ptI]];
        }

        collapseToPt = averagePt/maxPriorityPts.size();

//        collapseToPt = pts[maxPriorityPts.first()];
    }

//    DynamicList<label> faceBoundaryPts(f.size());
//    DynamicList<label> faceFeaturePts(f.size());
//
//    forAll(facePts, fPtI)
//    {
//        if (pointPriority[facePts[fPtI]] == 1)
//        {
//            faceFeaturePts.append(facePts[fPtI]);
//        }
//        else if (pointPriority[facePts[fPtI]] == 0)
//        {
//            faceBoundaryPts.append(facePts[fPtI]);
//        }
//    }
//
//    if (!faceBoundaryPts.empty() || !faceFeaturePts.empty())
//    {
//        if (!faceFeaturePts.empty())
//        {
//            collapseToPt = pts[faceFeaturePts.first()];
//        }
//        else if (faceBoundaryPts.size() == 2)
//        {
//            collapseToPt =
//                0.5
//               *(
//                    pts[faceBoundaryPts[0]]
//                  + pts[faceBoundaryPts[1]]
//                );
//        }
//        else if (faceBoundaryPts.size() <= f.size())
//        {
//            face bFace(faceBoundaryPts);
//
//            collapseToPt = bFace.centre(pts);
//        }
//    }

    const labelList& faceEdges = mesh_.faceEdges()[faceI];

    forAll(faceEdges, eI)
    {
        const label edgeI = faceEdges[eI];
        collapseEdge[edgeI] = true;
    }

    forAll(f, pI)
    {
        collapsePointToLocation.set(f[pI], collapseToPt);
    }
}


void Foam::edgeCollapser::faceCollapseAxisAndAspectRatio
(
    const face& f,
    const point& fC,
    vector& collapseAxis,
    scalar& aspectRatio
) const
{
    const pointField& pts = mesh_.points();

    tensor J = f.inertia(pts, fC);

    // Find the dominant collapse direction by finding the eigenvector
    // that corresponds to the normal direction, discarding it.  The
    // eigenvector corresponding to the smaller of the two remaining
    // eigenvalues is the dominant axis in a high aspect ratio face.

    scalar magJ = mag(J);

    scalar detJ = SMALL;

    if (magJ > VSMALL)
    {
        // Normalise inertia tensor to remove problems with small values

        J /= mag(J);
        // J /= cmptMax(J);
        // J /= max(eigenValues(J).x(), SMALL);

        // Calculating determinant, including stabilisation for zero or
        // small negative values

        detJ = max(det(J), SMALL);
    }

    if (detJ < 1e-5)
    {
        collapseAxis = f.edges()[longestEdge(f, pts)].vec(pts);

        // It is possible that all the points of a face are the same
        if (magSqr(collapseAxis) > VSMALL)
        {
            collapseAxis /= mag(collapseAxis);
        }

        // Empirical correlation for high aspect ratio faces

        aspectRatio = Foam::sqrt(0.35/detJ);
    }
    else
    {
        vector eVals = eigenValues(J);

        if (mag(eVals.y() - eVals.x()) < 100*SMALL)
        {
            // First two eigenvalues are the same: i.e. a square face

            // Cannot necessarily determine linearly independent
            // eigenvectors, or any at all, use longest edge direction.

            collapseAxis = f.edges()[longestEdge(f, pts)].vec(pts);

            collapseAxis /= mag(collapseAxis);

            aspectRatio = 1.0;
        }
        else
        {
            // The maximum eigenvalue (z()) must be the direction of the
            // normal, as it has the greatest value.  The minimum eigenvalue
            // is the dominant collapse axis for high aspect ratio faces.

            collapseAxis = eigenVector(J, eVals.x());

            // The inertia calculation describes the mass distribution as a
            // function of distance squared to the axis, so the square root of
            // the ratio of face-plane moments gives a good indication of the
            // aspect ratio.

            aspectRatio = Foam::sqrt(eVals.y()/max(eVals.x(), SMALL));
        }
    }
}


Foam::scalarField Foam::edgeCollapser::calcTargetFaceSizes() const
{
    scalarField targetFaceSizes(mesh_.nFaces(), -1);

    const scalarField& V = mesh_.cellVolumes();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const labelList& cellOwner = mesh_.faceOwner();
    const labelList& cellNeighbour = mesh_.faceNeighbour();

    const label nBoundaryFaces = mesh_.nFaces() - mesh_.nInternalFaces();

    // Calculate face size from cell volumes for internal faces
    for (label intFaceI = 0; intFaceI < mesh_.nInternalFaces(); ++intFaceI)
    {
        const scalar cellOwnerVol = max(0.0, V[cellOwner[intFaceI]]);
        const scalar cellNeighbourVol = max(0.0, V[cellNeighbour[intFaceI]]);

        scalar targetFaceSizeA = Foam::pow(cellOwnerVol, 1.0/3.0);
        scalar targetFaceSizeB = Foam::pow(cellNeighbourVol, 1.0/3.0);

        targetFaceSizes[intFaceI] = 0.5*(targetFaceSizeA + targetFaceSizeB);
    }

    scalarField neiCellVolumes(nBoundaryFaces, -1);

    // Now do boundary faces
    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        label bFaceI = patch.start() - mesh_.nInternalFaces();

        if (patch.coupled())
        {
            // Processor boundary face: Need to get the cell volume on the other
            // processor
            const labelUList& faceCells = patch.faceCells();

            forAll(faceCells, facei)
            {
                neiCellVolumes[bFaceI++] = max(0.0, V[faceCells[facei]]);
            }
        }
        else
        {
            // Normal boundary face: Just use owner cell volume to calculate
            // the target face size
            forAll(patch, patchFaceI)
            {
                const label extFaceI = patchFaceI + patch.start();
                const scalar cellOwnerVol = max(0.0, V[cellOwner[extFaceI]]);

                targetFaceSizes[extFaceI] = Foam::pow(cellOwnerVol, 1.0/3.0);
            }
        }
    }

    syncTools::swapBoundaryFaceList(mesh_, neiCellVolumes);

    forAll(patches, patchI)
    {
        const polyPatch& patch = patches[patchI];

        label bFaceI = patch.start() - mesh_.nInternalFaces();

        if (patch.coupled())
        {
            forAll(patch, patchFaceI)
            {
                const label localFaceI = patchFaceI + patch.start();
                const scalar cellOwnerVol = max(0.0, V[cellOwner[localFaceI]]);
                const scalar cellNeighbourVol = neiCellVolumes[bFaceI++];

                scalar targetFaceSizeA = Foam::pow(cellOwnerVol, 1.0/3.0);
                scalar targetFaceSizeB = Foam::pow(cellNeighbourVol, 1.0/3.0);

                targetFaceSizes[localFaceI]
                    = 0.5*(targetFaceSizeA + targetFaceSizeB);
            }
        }
    }

    // Returns a characteristic length, not an area
    return targetFaceSizes;
}


Foam::edgeCollapser::collapseType Foam::edgeCollapser::collapseFace
(
    const labelList& pointPriority,
    const face& f,
    const label faceI,
    const scalar targetFaceSize,
    PackedBoolList& collapseEdge,
    Map<point>& collapsePointToLocation,
    const scalarField& faceFilterFactor
) const
{
    const scalar collapseSizeLimitCoeff = faceFilterFactor[faceI];

    const pointField& pts = mesh_.points();

    labelList facePts(f);

    const Foam::point fC = f.centre(pts);

    const scalar fA = f.mag(pts);

    vector collapseAxis = vector::zero;
    scalar aspectRatio = 1.0;

    faceCollapseAxisAndAspectRatio(f, fC, collapseAxis, aspectRatio);

    // The signed distance along the collapse axis passing through the
    // face centre that each vertex projects to.

    scalarField d(f.size());

    forAll(f, fPtI)
    {
        const Foam::point& pt = pts[f[fPtI]];

        d[fPtI] = (collapseAxis & (pt - fC));
    }

    // Sort the projected distances and the corresponding vertex
    // indices along the collapse axis

    labelList oldToNew;

    sortedOrder(d, oldToNew);

    oldToNew = invert(oldToNew.size(), oldToNew);

    inplaceReorder(oldToNew, d);

    inplaceReorder(oldToNew, facePts);

    // Shift the points so that they are relative to the centre of the
    // collapse line.

    scalar dShift = -0.5*(d.first() + d.last());

    d += dShift;

    // Form two lists, one for each half of the set of points
    // projected along the collapse axis.

    // Middle value, index of first entry in the second half
    label middle = -1;

    forAll(d, dI)
    {
        if (d[dI] > 0)
        {
            middle = dI;

            break;
        }
    }

    if (middle == -1)
    {
//        SeriousErrorIn("collapseFace")
//            << "middle == -1, " << f << " " << d
//            << endl;//abort(FatalError);

        return noCollapse;
    }

    // Negative half
    SubList<scalar> dNeg(d, middle, 0);
    SubList<label> facePtsNeg(facePts, middle, 0);

    // Positive half
    SubList<scalar> dPos(d, d.size() - middle, middle);
    SubList<label> facePtsPos(facePts, d.size() - middle, middle);

    // Defining how close to the midpoint (M) of the projected
    // vertices line a projected vertex (X) can be before making this
    // an invalid edge collapse
    //
    // X---X-g----------------M----X-----------g----X--X
    //
    // Only allow a collapse if all projected vertices are outwith
    // guardFraction (g) of the distance form the face centre to the
    // furthest vertex in the considered direction

    if (dNeg.size() == 0 || dPos.size() == 0)
    {
        WarningIn
        (
            "Foam::conformalVoronoiMesh::collapseFace"
        )
            << "All points on one side of face centre, not collapsing."
            << endl;
    }

//    Info<< "Face : " << f << nl
//        << "    Collapse Axis: " << collapseAxis << nl
//        << "    Aspect Ratio : " << aspectRatio << endl;

    collapseType typeOfCollapse = noCollapse;

    if (magSqr(collapseAxis) < VSMALL)
    {
        typeOfCollapse = toPoint;
    }
    else if (fA < aspectRatio*sqr(targetFaceSize*collapseSizeLimitCoeff))
    {
        if
        (
            allowEarlyCollapseToPoint_
         && (d.last() - d.first())
          < targetFaceSize
           *allowEarlyCollapseCoeff_*maxCollapseFaceToPointSideLengthCoeff_
        )
        {
            typeOfCollapse = toPoint;
        }
        else if
        (
            (dNeg.last() < guardFraction_*dNeg.first())
         && (dPos.first() > guardFraction_*dPos.last())
        )
        {
            typeOfCollapse = toEdge;
        }
        else if
        (
            (d.last() - d.first())
          < targetFaceSize
           *maxCollapseFaceToPointSideLengthCoeff_
        )
        {
            // If the face can't be collapsed to an edge, and it has a
            // small enough span, collapse it to a point.
            typeOfCollapse = toPoint;
        }
    }

    if (typeOfCollapse == toPoint)
    {
        collapseToPoint
        (
            faceI,
            pts,
            pointPriority,
            fC,
            facePts,
            collapseEdge,
            collapsePointToLocation
        );
    }
    else if (typeOfCollapse == toEdge)
    {
        collapseToEdge
        (
            faceI,
            pts,
            pointPriority,
            collapseAxis,
            fC,
            facePtsNeg,
            facePtsPos,
            dNeg,
            dPos,
            dShift,
            collapseEdge,
            collapsePointToLocation
        );
    }

    return typeOfCollapse;
}


Foam::label Foam::edgeCollapser::edgeMaster
(
    const labelList& pointPriority,
    const edge& e
) const
{
    label masterPoint = -1;

    const label e0 = e.start();
    const label e1 = e.end();

    const label e0Priority = pointPriority[e0];
    const label e1Priority = pointPriority[e1];

    if (e0Priority > e1Priority)
    {
        masterPoint = e0;
    }
    else if (e0Priority < e1Priority)
    {
        masterPoint = e1;
    }
    else if (e0Priority == e1Priority)
    {
        masterPoint = e0;
    }

//    // Collapse edge to point with higher priority.
//    if (pointPriority[e0] >= 0)
//    {
//        if (pointPriority[e1] >= 0)
//        {
//            // Both points have high priority. Choose one to collapse to.
//            // Note: should look at feature edges/points!
//            masterPoint = e0;
//        }
//        else
//        {
//            masterPoint = e0;
//        }
//    }
//    else
//    {
//        if (pointPriority[e1] >= 0)
//        {
//            masterPoint = e1;
//        }
//        else
//        {
//            // None on boundary. Neither is a master.
//            return -1;
//        }
//    }

    return masterPoint;
}


void Foam::edgeCollapser::checkBoundaryPointMergeEdges
(
    const label pointI,
    const label otherPointI,
    const labelList& pointPriority,
    Map<point>& collapsePointToLocation
) const
{
   const pointField& points = mesh_.points();

   const label e0Priority = pointPriority[pointI];
   const label e1Priority = pointPriority[otherPointI];

   if (e0Priority > e1Priority)
   {
       collapsePointToLocation.set
       (
           otherPointI,
           points[pointI]
       );
   }
   else if (e0Priority < e1Priority)
   {
       collapsePointToLocation.set
       (
           pointI,
           points[otherPointI]
       );
   }
   else // e0Priority == e1Priority
   {
       collapsePointToLocation.set
       (
           pointI,
           points[otherPointI]
       );

//       Foam::point averagePt
//       (
//           0.5*(points[otherPointI] + points[pointI])
//       );
//
//       collapsePointToLocation.set(pointI, averagePt);
//       collapsePointToLocation.set(otherPointI, averagePt);
   }
}


Foam::label Foam::edgeCollapser::breakStringsAtEdges
(
    const PackedBoolList& markedEdges,
    PackedBoolList& collapseEdge,
    List<pointEdgeCollapse>& allPointInfo
) const
{
    const edgeList& edges = mesh_.edges();
    const labelListList& pointEdges = mesh_.pointEdges();

    label nUncollapsed = 0;

    forAll(edges, eI)
    {
        if (markedEdges[eI])
        {
            const edge& e = edges[eI];

            const label startCollapseIndex
                = allPointInfo[e.start()].collapseIndex();

            if (startCollapseIndex != -1 && startCollapseIndex != -2)
            {
                const label endCollapseIndex
                    = allPointInfo[e.end()].collapseIndex();

                if
                (
                    !collapseEdge[eI]
                 && startCollapseIndex == endCollapseIndex
                )
                {
                    const labelList& ptEdgesStart = pointEdges[e.start()];

                    forAll(ptEdgesStart, ptEdgeI)
                    {
                        const label edgeI = ptEdgesStart[ptEdgeI];

                        const label nbrPointI
                            = edges[edgeI].otherVertex(e.start());
                        const label nbrIndex
                            = allPointInfo[nbrPointI].collapseIndex();

                        if
                        (
                            collapseEdge[edgeI]
                         && nbrIndex == startCollapseIndex
                        )
                        {
                            collapseEdge[edgeI] = false;
                            nUncollapsed++;
                        }
                    }
                }
            }
        }
    }

    return nUncollapsed;
}


void Foam::edgeCollapser::determineDuplicatePointsOnFace
(
    const face& f,
    PackedBoolList& markedPoints,
    labelHashSet& uniqueCollapses,
    labelHashSet& duplicateCollapses,
    List<pointEdgeCollapse>& allPointInfo
) const
{
    uniqueCollapses.clear();
    duplicateCollapses.clear();

    forAll(f, fpI)
    {
        label index = allPointInfo[f[fpI]].collapseIndex();

        // Check for consecutive duplicate
        if (index != allPointInfo[f.prevLabel(fpI)].collapseIndex())
        {
            if (!uniqueCollapses.insert(index))
            {
                // Failed inserting so duplicate
                duplicateCollapses.insert(index);
            }
        }
    }

    // Now duplicateCollapses contains duplicate collapse indices.
    // Convert to points.
    forAll(f, fpI)
    {
        label index = allPointInfo[f[fpI]].collapseIndex();
        if (duplicateCollapses.found(index))
        {
            markedPoints[f[fpI]] = true;
        }
    }
}


Foam::label Foam::edgeCollapser::countEdgesOnFace
(
    const face& f,
    List<pointEdgeCollapse>& allPointInfo
) const
{
    label nEdges = 0;

    forAll(f, fpI)
    {
        const label pointI = f[fpI];
        const label newPointI = allPointInfo[pointI].collapseIndex();

        if (newPointI == -2)
        {
            nEdges++;
        }
        else
        {
            const label prevPointI = f[f.fcIndex(fpI)];
            const label prevNewPointI
                = allPointInfo[prevPointI].collapseIndex();

            if (newPointI != prevNewPointI)
            {
                nEdges++;
            }
        }
    }

    return nEdges;
}


bool Foam::edgeCollapser::isFaceCollapsed
(
    const face& f,
    List<pointEdgeCollapse>& allPointInfo
) const
{
    label nEdges = countEdgesOnFace(f, allPointInfo);

    // Polygons must have 3 or more edges to be valid
    if (nEdges < 3)
    {
        return true;
    }

    return false;
}


// Create consistent set of collapses.
//  collapseEdge : per edge:
//      -1 : do not collapse
//       0 : collapse to start
//       1 : collapse to end
//  Note: collapseEdge has to be parallel consistent (in orientation)
Foam::label Foam::edgeCollapser::syncCollapse
(
    const globalIndex& globalPoints,
    const labelList& pointPriority,
    const PackedBoolList& collapseEdge,
    const Map<point>& collapsePointToLocation,
    List<pointEdgeCollapse>& allPointInfo
) const
{
    const edgeList& edges = mesh_.edges();

    label nCollapsed = 0;

    DynamicList<label> initPoints(mesh_.nPoints());
    DynamicList<pointEdgeCollapse> initPointInfo(mesh_.nPoints());

    allPointInfo.clear();
    allPointInfo.setSize(mesh_.nPoints());

    // Initialise edges to no collapse
    List<pointEdgeCollapse> allEdgeInfo
    (
        mesh_.nEdges(),
        pointEdgeCollapse(vector::zero, -1, -1)
    );

    // Mark selected edges for collapse
    forAll(edges, edgeI)
    {
        if (collapseEdge[edgeI])
        {
            const edge& e = edges[edgeI];

            label masterPointI = e.start();

            // Choose the point on the edge with the highest priority.
            if (pointPriority[e.end()] > pointPriority[e.start()])
            {
                masterPointI = e.end();
            }

            label masterPointPriority = pointPriority[masterPointI];

            label index = globalPoints.toGlobal(masterPointI);

            if (!collapsePointToLocation.found(masterPointI))
            {
                const label otherVertex = e.otherVertex(masterPointI);

                if (!collapsePointToLocation.found(otherVertex))
                {
                    FatalErrorIn
                    (
                        "syncCollapse\n"
                        "(\n"
                        "   const polyMesh&,\n"
                        "   const globalIndex&,\n"
                        "   const labelList&,\n"
                        "   const PackedBoolList&,\n"
                        "   Map<point>&,\n"
                        "   List<pointEdgeCollapse>&\n"
                        ")\n"
                    )   << masterPointI << " on edge " << edgeI << " " << e
                        << " is not marked for collapse."
                        << abort(FatalError);
                }
                else
                {
                    masterPointI = otherVertex;
                    masterPointPriority = pointPriority[masterPointI];
                    index = globalPoints.toGlobal(masterPointI);
                }
            }

            const point& collapsePoint = collapsePointToLocation[masterPointI];

            const pointEdgeCollapse pec
            (
                collapsePoint,
                index,
                masterPointPriority
            );

            // Mark as collapsable but with nonsense master so it gets
            // overwritten and starts an update wave
            allEdgeInfo[edgeI] = pointEdgeCollapse
            (
                collapsePoint,
                labelMax,
                labelMin
            );

            initPointInfo.append(pec);
            initPoints.append(e.start());

            initPointInfo.append(pec);
            initPoints.append(e.end());

            nCollapsed++;
        }
    }

    PointEdgeWave<pointEdgeCollapse> collapsePropagator
    (
        mesh_,
        initPoints,
        initPointInfo,
        allPointInfo,
        allEdgeInfo,
        mesh_.globalData().nTotalPoints()  // Maximum number of iterations
    );

    return nCollapsed;
}


void Foam::edgeCollapser::filterFace
(
    const Map<DynamicList<label> >& collapseStrings,
    const List<pointEdgeCollapse>& allPointInfo,
    face& f
) const
{
    label newFp = 0;

    face oldFace = f;

    forAll(f, fp)
    {
        label pointI = f[fp];

        label collapseIndex = allPointInfo[pointI].collapseIndex();

        // Do we have a local point for this index?
        if (collapseStrings.found(collapseIndex))
        {
            label localPointI = collapseStrings[collapseIndex][0];

            if (findIndex(SubList<label>(f, newFp), localPointI) == -1)
            {
                f[newFp++] = localPointI;
            }
        }
        else if (collapseIndex == -1)
        {
            WarningIn
            (
                "filterFace"
                "(const label, const Map<DynamicList<label> >&, face&)"
            )   << "Point " << pointI << " was not visited by PointEdgeWave"
                << endl;
        }
        else
        {
            f[newFp++] = pointI;
        }
    }


    // Check for pinched face. Tries to correct
    // - consecutive duplicate vertex. Removes duplicate vertex.
    // - duplicate vertex with one other vertex in between (spike).
    // Both of these should not really occur! and should be checked before
    // collapsing edges.

    const label size = newFp;

    newFp = 2;

    for (label fp = 2; fp < size; fp++)
    {
        label fp1 = fp-1;
        label fp2 = fp-2;

        label pointI = f[fp];

        // Search for previous occurrence.
        label index = findIndex(SubList<label>(f, fp), pointI);

        if (index == fp1)
        {
            WarningIn
            (
                "Foam::edgeCollapser::filterFace(const label faceI, "
                "face& f) const"
            )   << "Removing consecutive duplicate vertex in face "
                << f << endl;
            // Don't store current pointI
        }
        else if (index == fp2)
        {
            WarningIn
            (
                "Foam::edgeCollapser::filterFace(const label faceI, "
                "face& f) const"
            )   << "Removing non-consecutive duplicate vertex in face "
                << f << endl;
            // Don't store current pointI and remove previous
            newFp--;
        }
        else if (index != -1)
        {
            WarningIn
            (
                "Foam::edgeCollapser::filterFace(const label faceI, "
                "face& f) const"
            )   << "Pinched face " << f << endl;
            f[newFp++] = pointI;
        }
        else
        {
            f[newFp++] = pointI;
        }
    }

    f.setSize(newFp);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeCollapser::edgeCollapser(const polyMesh& mesh)
:
    mesh_(mesh),
    guardFraction_(0),
    maxCollapseFaceToPointSideLengthCoeff_(0),
    allowEarlyCollapseToPoint_(false),
    allowEarlyCollapseCoeff_(0)
{}


Foam::edgeCollapser::edgeCollapser
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    guardFraction_
    (
        dict.lookupOrDefault<scalar>("guardFraction", 0)
    ),
    maxCollapseFaceToPointSideLengthCoeff_
    (
        dict.lookupOrDefault<scalar>("maxCollapseFaceToPointSideLengthCoeff", 0)
    ),
    allowEarlyCollapseToPoint_
    (
        dict.lookupOrDefault<Switch>("allowEarlyCollapseToPoint", true)
    ),
    allowEarlyCollapseCoeff_
    (
        dict.lookupOrDefault<scalar>("allowEarlyCollapseCoeff", 0)
    )
{
    if (debug)
    {
        Info<< "Edge Collapser Settings:" << nl
            << "    Guard Fraction = " << guardFraction_ << nl
            << "    Max collapse face to point side length = "
            << maxCollapseFaceToPointSideLengthCoeff_ << nl
            << "    " << (allowEarlyCollapseToPoint_ ? "Allow" : "Do not allow")
            << " early collapse to point" << nl
            << "    Early collapse coeff = " << allowEarlyCollapseCoeff_
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::edgeCollapser::setRefinement
(
    const List<pointEdgeCollapse>& allPointInfo,
    polyTopoChange& meshMod
) const
{
    const cellList& cells = mesh_.cells();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();
    const labelListList& pointFaces = mesh_.pointFaces();
    const pointZoneMesh& pointZones = mesh_.pointZones();




//    // Dump point collapses
//    label count = 0;
//    forAll(allPointInfo, ptI)
//    {
//        const pointEdgeCollapse& pec = allPointInfo[ptI];
//
//        if (mesh_.points()[ptI] != pec.collapsePoint())
//        {
//            count++;
//        }
//    }
//
//    OFstream str("collapses_" + name(count) + ".obj");
//    // Dump point collapses
//    forAll(allPointInfo, ptI)
//    {
//        const pointEdgeCollapse& pec = allPointInfo[ptI];
//
//        if
//        (
//            mesh_.points()[ptI] != pec.collapsePoint()
//         && pec.collapsePoint() != vector(GREAT, GREAT, GREAT)
//        )
//        {
//            meshTools::writeOBJ
//            (
//                str,
//                mesh_.points()[ptI],
//                pec.collapsePoint()
//            );
//        }
//    }



    bool meshChanged = false;

    PackedBoolList removedPoints(mesh_.nPoints());

    // Create strings of edges.
    // Map from collapseIndex(=global master point) to set of points
    Map<DynamicList<label> > collapseStrings;
    {
        // 1. Count elements per collapseIndex
        Map<label> nPerIndex(mesh_.nPoints()/10);
        forAll(allPointInfo, pointI)
        {
            label collapseIndex = allPointInfo[pointI].collapseIndex();

            if (collapseIndex != -1 && collapseIndex != -2)
            {
                Map<label>::iterator fnd = nPerIndex.find(collapseIndex);
                if (fnd != nPerIndex.end())
                {
                    fnd()++;
                }
                else
                {
                    nPerIndex.insert(collapseIndex, 1);
                }
            }
        }

        // 2. Size
        collapseStrings.resize(2*nPerIndex.size());
        forAllConstIter(Map<label>, nPerIndex, iter)
        {
            collapseStrings.insert(iter.key(), DynamicList<label>(iter()));
        }

        // 3. Fill
        forAll(allPointInfo, pointI)
        {
            const label collapseIndex = allPointInfo[pointI].collapseIndex();

            if (collapseIndex != -1 && collapseIndex != -2)
            {
                collapseStrings[collapseIndex].append(pointI);
            }
        }
    }




//    OFstream str2("collapseStrings_" + name(count) + ".obj");
//    // Dump point collapses
//    forAllConstIter(Map<DynamicList<label> >, collapseStrings, iter)
//    {
//        const label masterPoint = iter.key();
//        const DynamicList<label>& edgeCollapses = iter();
//
//        forAll(edgeCollapses, eI)
//        {
//            meshTools::writeOBJ
//            (
//                str2,
//                mesh_.points()[edgeCollapses[eI]],
//                mesh_.points()[masterPoint]
//            );
//        }
//    }



    // Current faces (is also collapseStatus: f.size() < 3)
    faceList newFaces(mesh_.faces());

    // Current cellCollapse status
    boolList cellRemoved(mesh_.nCells(), false);

    label nUnvisited = 0;
    label nUncollapsed = 0;
    label nCollapsed = 0;

    forAll(allPointInfo, pI)
    {
        const pointEdgeCollapse& pec = allPointInfo[pI];

        if (pec.collapseIndex() == -1)
        {
            nUnvisited++;
        }
        else if (pec.collapseIndex() == -2)
        {
            nUncollapsed++;
        }
        else
        {
            nCollapsed++;
        }
    }

    label nPoints = allPointInfo.size();

    reduce(nPoints, sumOp<label>());
    reduce(nUnvisited, sumOp<label>());
    reduce(nUncollapsed, sumOp<label>());
    reduce(nCollapsed, sumOp<label>());

    Info<< incrIndent;
    Info<< indent << "Number of points : " << nPoints << nl
        << indent << "Not visited      : " << nUnvisited << nl
        << indent << "Not collapsed    : " << nUncollapsed << nl
        << indent << "Collapsed        : " << nCollapsed << nl
        << endl;
    Info<< decrIndent;

    do
    {
        forAll(newFaces, faceI)
        {
            filterFace(collapseStrings, allPointInfo, newFaces[faceI]);
        }

        // Check if faces to be collapsed cause cells to become collapsed.
        label nCellCollapsed = 0;

        forAll(cells, cellI)
        {
            if (!cellRemoved[cellI])
            {
                const cell& cFaces = cells[cellI];

                label nFaces = cFaces.size();

                forAll(cFaces, i)
                {
                    label faceI = cFaces[i];

                    if (newFaces[faceI].size() < 3)
                    {
                        --nFaces;

                        if (nFaces < 4)
                        {
                            Pout<< "Cell:" << cellI
                                << " uses faces:" << cFaces
                                << " of which too many are marked for removal:"
                                << endl
                                << "   ";


                            forAll(cFaces, j)
                            {
                                if (newFaces[cFaces[j]].size() < 3)
                                {
                                    Pout<< ' '<< cFaces[j];
                                }
                            }
                            Pout<< endl;

                            cellRemoved[cellI] = true;

                            // Collapse all edges of cell to nothing
//                            collapseEdges(cellEdges[cellI]);

                            nCellCollapsed++;

                            break;
                        }
                    }
                }
            }
        }

        reduce(nCellCollapsed, sumOp<label>());
        Info<< indent << "Collapsing " << nCellCollapsed << " cells" << endl;

        if (nCellCollapsed == 0)
        {
            break;
        }

    } while (true);


    // Keep track of faces that have been done already.
    boolList doneFace(mesh_.nFaces(), false);

    {
        // Mark points used.
        boolList usedPoint(mesh_.nPoints(), false);

        forAll(cellRemoved, cellI)
        {
            if (cellRemoved[cellI])
            {
                meshMod.removeCell(cellI, -1);
            }
        }

        // Remove faces
        forAll(newFaces, faceI)
        {
            const face& f = newFaces[faceI];

            if (f.size() < 3)
            {
                meshMod.removeFace(faceI, -1);
                meshChanged = true;

                // Mark face as been done.
                doneFace[faceI] = true;
            }
            else
            {
                // Kept face. Mark vertices
                forAll(f, fp)
                {
                    usedPoint[f[fp]] = true;
                }
            }
        }

        // Remove unused vertices that have not been marked for removal already
        forAll(usedPoint, pointI)
        {
            if (!usedPoint[pointI])
            {
                removedPoints[pointI] = true;
                meshMod.removePoint(pointI, -1);
                meshChanged = true;
            }
        }
    }

    // Modify the point location of the remaining points
    forAll(allPointInfo, pointI)
    {
        const label collapseIndex = allPointInfo[pointI].collapseIndex();
        const point& collapsePoint = allPointInfo[pointI].collapsePoint();

        if
        (
            removedPoints[pointI] == false
         && collapseIndex != -1
         && collapseIndex != -2
        )
        {
            meshMod.modifyPoint
            (
                pointI,
                collapsePoint,
                pointZones.whichZone(pointI),
                false
            );
        }
    }


    const polyBoundaryMesh& boundaryMesh = mesh_.boundaryMesh();
    const faceZoneMesh& faceZones = mesh_.faceZones();

    // Renumber faces that use points
    forAll(allPointInfo, pointI)
    {
        if (removedPoints[pointI] == true)
        {
            const labelList& changedFaces = pointFaces[pointI];

            forAll(changedFaces, changedFaceI)
            {
                label faceI = changedFaces[changedFaceI];

                if (!doneFace[faceI])
                {
                    doneFace[faceI] = true;

                    // Get current zone info
                    label zoneID = faceZones.whichZone(faceI);

                    bool zoneFlip = false;

                    if (zoneID >= 0)
                    {
                        const faceZone& fZone = faceZones[zoneID];

                        zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                    }

                    // Get current connectivity
                    label own = faceOwner[faceI];
                    label nei = -1;
                    label patchID = -1;

                    if (mesh_.isInternalFace(faceI))
                    {
                        nei = faceNeighbour[faceI];
                    }
                    else
                    {
                        patchID = boundaryMesh.whichPatch(faceI);
                    }

                    meshMod.modifyFace
                    (
                        newFaces[faceI],            // face
                        faceI,                      // faceI to change
                        own,                        // owner
                        nei,                        // neighbour
                        false,                      // flipFaceFlux
                        patchID,                    // patch
                        zoneID,
                        zoneFlip
                    );

                    meshChanged = true;
                }
            }
        }
    }

    return meshChanged;
}


void Foam::edgeCollapser::consistentCollapse
(
    const globalIndex& globalPoints,
    const labelList& pointPriority,
    const Map<point>& collapsePointToLocation,
    PackedBoolList& collapseEdge,
    List<pointEdgeCollapse>& allPointInfo,
    const bool allowCellCollapse
) const
{
    // Make sure we don't collapse cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const faceList faces = mesh_.faces();
    const edgeList& edges = mesh_.edges();
    const labelListList& faceEdges = mesh_.faceEdges();
    const labelListList& pointEdges = mesh_.pointEdges();
    const cellList& cells = mesh_.cells();

    labelHashSet uniqueCollapses;
    labelHashSet duplicateCollapses;

    while (true)
    {
        label nUncollapsed = 0;

        syncTools::syncEdgeList
        (
            mesh_,
            collapseEdge,
            minEqOp<unsigned int>(),
            0
        );

        // Create consistent set of collapses
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Note: requires collapseEdge to be synchronised.
        syncCollapse
        (
            globalPoints,
            pointPriority,
            collapseEdge,
            collapsePointToLocation,
            allPointInfo
        );

        // Get collapsed faces

        PackedBoolList isCollapsedFace(mesh_.nFaces());
        PackedBoolList markedPoints(mesh_.nPoints());

        forAll(faces, faceI)
        {
            const face& f = faces[faceI];

            isCollapsedFace[faceI] = isFaceCollapsed(f, allPointInfo);

            if (isCollapsedFace[faceI] < 1)
            {
                determineDuplicatePointsOnFace
                (
                    f,
                    markedPoints,
                    uniqueCollapses,
                    duplicateCollapses,
                    allPointInfo
                );
            }
        }

        // Synchronise the marked points
        syncTools::syncPointList
        (
            mesh_,
            markedPoints,
            orEqOp<unsigned int>(),
            0
        );

        // Mark all edges attached to the point for collapse
        forAll(markedPoints, pointI)
        {
            if (markedPoints[pointI])
            {
                const label index = allPointInfo[pointI].collapseIndex();

                const labelList& ptEdges = pointEdges[pointI];

                forAll(ptEdges, ptEdgeI)
                {
                    const label edgeI = ptEdges[ptEdgeI];
                    const label nbrPointI = edges[edgeI].otherVertex(pointI);
                    const label nbrIndex
                        = allPointInfo[nbrPointI].collapseIndex();

                    if (collapseEdge[edgeI] && nbrIndex == index)
                    {
                        collapseEdge[edgeI] = false;
                        nUncollapsed++;
                    }
                }
            }
        }

        PackedBoolList markedEdges(mesh_.nEdges());

        if (!allowCellCollapse)
        {
            // Check collapsed cells
            forAll(cells, cellI)
            {
                const cell& cFaces = cells[cellI];

                label nFaces = cFaces.size();

                forAll(cFaces, fI)
                {
                    label faceI = cFaces[fI];

                    if (isCollapsedFace[faceI])
                    {
                        nFaces--;
                    }
                }

                if (nFaces < 4)
                {
                    forAll(cFaces, fI)
                    {
                        label faceI = cFaces[fI];

                        const labelList& fEdges = faceEdges[faceI];

                        // Unmark this face for collapse
                        forAll(fEdges, fEdgeI)
                        {
                            label edgeI = fEdges[fEdgeI];

                            if (collapseEdge[edgeI])
                            {
                                collapseEdge[edgeI] = false;
                                nUncollapsed++;
                            }

                            markedEdges[edgeI] = true;
                        }

                        // Uncollapsed this face.
                        isCollapsedFace[faceI] = false;
                        nFaces++;
                    }
                }

                if (nFaces < 4)
                {
                    FatalErrorIn
                    (
                        "consistentCollapse\n"
                        "(\n"
                        "   labelList&,\n"
                        "   Map<point>&,\n"
                        "   bool&,\n"
                        ")"
                    )   << "Cell " << cellI << " " << cFaces << nl
                        << "is " << nFaces << ", "
                        << "but cell collapse has been disabled."
                        << abort(FatalError);
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh_,
            markedEdges,
            orEqOp<unsigned int>(),
            0
        );

        nUncollapsed += breakStringsAtEdges
        (
            markedEdges,
            collapseEdge,
            allPointInfo
        );

        reduce(nUncollapsed, sumOp<label>());

        Info<< "            Uncollapsed edges = " << nUncollapsed << " / "
            << returnReduce(mesh_.nEdges(), sumOp<label>()) << endl;

        if (nUncollapsed == 0)
        {
            break;
        }
    }
}


Foam::label Foam::edgeCollapser::markSmallEdges
(
    const scalarField& minEdgeLen,
    const labelList& pointPriority,
    PackedBoolList& collapseEdge,
    Map<point>& collapsePointToLocation
) const
{
    // Work out which edges to collapse
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const pointField& points = mesh_.points();
    const edgeList& edges = mesh_.edges();

    label nCollapsed = 0;

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if (!collapseEdge[edgeI])
        {
            if (e.mag(points) < minEdgeLen[edgeI])
            {
                collapseEdge[edgeI] = true;

                label masterPointI = edgeMaster(pointPriority, e);

                if (masterPointI == -1)
                {
                    const point average
                        = 0.5*(points[e.start()] + points[e.end()]);

                    collapsePointToLocation.set(e.start(), average);
                }
                else
                {
                    const point& collapsePt = points[masterPointI];

                    collapsePointToLocation.set(masterPointI, collapsePt);
                }


                nCollapsed++;
            }
        }
    }

    return nCollapsed;
}


Foam::label Foam::edgeCollapser::markMergeEdges
(
    const scalar maxCos,
    const labelList& pointPriority,
    PackedBoolList& collapseEdge,
    Map<point>& collapsePointToLocation
) const
{
    const edgeList& edges = mesh_.edges();
    const pointField& points = mesh_.points();
    const labelListList& pointEdges = mesh_.pointEdges();

    // Point removal engine
    removePoints pointRemover(mesh_, false);

    // Find out points that can be deleted
    boolList pointCanBeDeleted;
    label nTotRemove = pointRemover.countPointUsage(maxCos, pointCanBeDeleted);

    // Rework point-to-remove into edge-to-collapse.

    label nCollapsed = 0;

    if (nTotRemove > 0)
    {
        forAll(pointEdges, pointI)
        {
            if (pointCanBeDeleted[pointI])
            {
                const labelList& pEdges = pointEdges[pointI];

                if (pEdges.size() == 2)
                {
                    // Always the case?

                    label e0 = pEdges[0];
                    label e1 = pEdges[1];

                    if (!collapseEdge[e0] && !collapseEdge[e1])
                    {
                        // Get lengths of both edges and choose the smallest
                        scalar e0length = mag
                        (
                            points[edges[e0][0]] - points[edges[e0][1]]
                        );

                        scalar e1length = mag
                        (
                            points[edges[e1][0]] - points[edges[e1][1]]
                        );

                        if (e0length <= e1length)
                        {
                            collapseEdge[e0] = true;

                            checkBoundaryPointMergeEdges
                            (
                                pointI,
                                edges[e0].otherVertex(pointI),
                                pointPriority,
                                collapsePointToLocation
                            );
                        }
                        else
                        {
                            collapseEdge[e1] = true;

                            checkBoundaryPointMergeEdges
                            (
                                pointI,
                                edges[e1].otherVertex(pointI),
                                pointPriority,
                                collapsePointToLocation
                            );
                        }

                        nCollapsed++;
                    }
                }
            }
        }
    }

    return nCollapsed;
}


Foam::labelPair Foam::edgeCollapser::markSmallSliverFaces
(
    const scalarField& faceFilterFactor,
    const labelList& pointPriority,
    PackedBoolList& collapseEdge,
    Map<point>& collapsePointToLocation
) const
{
    const faceList& faces = mesh_.faces();

    const scalarField targetFaceSizes = calcTargetFaceSizes();

    // Calculate number of faces that will be collapsed to a point or an edge
    label nCollapseToPoint = 0;
    label nCollapseToEdge = 0;

    forAll(faces, fI)
    {
        const face& f = faces[fI];

        if (faceFilterFactor[fI] <= 0)
        {
            continue;
        }

        collapseType flagCollapseFace = collapseFace
        (
            pointPriority,
            f,
            fI,
            targetFaceSizes[fI],
            collapseEdge,
            collapsePointToLocation,
            faceFilterFactor
        );

        if (flagCollapseFace == noCollapse)
        {
            continue;
        }
        else if (flagCollapseFace == toPoint)
        {
            nCollapseToPoint++;
        }
        else if (flagCollapseFace == toEdge)
        {
            nCollapseToEdge++;
        }
        else
        {
            FatalErrorIn("collapseFaces(const polyMesh&, List<labelPair>&)")
                << "Face is marked to be collapsed to " << flagCollapseFace
                << ". Currently can only collapse to point/edge."
                << abort(FatalError);
        }
    }

    return labelPair(nCollapseToPoint, nCollapseToEdge);
}


Foam::labelPair Foam::edgeCollapser::markFaceZoneEdges
(
    const faceZone& fZone,
    const scalarField& faceFilterFactor,
    const labelList& pointPriority,
    PackedBoolList& collapseEdge,
    Map<point>& collapsePointToLocation
) const
{
    const faceList& faces = mesh_.faces();

    const scalarField targetFaceSizes = calcTargetFaceSizes();

    // Calculate number of faces that will be collapsed to a point or an edge
    label nCollapseToPoint = 0;
    label nCollapseToEdge = 0;

    forAll(faces, fI)
    {
        if (fZone.whichFace(fI) == -1)
        {
            continue;
        }

        const face& f = faces[fI];

        if (faceFilterFactor[fI] <= 0)
        {
            continue;
        }

        collapseType flagCollapseFace = collapseFace
        (
            pointPriority,
            f,
            fI,
            targetFaceSizes[fI],
            collapseEdge,
            collapsePointToLocation,
            faceFilterFactor
        );

        if (flagCollapseFace == noCollapse)
        {
            continue;
        }
        else if (flagCollapseFace == toPoint)
        {
            nCollapseToPoint++;
        }
        else if (flagCollapseFace == toEdge)
        {
            nCollapseToEdge++;
        }
        else
        {
            FatalErrorIn("collapseFaces(const polyMesh&, List<labelPair>&)")
                << "Face is marked to be collapsed to " << flagCollapseFace
                << ". Currently can only collapse to point/edge."
                << abort(FatalError);
        }
    }

    return labelPair(nCollapseToPoint, nCollapseToEdge);

//    const edgeList& edges = mesh_.edges();
//    const pointField& points = mesh_.points();
//    const labelListList& edgeFaces = mesh_.edgeFaces();
//    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
//
//    forAll(edges, eI)
//    {
//        const edge& e = edges[eI];
//
//        const labelList& eFaces = edgeFaces[eI];
//
//        bool keepEdge = false;
//
//        label nInternalFaces = 0;
//        label nPatchFaces = 0;
//        label nIndirectFaces = 0;
//
//        bool coupled = false;
//
//        forAll(eFaces, eFaceI)
//        {
//            const label eFaceIndex = eFaces[eFaceI];
//
//            if (mesh_.isInternalFace(eFaceIndex))
//            {
//                nInternalFaces++;
//            }
//            else
//            {
//                const label patchIndex = bMesh.whichPatch(eFaceIndex);
//                const polyPatch& pPatch = bMesh[patchIndex];
//
//                if (pPatch.coupled())
//                {
//                    coupled = true;
//                    nInternalFaces++;
//                }
//                else
//                {
//                    // Keep the edge if an attached face is not in the zone
//                    if (fZone.whichFace(eFaceIndex) == -1)
//                    {
//                        nPatchFaces++;
//                    }
//                    else
//                    {
//                        nIndirectFaces++;
//                    }
//                }
//            }
//        }
//
//        if (eFaces.size() != nInternalFaces + nPatchFaces + nIndirectFaces)
//        {
//            Pout<< eFaces.size() << " ("
//                << nInternalFaces << "/" << nPatchFaces << "/"
//                << nIndirectFaces << ")" << endl;
//        }
//
//        if
//        (
//            eFaces.size() == nInternalFaces
//         || nIndirectFaces < (coupled ? 1 : 2)
//        )
//        {
//            keepEdge = true;
//        }
//
//        if (!keepEdge)
//        {
//            collapseEdge[eI] = true;
//
//            const Foam::point collapsePoint =
//                0.5*(points[e.end()] + points[e.start()]);
//
//            collapsePointToLocation.insert(e.start(), collapsePoint);
//            collapsePointToLocation.insert(e.end(), collapsePoint);
//        }
//    }

//    OFstream str
//    (
//        mesh_.time().path()
//       /"markedEdges_" + name(collapseEdge.count()) + ".obj"
//    );
//    label count = 0;
//
//    forAll(collapseEdge, eI)
//    {
//        if (collapseEdge[eI])
//        {
//            const edge& e = edges[eI];
//
//            meshTools::writeOBJ
//            (
//                str,
//                points[e.start()],
//                points[e.end()],
//                count
//            );
//        }
//    }
}


// ************************************************************************* //
