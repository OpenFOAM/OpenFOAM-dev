/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "meshCheck.H"
#include "PatchTools.H"
#include "polyMeshCheck.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "EdgeMap.H"
#include "wedgePolyPatch.H"
#include "polyMeshTetDecomposition.H"

#include "vtkSurfaceWriter.H"
#include "setWriter.H"
#include "writeFile.H"
#include "nonConformalCyclicPolyPatch.H"

#include "mergeAndWrite.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::meshCheck::findOppositeWedge
(
    const polyMesh& mesh,
    const wedgePolyPatch& wpp
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    scalar wppCosAngle = wpp.cosAngle();

    forAll(patches, patchi)
    {
        if
        (
            patchi != wpp.index()
         && patches[patchi].size()
         && isA<wedgePolyPatch>(patches[patchi])
        )
        {
            const wedgePolyPatch& pp =
                refCast<const wedgePolyPatch>(patches[patchi]);

            // Calculate (cos of) angle to wpp (not pp!) centre normal
            scalar ppCosAngle = wpp.centreNormal() & pp.n();

            if
            (
                pp.size() == wpp.size()
             && mag(pp.axis() & wpp.axis()) >= (1-1e-3)
             && mag(ppCosAngle - wppCosAngle) >= 1e-3
            )
            {
                return patchi;
            }
        }
    }
    return -1;
}


bool Foam::meshCheck::checkWedges
(
    const polyMesh& mesh,
    const bool report,
    const Vector<label>& directions,
    labelHashSet* setPtr
)
{
    // To mark edges without calculating edge addressing
    EdgeMap<label> edgesInError;

    const pointField& p = mesh.points();
    const faceList& fcs = mesh.faces();


    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(patches, patchi)
    {
        if (patches[patchi].size() && isA<wedgePolyPatch>(patches[patchi]))
        {
            const wedgePolyPatch& pp =
                refCast<const wedgePolyPatch>(patches[patchi]);

            scalar wedgeAngle = acos(pp.cosAngle());

            if (report)
            {
                Info<< "    Wedge " << pp.name() << " with angle "
                    << radToDeg(wedgeAngle) << " degrees"
                    << endl;
            }

            // Find opposite
            label oppositePatchi = findOppositeWedge(mesh, pp);

            if (oppositePatchi == -1)
            {
                if (report)
                {
                    Info<< " ***Cannot find opposite wedge for wedge "
                        << pp.name() << endl;
                }
                return true;
            }

            const wedgePolyPatch& opp =
                refCast<const wedgePolyPatch>(patches[oppositePatchi]);


            if (mag(opp.axis() & pp.axis()) < (1-1e-3))
            {
                if (report)
                {
                    Info<< " ***Wedges do not have the same axis."
                        << " Encountered " << pp.axis()
                        << " on patch " << pp.name()
                        << " which differs from " << opp.axis()
                        << " on opposite wedge patch" << opp.axis()
                        << endl;
                }
                return true;
            }



            // Mark edges on wedgePatches
            forAll(pp, i)
            {
                const face& f = pp[i];
                forAll(f, fp)
                {
                    label p0 = f[fp];
                    label p1 = f.nextLabel(fp);
                    edgesInError.insert(edge(p0, p1), -1);  // non-error value
                }
            }


            // Check that wedge patch is flat
            const point& p0 = p[pp.meshPoints()[0]];
            forAll(pp.meshPoints(), i)
            {
                const point& pt = p[pp.meshPoints()[i]];
                scalar d = mag((pt - p0) & pp.n());

                if (d > rootSmall)
                {
                    if (report)
                    {
                        Info<< " ***Wedge patch " << pp.name() << " not planar."
                            << " Point " << pt << " is not in patch plane by "
                            << d << " metre."
                            << endl;
                    }
                    return true;
                }
            }
        }
    }



    // Check all non-wedge faces
    label nEdgesInError = 0;

    forAll(fcs, facei)
    {
        const face& f = fcs[facei];

        forAll(f, fp)
        {
            label p0 = f[fp];
            label p1 = f.nextLabel(fp);
            if (p0 < p1)
            {
                vector d(p[p1]-p[p0]);
                scalar magD = mag(d);

                if (magD > rootVSmall)
                {
                    d /= magD;

                    // Check how many empty directions are used by the edge.
                    label nEmptyDirs = 0;
                    label nNonEmptyDirs = 0;
                    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
                    {
                        if (mag(d[cmpt]) > 1e-6)
                        {
                            if (directions[cmpt] == 0)
                            {
                                nEmptyDirs++;
                            }
                            else
                            {
                                nNonEmptyDirs++;
                            }
                        }
                    }

                    if (nEmptyDirs == 0)
                    {
                        // Purely in ok directions.
                    }
                    else if (nEmptyDirs == 1)
                    {
                        // Ok if purely in empty directions.
                        if (nNonEmptyDirs > 0)
                        {
                            if (edgesInError.insert(edge(p0, p1), facei))
                            {
                                nEdgesInError++;
                            }
                        }
                    }
                    else if (nEmptyDirs > 1)
                    {
                        // Always an error
                        if (edgesInError.insert(edge(p0, p1), facei))
                        {
                            nEdgesInError++;
                        }
                    }
                }
            }
        }
    }

    label nErrorEdges = returnReduce(nEdgesInError, sumOp<label>());

    if (nErrorEdges > 0)
    {
        if (report)
        {
            Info<< " ***Number of edges not aligned with or perpendicular to "
                << "non-empty directions: " << nErrorEdges << endl;
        }

        if (setPtr)
        {
            setPtr->resize(2*nEdgesInError);
            forAllConstIter(EdgeMap<label>, edgesInError, iter)
            {
                if (iter() >= 0)
                {
                    setPtr->insert(iter.key()[0]);
                    setPtr->insert(iter.key()[1]);
                }
            }
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    All edges aligned with or perpendicular to "
                << "non-empty directions." << endl;
        }
        return false;
    }
}


namespace Foam
{
    //- Default transformation behaviour for position
    class transformPositionList
    {
    public:

        //- Transform patch-based field
        void operator()
        (
            const coupledPolyPatch& cpp,
            List<pointField>& pts
        ) const
        {
            // Each element of pts is all the points in the face. Convert into
            // lists of size cpp to transform.

            List<pointField> newPts(pts.size());
            forAll(pts, facei)
            {
                newPts[facei].setSize(pts[facei].size());
            }

            label index = 0;
            while (true)
            {
                label n = 0;

                // Extract for every face the i'th position
                pointField ptsAtIndex(pts.size(), Zero);
                forAll(cpp, facei)
                {
                    const pointField& facePts = pts[facei];
                    if (facePts.size() > index)
                    {
                        ptsAtIndex[facei] = facePts[index];
                        n++;
                    }
                }

                if (n == 0)
                {
                    break;
                }

                // Now ptsAtIndex will have for every face either zero or
                // the position of the i'th vertex. Transform.
                cpp.transform().transformPosition(ptsAtIndex, ptsAtIndex);

                // Extract back from ptsAtIndex into newPts
                forAll(cpp, facei)
                {
                    pointField& facePts = newPts[facei];
                    if (facePts.size() > index)
                    {
                        facePts[index] = ptsAtIndex[facei];
                    }
                }

                index++;
            }

            pts.transfer(newPts);
        }
    };
}


bool Foam::meshCheck::checkCoupledPoints
(
    const polyMesh& mesh,
    const bool report,
    labelHashSet* setPtr
)
{
    const pointField& p = mesh.points();
    const faceList& fcs = mesh.faces();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Zero'th point on coupled faces
    // pointField nbrZeroPoint(fcs.size()-mesh.nInternalFaces(), vector::max);
    List<pointField> nbrPoints(fcs.size() - mesh.nInternalFaces());

    // Exchange zero point
    forAll(patches, patchi)
    {
        if (patches[patchi].coupled())
        {
            const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>
            (
                patches[patchi]
            );

            forAll(cpp, i)
            {
                label bFacei = cpp.start() + i - mesh.nInternalFaces();
                const face& f = cpp[i];
                nbrPoints[bFacei].setSize(f.size());
                forAll(f, fp)
                {
                    const point& p0 = p[f[fp]];
                    nbrPoints[bFacei][fp] = p0;
                }
            }
        }
    }
    syncTools::syncBoundaryFaceList
    (
        mesh,
        nbrPoints,
        eqOp<pointField>(),
        transformPositionList()
    );

    // Compare to local ones. Use same tolerance as for matching
    label nErrorFaces = 0;
    scalar avgMismatch = 0;
    label nCoupledPoints = 0;

    forAll(patches, patchi)
    {
        if (patches[patchi].coupled())
        {
            const coupledPolyPatch& cpp =
                refCast<const coupledPolyPatch>(patches[patchi]);

            if (cpp.owner())
            {
                scalarField smallDist
                (
                    cpp.calcFaceTol
                    (
                        // cpp.matchTolerance(),
                        cpp,
                        cpp.points(),
                        cpp.faceCentres()
                    )
                );

                forAll(cpp, i)
                {
                    label bFacei = cpp.start() + i - mesh.nInternalFaces();
                    const face& f = cpp[i];

                    if (f.size() != nbrPoints[bFacei].size())
                    {
                        FatalErrorInFunction
                            << "Local face size : " << f.size()
                            << " does not equal neighbour face size : "
                            << nbrPoints[bFacei].size()
                            << abort(FatalError);
                    }

                    label fp = 0;
                    forAll(f, j)
                    {
                        const point& p0 = p[f[fp]];
                        scalar d = mag(p0 - nbrPoints[bFacei][j]);

                        if (d > smallDist[i])
                        {
                            if (setPtr)
                            {
                                setPtr->insert(cpp.start()+i);
                            }
                            nErrorFaces++;

                            break;
                        }

                        avgMismatch += d;
                        nCoupledPoints++;

                        fp = f.rcIndex(fp);
                    }
                }
            }
        }
    }

    reduce(nErrorFaces, sumOp<label>());
    reduce(avgMismatch, maxOp<scalar>());
    reduce(nCoupledPoints, sumOp<label>());

    if (nCoupledPoints > 0)
    {
        avgMismatch /= nCoupledPoints;
    }

    if (nErrorFaces > 0)
    {
        if (report)
        {
            Info<< "  **Error in coupled point location: "
                << nErrorFaces
                << " faces have their 0th or consecutive vertex not opposite"
                << " their coupled equivalent. Average mismatch "
                << avgMismatch << "."
                << endl;
        }

        return true;
    }
    else
    {
        if (report)
        {
            Info<< "    Coupled point location match (average "
                << avgMismatch << ") OK." << endl;
        }

        return false;
    }
}


Foam::label Foam::meshCheck::checkGeometry
(
    const polyMesh& mesh,
    const bool allGeometry,
    const scalar nonOrthThreshold,
    const scalar skewThreshold,
    const autoPtr<surfaceWriter>& surfWriter,
    const autoPtr<Foam::setWriter>& setWriter
)
{
    const scalar closedThreshold = 1.0e-6;
    const scalar aspectThreshold = 1000;
    const scalar planarCosAngle = 1.0e-6;

    label noFailedChecks = 0;

    Info<< "\nChecking geometry..." << endl;

    // Get a small relative length from the bounding box
    const boundBox& globalBb = mesh.bounds();

    Info<< "    Overall domain bounding box "
        << globalBb.min() << " " << globalBb.max() << endl;

    if (allGeometry)
    {
        Info<< "    Patch bounding boxes" << endl;

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            if (!isA<processorPolyPatch>(pp))
            {
                Info<< "        " << setw(20) << pp.name();

                const pointField& pts = pp.points();
                const labelList& mp = pp.meshPoints();

                if (returnReduce(mp.size(), sumOp<label>()) > 0)
                {
                    boundBox bb(point::max, point::min);
                    forAll(mp, i)
                    {
                        bb.min() = min(bb.min(), pts[mp[i]]);
                        bb.max() = max(bb.max(), pts[mp[i]]);
                    }
                    reduce(bb.min(), minOp<vector>());
                    reduce(bb.max(), maxOp<vector>());
                    Info<< ' ' << bb;
                }
            }
            Info<< endl;
        }
    }

    // Min length
    scalar minDistSqr = magSqr(1e-6 * globalBb.span());

    // Geometric directions
    const Vector<label> validDirs = (mesh.geometricD() + Vector<label>::one)/2;
    Info<< "    Mesh has " << mesh.nGeometricD()
        << " geometric (non-empty/wedge) directions " << validDirs << endl;

    // Solution directions
    const Vector<label> solDirs = (mesh.solutionD() + Vector<label>::one)/2;
    Info<< "    Mesh has " << mesh.nSolutionD()
        << " solution (non-empty) directions " << solDirs << endl;

    if (mesh.nGeometricD() < 3)
    {
        pointSet nonAlignedPoints(mesh, "nonAlignedEdges", mesh.nPoints()/100);

        if
        (
            (
                validDirs != solDirs
             && checkWedges(mesh, true, validDirs, &nonAlignedPoints)
            )
         || (
                validDirs == solDirs
             && meshCheck::checkEdgeAlignment
                (
                    mesh,
                    true,
                    validDirs,
                    &nonAlignedPoints
                )
            )
        )
        {
            noFailedChecks++;
            label nNonAligned = returnReduce
            (
                nonAlignedPoints.size(),
                sumOp<label>()
            );

            if (nNonAligned > 0)
            {
                if (setWriter.valid())
                {
                    Info<< "  <<Writing " << nNonAligned
                        << " points on non-aligned edges to set "
                        << nonAlignedPoints.name() << endl;
                    nonAlignedPoints.instance() = mesh.pointsInstance();
                    nonAlignedPoints.write();
                    meshCheck::mergeAndWrite(setWriter, nonAlignedPoints);
                }
            }
        }
    }

    if (meshCheck::checkClosedBoundary(mesh, true)) noFailedChecks++;

    {
        cellSet cells(mesh, "nonClosedCells", mesh.nCells()/100+1);
        cellSet aspectCells(mesh, "highAspectRatioCells", mesh.nCells()/100+1);
        if
        (
            meshCheck::checkClosedCells
            (
                mesh,
                closedThreshold,
                aspectThreshold,
                true,
                &cells,
                &aspectCells,
                mesh.geometricD()
            )
        )
        {
            noFailedChecks++;

            label nNonClosed = returnReduce(cells.size(), sumOp<label>());

            if (nNonClosed > 0)
            {
                if (setWriter.valid())
                {
                    Info<< "  <<Writing " << nNonClosed
                        << " non closed cells to set " << cells.name() << endl;
                    cells.instance() = mesh.pointsInstance();
                    cells.write();
                }

                if (surfWriter.valid())
                {
                    meshCheck::mergeAndWrite(surfWriter(), cells);
                }
            }
        }

        label nHighAspect = returnReduce(aspectCells.size(), sumOp<label>());

        if (nHighAspect > 0)
        {
            if (setWriter.valid())
            {
                Info<< "  <<Writing " << nHighAspect
                    << " cells with high aspect ratio to set "
                    << aspectCells.name() << endl;
                aspectCells.instance() = mesh.pointsInstance();
                aspectCells.write();
            }

            if (surfWriter.valid())
            {
                meshCheck::mergeAndWrite(surfWriter(), aspectCells);
            }
        }
    }

    {
        faceSet faces(mesh, "zeroAreaFaces", mesh.nFaces()/100+1);
        if (meshCheck::checkFaceAreas(mesh, true, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                if (setWriter.valid())
                {
                    Info<< "  <<Writing " << nFaces
                        << " zero area faces to set " << faces.name() << endl;
                    faces.instance() = mesh.pointsInstance();
                    faces.write();
                }

                if (surfWriter.valid())
                {
                    meshCheck::mergeAndWrite(surfWriter(), faces);
                }
            }
        }
    }

    {
        cellSet cells(mesh, "zeroVolumeCells", mesh.nCells()/100+1);
        if (meshCheck::checkCellVolumes(mesh, true, &cells))
        {
            noFailedChecks++;

            label nCells = returnReduce(cells.size(), sumOp<label>());

            if (nCells > 0)
            {
                if (setWriter.valid())
                {
                    Info<< "  <<Writing " << nCells
                        << " zero volume cells to set " << cells.name() << endl;
                    cells.instance() = mesh.pointsInstance();
                    cells.write();
                }

                if (surfWriter.valid())
                {
                    meshCheck::mergeAndWrite(surfWriter(), cells);
                }
            }
        }
    }

    {
        faceSet faces(mesh, "nonOrthoFaces", mesh.nFaces()/100+1);
        if
        (
            meshCheck::checkFaceOrthogonality
            (
                mesh,
                nonOrthThreshold,
                true,
                &faces
            )
        )
        {
            noFailedChecks++;
        }

        label nFaces = returnReduce(faces.size(), sumOp<label>());

        if (nFaces > 0)
        {
            if (setWriter.valid())
            {
                Info<< "  <<Writing " << nFaces
                    << " non-orthogonal faces to set " << faces.name() << endl;
                faces.instance() = mesh.pointsInstance();
                faces.write();
            }

            if (surfWriter.valid())
            {
                meshCheck::mergeAndWrite(surfWriter(), faces);
            }
        }
    }

    {
        faceSet faces(mesh, "wrongOrientedFaces", mesh.nFaces()/100 + 1);
        if (meshCheck::checkFacePyramids(mesh, true, -small, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                if (setWriter.valid())
                {
                    Info<< "  <<Writing " << nFaces
                        << " faces with incorrect orientation to set "
                        << faces.name() << endl;
                    faces.instance() = mesh.pointsInstance();
                    faces.write();
                }

                if (surfWriter.valid())
                {
                    meshCheck::mergeAndWrite(surfWriter(), faces);
                }
            }
        }
    }

    {
        faceSet faces(mesh, "skewFaces", mesh.nFaces()/100+1);
        if
        (
            meshCheck::checkFaceSkewness
            (
                mesh,
                skewThreshold,
                true,
                &faces
            )
        )
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                if (setWriter.valid())
                {
                    Info<< "  <<Writing " << nFaces
                        << " skew faces to set " << faces.name() << endl;
                    faces.instance() = mesh.pointsInstance();
                    faces.write();
                }

                if (surfWriter.valid())
                {
                    meshCheck::mergeAndWrite(surfWriter(), faces);
                }
            }
        }
    }

    {
        faceSet faces(mesh, "coupledFaces", mesh.nFaces()/100 + 1);
        if (checkCoupledPoints(mesh, true, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                if (setWriter.valid())
                {
                    Info<< "  <<Writing " << nFaces
                        << " faces with incorrectly matched 0th "
                           "(or consecutive) vertex to set "
                        << faces.name() << endl;
                    faces.instance() = mesh.pointsInstance();
                    faces.write();
                }

                if (surfWriter.valid())
                {
                    meshCheck::mergeAndWrite(surfWriter(), faces);
                }
            }
        }
    }

    if (allGeometry)
    {
        faceSet faces(mesh, "lowQualityTetFaces", mesh.nFaces()/100+1);
        if
        (
            polyMeshTetDecomposition::checkFaceTets
            (
                mesh,
                polyMeshTetDecomposition::minTetQuality,
                true,
                &faces
            )
        )
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                if (setWriter.valid())
                {
                    Info<< "  <<Writing " << nFaces
                        << " faces with low quality or negative volume "
                        << "decomposition tets to set " << faces.name() << endl;
                    faces.instance() = mesh.pointsInstance();
                    faces.write();
                }

                if (surfWriter.valid())
                {
                    meshCheck::mergeAndWrite(surfWriter(), faces);
                }
            }
        }
    }

    if (allGeometry)
    {
        // Note use of nPoints since don't want edge construction.
        pointSet points(mesh, "shortEdges", mesh.nPoints()/1000 + 1);
        if (meshCheck::checkEdgeLength(mesh, true, minDistSqr, &points))
        {
            // noFailedChecks++;

            label nPoints = returnReduce(points.size(), sumOp<label>());

            if (nPoints > 0)
            {
                if (setWriter.valid())
                {
                    Info<< "  <<Writing " << nPoints
                        << " points on short edges to set " << points.name()
                        << endl;
                    points.instance() = mesh.pointsInstance();
                    points.write();
                    meshCheck::mergeAndWrite(setWriter, points);
                }
            }
        }

        label nEdgeClose = returnReduce(points.size(), sumOp<label>());

        if (meshCheck::checkPointNearness(mesh, false, minDistSqr, &points))
        {
            // noFailedChecks++;

            label nPoints = returnReduce(points.size(), sumOp<label>());

            if (nPoints > nEdgeClose)
            {
                if (setWriter.valid())
                {
                    pointSet nearPoints(mesh, "nearPoints", points);
                    Info<< "  <<Writing " << nPoints
                        << " near (closer than " << Foam::sqrt(minDistSqr)
                        << " apart) points to set " << nearPoints.name()
                        << endl;
                    nearPoints.instance() = mesh.pointsInstance();
                    nearPoints.write();
                    meshCheck::mergeAndWrite(setWriter, nearPoints);
                }
            }
        }
    }

    if (allGeometry)
    {
        faceSet faces(mesh, "concaveFaces", mesh.nFaces()/100 + 1);
        if (meshCheck::checkFaceAngles(mesh, true, degToRad(10), &faces))
        {
            // noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                if (setWriter.valid())
                {
                    Info<< "  <<Writing " << nFaces
                        << " faces with concave angles to set " << faces.name()
                        << endl;
                    faces.instance() = mesh.pointsInstance();
                    faces.write();
                }

                if (surfWriter.valid())
                {
                    meshCheck::mergeAndWrite(surfWriter(), faces);
                }
            }
        }
    }

    if (allGeometry)
    {
        faceSet faces(mesh, "warpedFaces", mesh.nFaces()/100 + 1);
        if (meshCheck::checkFaceFlatness(mesh, true, 0.8, &faces))
        {
            // noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                if (setWriter.valid())
                {
                    Info<< "  <<Writing " << nFaces
                        << " warped faces to set " << faces.name() << endl;
                    faces.instance() = mesh.pointsInstance();
                    faces.write();
                }

                if (surfWriter.valid())
                {
                    meshCheck::mergeAndWrite(surfWriter(), faces);
                }
            }
        }
    }

    if (allGeometry)
    {
        cellSet cells(mesh, "underdeterminedCells", mesh.nCells()/100);
        if (meshCheck::checkCellDeterminant(mesh, true, &cells))
        {
            noFailedChecks++;

            if (setWriter.valid())
            {
                label nCells = returnReduce(cells.size(), sumOp<label>());

                Info<< "  <<Writing " << nCells
                    << " under-determined cells to set " << cells.name()
                    << endl;
                cells.instance() = mesh.pointsInstance();
                cells.write();
            }

            if (surfWriter.valid())
            {
                meshCheck::mergeAndWrite(surfWriter(), cells);
            }
        }
    }

    if (allGeometry)
    {
        cellSet cells(mesh, "concaveCells", mesh.nCells()/100);
        if (meshCheck::checkConcaveCells(mesh, planarCosAngle, true, &cells))
        {
            noFailedChecks++;

            if (setWriter.valid())
            {
                label nCells = returnReduce(cells.size(), sumOp<label>());

                Info<< "  <<Writing " << nCells
                    << " concave cells to set " << cells.name() << endl;
                cells.instance() = mesh.pointsInstance();
                cells.write();
            }

            if (surfWriter.valid())
            {
                meshCheck::mergeAndWrite(surfWriter(), cells);
            }
        }
    }

    if (allGeometry)
    {
        faceSet faces(mesh, "lowWeightFaces", mesh.nFaces()/100);
        if (meshCheck::checkFaceWeight(mesh, true, 0.05, &faces))
        {
            noFailedChecks++;

            if (setWriter.valid())
            {
                label nFaces = returnReduce(faces.size(), sumOp<label>());

                Info<< "  <<Writing " << nFaces
                    << " faces with low interpolation weights to set "
                    << faces.name() << endl;
                faces.instance() = mesh.pointsInstance();
                faces.write();
            }

            if (surfWriter.valid())
            {
                meshCheck::mergeAndWrite(surfWriter(), faces);
            }
        }
    }

    if (allGeometry)
    {
        faceSet faces(mesh, "lowVolRatioFaces", mesh.nFaces()/100);
        if (meshCheck::checkVolRatio(mesh, true, 0.01, &faces))
        {
            noFailedChecks++;

            if (setWriter.valid())
            {
                label nFaces = returnReduce(faces.size(), sumOp<label>());

                Info<< "  <<Writing " << nFaces
                    << " faces with low volume ratio cells to set "
                    << faces.name() << endl;
                faces.instance() = mesh.pointsInstance();
                faces.write();
            }

            if (surfWriter.valid())
            {
                meshCheck::mergeAndWrite(surfWriter(), faces);
            }
        }
    }

    if (allGeometry)
    {
        const fileName outputPath =
            mesh.time().globalPath()
           /functionObjects::writeFile::outputPrefix
           /(mesh.name() != polyMesh::defaultRegion ? mesh.name() : word())
           /"checkMesh"
           /mesh.time().name();

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Compute coverage for all orig patches
        PtrList<scalarField> patchCoverage(patches.size());
        forAll(patches, nccPatchi)
        {
            if (isA<nonConformalCyclicPolyPatch>(patches[nccPatchi]))
            {
                const nonConformalCyclicPolyPatch& nccPp =
                    refCast<const nonConformalCyclicPolyPatch>
                    (patches[nccPatchi]);

                if (nccPp.owner())
                {
                    const polyPatch& origPp = nccPp.origPatch();
                    const polyPatch& nbrOrigPp = nccPp.nbrPatch().origPatch();

                    const patchToPatches::intersection& intersection =
                        nccPp.intersection();

                    if (!patchCoverage.set(origPp.index()))
                    {
                        patchCoverage.set
                        (
                            origPp.index(),
                            scalarField(origPp.size(), 0)
                        );
                    }

                    patchCoverage[origPp.index()] +=
                        intersection.srcCoverage();

                    if (!patchCoverage.set(nbrOrigPp.index()))
                    {
                        patchCoverage.set
                        (
                            nbrOrigPp.index(),
                            scalarField(nbrOrigPp.size(), 0)
                        );
                    }

                    patchCoverage[nbrOrigPp.index()] +=
                        intersection.tgtCoverage();
                }
            }
        }

        // Write out to surface files
        forAll(patches, patchi)
        {
            if (patchCoverage.set(patchi))
            {
                const polyPatch& patch = patches[patchi];

                // Collect the patch geometry
                labelList pointToGlobal;
                labelList uniqueMeshPointLabels;
                autoPtr<globalIndex> globalPoints;
                autoPtr<globalIndex> globalFaces;
                faceList mergedFaces;
                pointField mergedPoints;
                Foam::PatchTools::gatherAndMerge
                (
                    mesh,
                    patch.localFaces(),
                    patch.meshPoints(),
                    patch.meshPointMap(),
                    pointToGlobal,
                    uniqueMeshPointLabels,
                    globalPoints,
                    globalFaces,
                    mergedFaces,
                    mergedPoints
                );

                // Collect the patch coverage
                scalarField mergedCoverage;
                globalFaces().gather
                (
                    UPstream::worldComm,
                    labelList(UPstream::procID(UPstream::worldComm)),
                    patchCoverage[patchi],
                    mergedCoverage
                );

                // Write the surface
                if (Pstream::master())
                {
                    vtkSurfaceWriter
                    (
                        mesh.time().writeFormat(),
                        mesh.time().writeCompression()
                    ).write
                    (
                        outputPath,
                        patch.name() + "_coverage",
                        mergedPoints,
                        mergedFaces,
                        false,
                        "coverage",
                        mergedCoverage
                    );
                }
            }
        }
    }

    return noFailedChecks;
}


// ************************************************************************* //
