#include "checkGeometry.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "EdgeMap.H"
#include "wedgePolyPatch.H"
#include "unitConversion.H"
#include "polyMeshTetDecomposition.H"


// Find wedge with opposite orientation. Note: does not actually check that
// it is opposite, only that it has opposite normal and same axis
Foam::label Foam::findOppositeWedge
(
    const polyMesh& mesh,
    const wedgePolyPatch& wpp
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    scalar wppCosAngle = wpp.cosAngle();

    forAll(patches, patchI)
    {
        if
        (
            patchI != wpp.index()
         && patches[patchI].size()
         && isA<wedgePolyPatch>(patches[patchI])
        )
        {
            const wedgePolyPatch& pp =
                refCast<const wedgePolyPatch>(patches[patchI]);

            // Calculate (cos of) angle to wpp (not pp!) centre normal
            scalar ppCosAngle = wpp.centreNormal() & pp.n();

            if
            (
                pp.size() == wpp.size()
             && mag(pp.axis() & wpp.axis()) >= (1-1e-3)
             && mag(ppCosAngle - wppCosAngle) >= 1e-3
            )
            {
                return patchI;
            }
        }
    }
    return -1;
}


bool Foam::checkWedges
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
    forAll(patches, patchI)
    {
        if (patches[patchI].size() && isA<wedgePolyPatch>(patches[patchI]))
        {
            const wedgePolyPatch& pp =
                refCast<const wedgePolyPatch>(patches[patchI]);

            scalar wedgeAngle = acos(pp.cosAngle());

            if (report)
            {
                Info<< "    Wedge " << pp.name() << " with angle "
                    << radToDeg(wedgeAngle) << " degrees"
                    << endl;
            }

            // Find opposite
            label oppositePatchI = findOppositeWedge(mesh, pp);

            if (oppositePatchI == -1)
            {
                if (report)
                {
                    Info<< " ***Cannot find opposite wedge for wedge "
                        << pp.name() << endl;
                }
                return true;
            }

            const wedgePolyPatch& opp =
                refCast<const wedgePolyPatch>(patches[oppositePatchI]);


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

                if (d > sqrt(SMALL))
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

    forAll(fcs, faceI)
    {
        const face& f = fcs[faceI];

        forAll(f, fp)
        {
            label p0 = f[fp];
            label p1 = f.nextLabel(fp);
            if (p0 < p1)
            {
                vector d(p[p1]-p[p0]);
                scalar magD = mag(d);

                if (magD > ROOTVSMALL)
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
                            if (edgesInError.insert(edge(p0, p1), faceI))
                            {
                                nEdgesInError++;
                            }
                        }
                    }
                    else if (nEmptyDirs > 1)
                    {
                        // Always an error
                        if (edgesInError.insert(edge(p0, p1), faceI))
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
            forAll(pts, faceI)
            {
                newPts[faceI].setSize(pts[faceI].size());
            }

            label index = 0;
            while (true)
            {
                label n = 0;

                // Extract for every face the i'th position
                pointField ptsAtIndex(pts.size(), vector::zero);
                forAll(cpp, faceI)
                {
                    const pointField& facePts = pts[faceI];
                    if (facePts.size() > index)
                    {
                        ptsAtIndex[faceI] = facePts[index];
                        n++;
                    }
                }

                if (n == 0)
                {
                    break;
                }

                // Now ptsAtIndex will have for every face either zero or
                // the position of the i'th vertex. Transform.
                cpp.transformPosition(ptsAtIndex);

                // Extract back from ptsAtIndex into newPts
                forAll(cpp, faceI)
                {
                    pointField& facePts = newPts[faceI];
                    if (facePts.size() > index)
                    {
                        facePts[index] = ptsAtIndex[faceI];
                    }
                }

                index++;
            }

            pts.transfer(newPts);
        }
    };
}


bool Foam::checkCoupledPoints
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
    //pointField nbrZeroPoint(fcs.size()-mesh.nInternalFaces(), vector::max);
    List<pointField> nbrPoints(fcs.size() - mesh.nInternalFaces());

    // Exchange zero point
    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>
            (
                patches[patchI]
            );

            forAll(cpp, i)
            {
                label bFaceI = cpp.start() + i - mesh.nInternalFaces();
                const face& f = cpp[i];
                nbrPoints[bFaceI].setSize(f.size());
                forAll(f, fp)
                {
                    const point& p0 = p[f[fp]];
                    nbrPoints[bFaceI][fp] = p0;
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

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            const coupledPolyPatch& cpp =
                refCast<const coupledPolyPatch>(patches[patchI]);

            if (cpp.owner())
            {
                scalarField smallDist
                (
                    cpp.calcFaceTol
                    (
                        //cpp.matchTolerance(),
                        cpp,
                        cpp.points(),
                        cpp.faceCentres()
                    )
                );

                forAll(cpp, i)
                {
                    label bFaceI = cpp.start() + i - mesh.nInternalFaces();
                    const face& f = cpp[i];

                    if (f.size() != nbrPoints[bFaceI].size())
                    {
                        FatalErrorIn
                        (
                            "Foam::checkCoupledPoints\n"
                            "(\n"
                            "   const polyMesh&, const bool, labelHashSet*\n"
                            ")\n"
                        )   << "Local face size : " << f.size()
                            << " does not equal neighbour face size : "
                            << nbrPoints[bFaceI].size()
                            << abort(FatalError);
                    }

                    label fp = 0;
                    forAll(f, j)
                    {
                        const point& p0 = p[f[fp]];
                        scalar d = mag(p0 - nbrPoints[bFaceI][j]);

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


Foam::label Foam::checkGeometry(const polyMesh& mesh, const bool allGeometry)
{
    label noFailedChecks = 0;

    Info<< "\nChecking geometry..." << endl;

    // Get a small relative length from the bounding box
    const boundBox& globalBb = mesh.bounds();

    Info<< "    Overall domain bounding box "
        << globalBb.min() << " " << globalBb.max() << endl;


    // Min length
    scalar minDistSqr = magSqr(1e-6 * globalBb.span());

    // Non-empty directions
    const Vector<label> validDirs = (mesh.geometricD() + Vector<label>::one)/2;
    Info<< "    Mesh (non-empty, non-wedge) directions " << validDirs << endl;

    const Vector<label> solDirs = (mesh.solutionD() + Vector<label>::one)/2;
    Info<< "    Mesh (non-empty) directions " << solDirs << endl;

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
             && mesh.checkEdgeAlignment(true, validDirs, &nonAlignedPoints)
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
                Info<< "  <<Writing " << nNonAligned
                    << " points on non-aligned edges to set "
                    << nonAlignedPoints.name() << endl;
                nonAlignedPoints.instance() = mesh.pointsInstance();
                nonAlignedPoints.write();
            }
        }
    }

    if (mesh.checkClosedBoundary(true)) noFailedChecks++;

    {
        cellSet cells(mesh, "nonClosedCells", mesh.nCells()/100+1);
        cellSet aspectCells(mesh, "highAspectRatioCells", mesh.nCells()/100+1);
        if
        (
            mesh.checkClosedCells
            (
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
                Info<< "  <<Writing " << nNonClosed
                    << " non closed cells to set " << cells.name() << endl;
                cells.instance() = mesh.pointsInstance();
                cells.write();
            }
        }

        label nHighAspect = returnReduce(aspectCells.size(), sumOp<label>());

        if (nHighAspect > 0)
        {
            Info<< "  <<Writing " << nHighAspect
                << " cells with high aspect ratio to set "
                << aspectCells.name() << endl;
            aspectCells.instance() = mesh.pointsInstance();
            aspectCells.write();
        }
    }

    {
        faceSet faces(mesh, "zeroAreaFaces", mesh.nFaces()/100+1);
        if (mesh.checkFaceAreas(true, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                Info<< "  <<Writing " << nFaces
                    << " zero area faces to set " << faces.name() << endl;
                faces.instance() = mesh.pointsInstance();
                faces.write();
            }
        }
    }

    {
        cellSet cells(mesh, "zeroVolumeCells", mesh.nCells()/100+1);
        if (mesh.checkCellVolumes(true, &cells))
        {
            noFailedChecks++;

            label nCells = returnReduce(cells.size(), sumOp<label>());

            if (nCells > 0)
            {
                Info<< "  <<Writing " << nCells
                    << " zero volume cells to set " << cells.name() << endl;
                cells.instance() = mesh.pointsInstance();
                cells.write();
            }
        }
    }

    {
        faceSet faces(mesh, "nonOrthoFaces", mesh.nFaces()/100+1);
        if (mesh.checkFaceOrthogonality(true, &faces))
        {
            noFailedChecks++;
        }

        label nFaces = returnReduce(faces.size(), sumOp<label>());

        if (nFaces > 0)
        {
            Info<< "  <<Writing " << nFaces
                << " non-orthogonal faces to set " << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
        }
    }

    {
        faceSet faces(mesh, "wrongOrientedFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFacePyramids(true, -SMALL, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                Info<< "  <<Writing " << nFaces
                    << " faces with incorrect orientation to set "
                    << faces.name() << endl;
                faces.instance() = mesh.pointsInstance();
                faces.write();
            }
        }
    }

    {
        faceSet faces(mesh, "skewFaces", mesh.nFaces()/100+1);
        if (mesh.checkFaceSkewness(true, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                Info<< "  <<Writing " << nFaces
                    << " skew faces to set " << faces.name() << endl;
                faces.instance() = mesh.pointsInstance();
                faces.write();
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
                Info<< "  <<Writing " << nFaces
                    << " faces with incorrectly matched 0th (or consecutive)"
                    << " vertex to set "
                    << faces.name() << endl;
                faces.instance() = mesh.pointsInstance();
                faces.write();
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
                Info<< "  <<Writing " << nFaces
                    << " faces with low quality or negative volume "
                    << "decomposition tets to set " << faces.name() << endl;
                faces.instance() = mesh.pointsInstance();
                faces.write();
            }
        }
    }

    if (allGeometry)
    {
        // Note use of nPoints since don't want edge construction.
        pointSet points(mesh, "shortEdges", mesh.nPoints()/1000 + 1);
        if (mesh.checkEdgeLength(true, minDistSqr, &points))
        {
            //noFailedChecks++;

            label nPoints = returnReduce(points.size(), sumOp<label>());

            if (nPoints > 0)
            {
                Info<< "  <<Writing " << nPoints
                    << " points on short edges to set " << points.name()
                    << endl;
                points.instance() = mesh.pointsInstance();
                points.write();
            }
        }

        label nEdgeClose = returnReduce(points.size(), sumOp<label>());

        if (mesh.checkPointNearness(false, minDistSqr, &points))
        {
            //noFailedChecks++;

            label nPoints = returnReduce(points.size(), sumOp<label>());

            if (nPoints > nEdgeClose)
            {
                pointSet nearPoints(mesh, "nearPoints", points);
                Info<< "  <<Writing " << nPoints
                    << " near (closer than " << Foam::sqrt(minDistSqr)
                    << " apart) points to set " << nearPoints.name() << endl;
                nearPoints.instance() = mesh.pointsInstance();
                nearPoints.write();
            }
        }
    }

    if (allGeometry)
    {
        faceSet faces(mesh, "concaveFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceAngles(true, 10, &faces))
        {
            //noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                Info<< "  <<Writing " << nFaces
                    << " faces with concave angles to set " << faces.name()
                    << endl;
                faces.instance() = mesh.pointsInstance();
                faces.write();
            }
        }
    }

    if (allGeometry)
    {
        faceSet faces(mesh, "warpedFaces", mesh.nFaces()/100 + 1);
        if (mesh.checkFaceFlatness(true, 0.8, &faces))
        {
            //noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            if (nFaces > 0)
            {
                Info<< "  <<Writing " << nFaces
                    << " warped faces to set " << faces.name() << endl;
                faces.instance() = mesh.pointsInstance();
                faces.write();
            }
        }
    }

    if (allGeometry)
    {
        cellSet cells(mesh, "underdeterminedCells", mesh.nCells()/100);
        if (mesh.checkCellDeterminant(true, &cells))
        {
            noFailedChecks++;

            label nCells = returnReduce(cells.size(), sumOp<label>());

            Info<< "  <<Writing " << nCells
                << " under-determined cells to set " << cells.name() << endl;
            cells.instance() = mesh.pointsInstance();
            cells.write();
        }
    }

    if (allGeometry)
    {
        cellSet cells(mesh, "concaveCells", mesh.nCells()/100);
        if (mesh.checkConcaveCells(true, &cells))
        {
            noFailedChecks++;

            label nCells = returnReduce(cells.size(), sumOp<label>());

            Info<< "  <<Writing " << nCells
                << " concave cells to set " << cells.name() << endl;
            cells.instance() = mesh.pointsInstance();
            cells.write();
        }
    }

    if (allGeometry)
    {
        faceSet faces(mesh, "lowWeightFaces", mesh.nFaces()/100);
        if (mesh.checkFaceWeight(true, 0.05, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            Info<< "  <<Writing " << nFaces
                << " faces with low interpolation weights to set "
                << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
        }
    }

    if (allGeometry)
    {
        faceSet faces(mesh, "lowVolRatioFaces", mesh.nFaces()/100);
        if (mesh.checkVolRatio(true, 0.01, &faces))
        {
            noFailedChecks++;

            label nFaces = returnReduce(faces.size(), sumOp<label>());

            Info<< "  <<Writing " << nFaces
                << " faces with low volume ratio cells to set "
                << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
        }
    }

    return noFailedChecks;
}
