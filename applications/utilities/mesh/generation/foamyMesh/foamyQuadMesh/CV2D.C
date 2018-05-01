/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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

\*----------------------------------------------------------------------------*/

#include "CV2D.H"
#include "Random.H"
#include "transform.H"
#include "IFstream.H"
#include "uint.H"

namespace Foam
{
    defineTypeNameAndDebug(CV2D, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CV2D::insertBoundingBox()
{
    Info<< "insertBoundingBox: creating bounding mesh" << endl;
    scalar bigSpan = 10*meshControls().span();
    insertPoint(point2D(-bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point2D(-bigSpan, bigSpan), Vb::FAR_POINT);
    insertPoint(point2D(bigSpan, -bigSpan), Vb::FAR_POINT);
    insertPoint(point2D(bigSpan, bigSpan), Vb::FAR_POINT);
}


void Foam::CV2D::fast_restore_Delaunay(Vertex_handle vh)
{
    int i;
    Face_handle f = vh->face(), next, start(f);

    do
    {
        i=f->index(vh);
        if (!is_infinite(f))
        {
            if (!internal_flip(f, cw(i))) external_flip(f, i);
            if (f->neighbor(i) == start) start = f;
        }
        f = f->neighbor(cw(i));
    } while (f != start);
}


void Foam::CV2D::external_flip(Face_handle& f, int i)
{
    Face_handle n = f->neighbor(i);

    if
    (
        CGAL::ON_POSITIVE_SIDE
     != side_of_oriented_circle(n, f->vertex(i)->point())
    ) return;

    flip(f, i);
    i = n->index(f->vertex(i));
    external_flip(n, i);
}


bool Foam::CV2D::internal_flip(Face_handle& f, int i)
{
    Face_handle n = f->neighbor(i);

    if
    (
        CGAL::ON_POSITIVE_SIDE
     != side_of_oriented_circle(n, f->vertex(i)->point())
    )
    {
        return false;
    }

    flip(f, i);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CV2D::CV2D
(
    const Time& runTime,
    const dictionary& cvMeshDict
)
:
    Delaunay(),
    runTime_(runTime),
    rndGen_(64293*Pstream::myProcNo()),
    allGeometry_
    (
        IOobject
        (
            "cvSearchableSurfaces",
            runTime_.constant(),
            "triSurface",
            runTime_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        cvMeshDict.subDict("geometry"),
        cvMeshDict.lookupOrDefault("singleRegionName", true)
    ),
    qSurf_
    (
        runTime_,
        rndGen_,
        allGeometry_,
        cvMeshDict.subDict("surfaceConformation")
    ),
    controls_(cvMeshDict, qSurf_.globalBounds()),
    cellSizeControl_
    (
        runTime,
        cvMeshDict.subDict("motionControl").subDict("shapeControlFunctions"),
        qSurf_,
        controls_.minCellSize()
    ),
    relaxationModel_
    (
        relaxationModel::New
        (
            cvMeshDict.subDict("motionControl"),
            runTime
        )
    ),
    z_
    (
        point
        (
            cvMeshDict.subDict("surfaceConformation").lookup("locationInMesh")
        ).z()
    ),
    startOfInternalPoints_(0),
    startOfSurfacePointPairs_(0),
    startOfBoundaryConformPointPairs_(0),
    featurePoints_()
{
    Info<< meshControls() << endl;

    insertBoundingBox();
    insertFeaturePoints();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CV2D::~CV2D()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CV2D::insertPoints
(
    const point2DField& points,
    const scalar nearness
)
{
    Info<< "insertInitialPoints(const point2DField& points): ";

    startOfInternalPoints_ = number_of_vertices();
    label nVert = startOfInternalPoints_;

    // Add the points and index them
    forAll(points, i)
    {
        const point2D& p = points[i];

        if (qSurf_.wellInside(toPoint3D(p), nearness))
        {
            insert(toPoint(p))->index() = nVert++;
        }
        else
        {
            Warning
                << "Rejecting point " << p << " outside surface" << endl;
        }
    }

    Info<< nVert << " vertices inserted" << endl;

    if (meshControls().objOutput())
    {
        // Checking validity of triangulation
        assert(is_valid());

        writeTriangles("initial_triangles.obj", true);
        writeFaces("initial_faces.obj", true);
    }
}


void Foam::CV2D::insertPoints(const fileName& pointFileName)
{
    IFstream pointsFile(pointFileName);

    if (pointsFile.good())
    {
        insertPoints
        (
            point2DField(pointsFile),
            0.5*meshControls().minCellSize2()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Could not open pointsFile " << pointFileName
            << exit(FatalError);
    }
}


void Foam::CV2D::insertGrid()
{
    Info<< "insertInitialGrid: ";

    startOfInternalPoints_ = number_of_vertices();
    label nVert = startOfInternalPoints_;

    scalar x0 = qSurf_.globalBounds().min().x();
    scalar xR = qSurf_.globalBounds().max().x() - x0;
    int ni = int(xR/meshControls().minCellSize()) + 1;
    scalar deltax = xR/ni;

    scalar y0 = qSurf_.globalBounds().min().y();
    scalar yR = qSurf_.globalBounds().max().y() - y0;
    int nj = int(yR/meshControls().minCellSize()) + 1;
    scalar deltay = yR/nj;

    Random rndGen(1321);
    scalar pert = meshControls().randomPerturbation()*min(deltax, deltay);

    for (int i=0; i<ni; i++)
    {
        for (int j=0; j<nj; j++)
        {
            point p(x0 + i*deltax, y0 + j*deltay, 0);

            if (meshControls().randomiseInitialGrid())
            {
                p.x() += pert*(rndGen.scalar01() - 0.5);
                p.y() += pert*(rndGen.scalar01() - 0.5);
            }

            if (qSurf_.wellInside(p, 0.5*meshControls().minCellSize2()))
            {
                insert(Point(p.x(), p.y()))->index() = nVert++;
            }
        }
    }

    Info<< nVert << " vertices inserted" << endl;

    if (meshControls().objOutput())
    {
        // Checking validity of triangulation
        assert(is_valid());

        writeTriangles("initial_triangles.obj", true);
        writeFaces("initial_faces.obj", true);
    }
}


void Foam::CV2D::insertSurfacePointPairs()
{
    startOfSurfacePointPairs_ = number_of_vertices();

    if (meshControls().insertSurfaceNearestPointPairs())
    {
        insertSurfaceNearestPointPairs();
    }

    write("nearest");

    // Insertion of point-pairs for near-points may cause protrusions
    // so insertBoundaryConformPointPairs must be executed last
    if (meshControls().insertSurfaceNearPointPairs())
    {
        insertSurfaceNearPointPairs();
    }

    startOfBoundaryConformPointPairs_ = number_of_vertices();
}


void Foam::CV2D::boundaryConform()
{
    if (!meshControls().insertSurfaceNearestPointPairs())
    {
        markNearBoundaryPoints();
    }

    // Mark all the faces as SAVE_CHANGED
    for
    (
        Triangulation::Finite_faces_iterator fit = finite_faces_begin();
        fit != finite_faces_end();
        fit++
    )
    {
        fit->faceIndex() = Fb::SAVE_CHANGED;
    }

    for (label iter=1; iter<=meshControls().maxBoundaryConformingIter(); iter++)
    {
        label nIntersections = insertBoundaryConformPointPairs
        (
            "surfaceIntersections_" + Foam::name(iter) + ".obj"
        );

        if (nIntersections == 0)
        {
            break;
        }
        else
        {
            Info<< "BC iteration " << iter << ": "
                << nIntersections << " point-pairs inserted" << endl;
        }

        // Any faces changed by insertBoundaryConformPointPairs will now
        // be marked CHANGED, mark those as SAVE_CHANGED and those that
        // remained SAVE_CHANGED as UNCHANGED
        for
        (
            Triangulation::Finite_faces_iterator fit = finite_faces_begin();
            fit != finite_faces_end();
            fit++
        )
        {
            if (fit->faceIndex() == Fb::SAVE_CHANGED)
            {
                fit->faceIndex() = Fb::UNCHANGED;
            }
            else if (fit->faceIndex() == Fb::CHANGED)
            {
                fit->faceIndex() = Fb::SAVE_CHANGED;
            }
        }
    }

    Info<< nl;

    write("boundary");
}


void Foam::CV2D::removeSurfacePointPairs()
{
    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->index() >= startOfSurfacePointPairs_)
        {
            remove(vit);
        }
    }
}


void Foam::CV2D::newPoints()
{
    const scalar relaxation = relaxationModel_->relaxation();

    Info<< "Relaxation = " << relaxation << endl;

    Field<point2D> dualVertices(number_of_faces());

    label dualVerti = 0;

    // Find the dual point of each tetrahedron and assign it an index.
    for
    (
        Triangulation::Finite_faces_iterator fit = finite_faces_begin();
        fit != finite_faces_end();
        ++fit
    )
    {
        fit->faceIndex() = -1;

        if
        (
            fit->vertex(0)->internalOrBoundaryPoint()
         || fit->vertex(1)->internalOrBoundaryPoint()
         || fit->vertex(2)->internalOrBoundaryPoint()
        )
        {
            fit->faceIndex() = dualVerti;

            dualVertices[dualVerti] = toPoint2D(circumcenter(fit));

            dualVerti++;
        }
    }

    dualVertices.setSize(dualVerti);

    Field<vector2D> displacementAccumulator
    (
        startOfSurfacePointPairs_,
        vector2D::zero
    );

    // Calculate target size and alignment for vertices
    scalarField sizes
    (
        number_of_vertices(),
        meshControls().minCellSize()
    );

    Field<vector2D> alignments
    (
        number_of_vertices(),
        vector2D(1, 0)
    );

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            point2D vert = toPoint2D(vit->point());

            // alignment and size determination
            pointIndexHit pHit;
            label hitSurface = -1;

            qSurf_.findSurfaceNearest
            (
                toPoint3D(vert),
                meshControls().span2(),
                pHit,
                hitSurface
            );

            if (pHit.hit())
            {
                vectorField norm(1);
                allGeometry_[hitSurface].getNormal
                (
                    List<pointIndexHit>(1, pHit),
                    norm
                );

                alignments[vit->index()] = toPoint2D(norm[0]);

                sizes[vit->index()] =
                    cellSizeControl_.cellSize
                    (
                        toPoint3D(vit->point())
                    );
            }
        }
    }

    // Info<< "Calculated alignments" << endl;

    scalar cosAlignmentAcceptanceAngle = 0.68;

    // Upper and lower edge length ratios for weight
    scalar u = 1.0;
    scalar l = 0.7;

    PackedBoolList pointToBeRetained(startOfSurfacePointPairs_, true);

    std::list<Point> pointsToInsert;

    for
    (
        Triangulation::Finite_edges_iterator eit = finite_edges_begin();
        eit != finite_edges_end();
        eit++
    )
    {
        Vertex_handle vA = eit->first->vertex(cw(eit->second));
        Vertex_handle vB = eit->first->vertex(ccw(eit->second));

        if (!vA->internalOrBoundaryPoint() || !vB->internalOrBoundaryPoint())
        {
            continue;
        }

        const point2D& dualV1 = dualVertices[eit->first->faceIndex()];
        const point2D& dualV2 =
            dualVertices[eit->first->neighbor(eit->second)->faceIndex()];

        scalar dualEdgeLength = mag(dualV1 - dualV2);

        point2D dVA = toPoint2D(vA->point());
        point2D dVB = toPoint2D(vB->point());

        Field<vector2D> alignmentDirsA(2);

        alignmentDirsA[0] = alignments[vA->index()];
        alignmentDirsA[1] = vector2D
        (
           -alignmentDirsA[0].y(),
            alignmentDirsA[0].x()
        );

        Field<vector2D> alignmentDirsB(2);

        alignmentDirsB[0] = alignments[vB->index()];
        alignmentDirsB[1] = vector2D
        (
           -alignmentDirsB[0].y(),
            alignmentDirsB[0].x()
        );

        Field<vector2D> alignmentDirs(alignmentDirsA);

        forAll(alignmentDirsA, aA)
        {
            const vector2D& a(alignmentDirsA[aA]);

            scalar maxDotProduct = 0.0;

            forAll(alignmentDirsB, aB)
            {
                const vector2D& b(alignmentDirsB[aB]);

                scalar dotProduct = a & b;

                if (mag(dotProduct) > maxDotProduct)
                {
                    maxDotProduct = mag(dotProduct);

                    alignmentDirs[aA] = a + sign(dotProduct)*b;

                    alignmentDirs[aA] /= mag(alignmentDirs[aA]);
                }
            }
        }

        vector2D rAB = dVA - dVB;

        scalar rABMag = mag(rAB);

        forAll(alignmentDirs, aD)
        {
            vector2D& alignmentDir = alignmentDirs[aD];

            if ((rAB & alignmentDir) < 0)
            {
                // swap the direction of the alignment so that has the
                // same sense as rAB
                alignmentDir *= -1;
            }

            scalar alignmentDotProd = ((rAB/rABMag) & alignmentDir);

            if (alignmentDotProd > cosAlignmentAcceptanceAngle)
            {
                scalar targetFaceSize =
                    0.5*(sizes[vA->index()] + sizes[vB->index()]);

                // Test for changing aspect ratio on second alignment (first
                // alignment is neartest surface normal)
                // if (aD == 1)
                // {
                //     targetFaceSize *= 2.0;
                // }

                alignmentDir *= 0.5*targetFaceSize;

                vector2D delta = alignmentDir - 0.5*rAB;

                if (dualEdgeLength < 0.7*targetFaceSize)
                {
                    delta *= 0;
                }
                else if (dualEdgeLength < targetFaceSize)
                {
                    delta *=
                        (
                            dualEdgeLength
                           /(targetFaceSize*(u - l))
                          - 1/((u/l) - 1)
                        );
                }

                if
                (
                    vA->internalPoint()
                 && vB->internalPoint()
                 && rABMag > 1.75*targetFaceSize
                 && dualEdgeLength > 0.05*targetFaceSize
                 && alignmentDotProd > 0.93
                )
                {
                    // Point insertion
                    pointsToInsert.push_back(toPoint(0.5*(dVA + dVB)));
                }
                else if
                (
                    (vA->internalPoint() || vB->internalPoint())
                 && rABMag < 0.65*targetFaceSize
                )
                {
                    // Point removal

                    // Only insert a point at the midpoint of the short edge
                    // if neither attached point has already been identified
                    // to be removed.
                    if
                    (
                        pointToBeRetained[vA->index()] == true
                     && pointToBeRetained[vB->index()] == true
                    )
                    {
                        pointsToInsert.push_back(toPoint(0.5*(dVA + dVB)));
                    }

                    if (vA->internalPoint())
                    {
                        pointToBeRetained[vA->index()] = false;
                    }

                    if (vB->internalPoint())
                    {
                        pointToBeRetained[vB->index()] = false;
                    }
                }
                else
                {
                    if (vA->internalPoint())
                    {
                        displacementAccumulator[vA->index()] += delta;
                    }

                    if (vB->internalPoint())
                    {
                        displacementAccumulator[vB->index()] += -delta;
                    }
                }
            }
        }
    }

    vector2D totalDisp = sum(displacementAccumulator);
    scalar totalDist = sum(mag(displacementAccumulator));

    // Relax the calculated displacement
    displacementAccumulator *= relaxation;

    label numberOfNewPoints = pointsToInsert.size();

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint())
        {
            if (pointToBeRetained[vit->index()])
            {
                pointsToInsert.push_front
                (
                    toPoint
                    (
                        toPoint2D(vit->point())
                      + displacementAccumulator[vit->index()]
                    )
                );
            }
        }
    }

    // Clear the triangulation and reinsert the bounding box and feature points.
    // This is faster than removing and moving points.
    this->clear();

    insertBoundingBox();

    reinsertFeaturePoints();

    startOfInternalPoints_ = number_of_vertices();

    label nVert = startOfInternalPoints_;

    Info<< "Inserting " << numberOfNewPoints << " new points" << endl;

    // Use the range insert as it is faster than individually inserting points.
    insert(pointsToInsert.begin(), pointsToInsert.end());

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if
        (
            vit->type() == Vb::INTERNAL_POINT
         && vit->index() == Vb::INTERNAL_POINT
        )
        {
            vit->index() = nVert++;
        }
    }

    Info<< "    Total displacement = " << totalDisp << nl
        << "    Total distance = " << totalDist << nl
        << "    Points added = " << pointsToInsert.size()
        << endl;

    write("internal");

    insertSurfacePointPairs();

    boundaryConform();


    // Old Method
    /*
    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint())
        {
            // Current dual-cell defining vertex ("centre")
            point2DFromPoint defVert0 = toPoint2D(vit->point());

            Triangulation::Edge_circulator ec = incident_edges(vit);
            Triangulation::Edge_circulator ecStart = ec;

            // Circulate around the edges to find the first which is not
            // infinite
            do
            {
                if (!is_infinite(ec)) break;
            } while (++ec != ecStart);

            // Store the start-end of the first non-infinte edge
            point2D de0 = toPoint2D(circumcenter(ec->first));

            // Keep track of the maximum edge length^2
            scalar maxEdgeLen2 = 0.0;

            // Keep track of the index of the longest edge
            label edgecd0i = -1;

            // Edge counter
            label edgei = 0;

            do
            {
                if (!is_infinite(ec))
                {
                    // Get the end of the current edge
                    point2D de1 = toPoint2D
                    (
                        circumcenter(ec->first->neighbor(ec->second))
                    );

                    // Store the current edge vector
                    edges[edgei] = de1 - de0;

                    // Store the edge mid-point in the vertices array
                    vertices[edgei] = 0.5*(de1 + de0);

                    // Move the current edge end into the edge start for the
                    // next iteration
                    de0 = de1;

                    // Keep track of the longest edge

                    scalar edgeLen2 = magSqr(edges[edgei]);

                    if (edgeLen2 > maxEdgeLen2)
                    {
                        maxEdgeLen2 = edgeLen2;
                        edgecd0i = edgei;
                    }

                    edgei++;
                }
            } while (++ec != ecStart);

            // Initialise cd0 such that the mesh will align
            // in in the x-y directions
            vector2D cd0(1, 0);

            if (meshControls().relaxOrientation())
            {
                // Get the longest edge from the array and use as the primary
                // direction of the coordinate system of the "square" cell
                cd0 = edges[edgecd0i];
            }

            if (meshControls().nearWallAlignedDist() > 0)
            {
                pointIndexHit pHit = qSurf_.tree().findNearest
                (
                    toPoint3D(defVert0),
                    meshControls().nearWallAlignedDist2()
                );

                if (pHit.hit())
                {
                    cd0 = toPoint2D(faceNormals[pHit.index()]);
                }
            }

            // Rotate by 45deg needed to create an averaging procedure which
            // encourages the cells to be square
            cd0 = vector2D(cd0.x() + cd0.y(), cd0.y() - cd0.x());

            // Normalise the primary coordinate direction
            cd0 /= mag(cd0);

            // Calculate the orthogonal coordinate direction
            vector2D cd1(-cd0.y(), cd0.x());


            // Restart the circulator
            ec = ecStart;

            // ... and the counter
            edgei = 0;

            // Initialise the displacement for the centre and sum-weights
            vector2D disp = Zero;
            scalar sumw = 0;

            do
            {
                if (!is_infinite(ec))
                {
                    // Pick up the current edge
                    const vector2D& ei = edges[edgei];

                    // Calculate the centre to edge-centre vector
                    vector2D deltai = vertices[edgei] - defVert0;

                    // Set the weight for this edge contribution
                    scalar w = 1;

                    if (meshControls().squares())
                    {
                        w = magSqr(deltai.x()*ei.y() - deltai.y()*ei.x());
                        // alternative weights
                        // w = mag(deltai.x()*ei.y() - deltai.y()*ei.x());
                        // w = magSqr(ei)*mag(deltai);

                        // Use the following for an ~square mesh
                        // Find the coordinate contributions for this edge delta
                        scalar cd0deltai = cd0 & deltai;
                        scalar cd1deltai = cd1 & deltai;

                        // Create a "square" displacement
                        if (mag(cd0deltai) > mag(cd1deltai))
                        {
                            disp += (w*cd0deltai)*cd0;
                        }
                        else
                        {
                            disp += (w*cd1deltai)*cd1;
                        }
                    }
                    else
                    {
                        // Use this for a hexagon/pentagon mesh
                        disp += w*deltai;
                    }

                    // Sum the weights
                    sumw += w;
                }
                else
                {
                    FatalErrorInFunction
                        << "Infinite triangle found in internal mesh"
                        << exit(FatalError);
                }

                edgei++;

            } while (++ec != ecStart);

            // Calculate the average displacement
            disp /= sumw;
            totalDisp += disp;
            totalDist += mag(disp);

            // Move the point by a fraction of the average displacement
            movePoint(vit, defVert0 + relaxation*disp);
        }
    }

    Info << "\nTotal displacement = " << totalDisp
         << " total distance = " << totalDist << endl;
    */
}

/*
void Foam::CV2D::moveInternalPoints(const point2DField& newPoints)
{
    label pointi = 0;

    for
    (
        Triangulation::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint())
        {
            movePoint(vit, newPoints[pointi++]);
        }
    }
}
*/

void Foam::CV2D::write() const
{
    if (meshControls().objOutput())
    {
        writeFaces("allFaces.obj", false);
        writeFaces("faces.obj", true);
        writeTriangles("allTriangles.obj", false);
        writeTriangles("triangles.obj", true);
        writePatch("patch.pch");
    }
}


void Foam::CV2D::write(const word& stage) const
{
    if (meshControls().objOutput())
    {
        Foam::mkDir(stage + "Faces");
        Foam::mkDir(stage + "Triangles");

        writeFaces
        (
            stage
          + "Faces/allFaces_"
          + runTime_.timeName()
          + ".obj",
            false
        );

        writeTriangles
        (
            stage
          + "Triangles/allTriangles_"
          + runTime_.timeName()
          + ".obj",
            false
        );
    }
}


// ************************************************************************* //
