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
#include "initialPointsMethod.H"
#include "relaxationModel.H"
#include "faceAreaWeightModel.H"
#include "meshSearch.H"
#include "vectorTools.H"
#include "IOmanip.H"
#include "indexedCellChecks.H"
#include "controlMeshRefinement.H"
#include "smoothAlignmentSolver.H"
#include "OBJstream.H"
#include "indexedVertexOps.H"
#include "DelaunayMeshTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conformalVoronoiMesh, 0);

    template<>
    const char* NamedEnum
    <
        conformalVoronoiMesh::dualMeshPointType,
        5
    >::names[] =
    {
        "internal",
        "surface",
        "featureEdge",
        "featurePoint",
        "constrained"
    };
}

const Foam::NamedEnum<Foam::conformalVoronoiMesh::dualMeshPointType, 5>
    Foam::conformalVoronoiMesh::dualMeshPointTypeNames_;


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::cellSizeMeshOverlapsBackground() const
{
    const cellShapeControlMesh& cellSizeMesh =
        cellShapeControl_.shapeControlMesh();

    DynamicList<Foam::point> pts(number_of_vertices());

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalOrBoundaryPoint() && !vit->referred())
        {
            pts.append(topoint(vit->point()));
        }
    }

    boundBox bb(pts);

    boundBox cellSizeMeshBb = cellSizeMesh.bounds();

    bool fullyContained = true;

    if (!cellSizeMeshBb.contains(bb))
    {
        Pout<< "Triangulation not fully contained in cell size mesh."
            << endl;

        Pout<< "Cell Size Mesh Bounds = " << cellSizeMesh.bounds() << endl;
        Pout<< "foamyHexMesh Bounds         = " << bb << endl;

        fullyContained = false;
    }

    reduce(fullyContained, andOp<unsigned int>());

    Info<< "Triangulation is "
        << (fullyContained ? "fully" : "not fully")
        << " contained in the cell size mesh"
        << endl;
}


void Foam::conformalVoronoiMesh::insertInternalPoints
(
    List<Point>& points,
    bool distribute
)
{
    label nPoints = points.size();

    if (Pstream::parRun())
    {
        reduce(nPoints, sumOp<label>());
    }

    Info<< "    " << nPoints << " points to insert..." << endl;

    if (Pstream::parRun() && distribute)
    {
        List<Foam::point> transferPoints(points.size());

        forAll(points, pI)
        {
            transferPoints[pI] = topoint(points[pI]);
        }

        // Send the points that are not on this processor to the appropriate
        // place
        Foam::autoPtr<Foam::mapDistribute> map
        (
            decomposition_().distributePoints(transferPoints)
        );

        transferPoints.clear();

        map().distribute(points);
    }

    label nVert = number_of_vertices();

    insert(points.begin(), points.end());

    label nInserted(number_of_vertices() - nVert);

    if (Pstream::parRun())
    {
        reduce(nInserted, sumOp<label>());
    }

    Info<< "    " << nInserted << " points inserted"
        << ", failed to insert " << nPoints - nInserted
        << " ("
        << 100.0*(nPoints - nInserted)/(nInserted + small)
        << " %)"<< endl;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (CGAL::indexedVertexOps::uninitialised(vit))
        {
            vit->index() = getNewVertexIndex();
            vit->type() = Vb::vtInternal;
        }
    }
}


Foam::Map<Foam::label> Foam::conformalVoronoiMesh::insertPointPairs
(
    List<Vb>& vertices,
    bool distribute,
    bool reIndex
)
{
    if (Pstream::parRun() && distribute)
    {
        autoPtr<mapDistribute> mapDist =
            decomposition_().distributePoints(vertices);

        // Re-index the point pairs if one or both have been distributed.
        // If both, remove

        // If added a point, then need to know its point pair
        // If one moved, then update procIndex locally

        forAll(vertices, vI)
        {
            vertices[vI].procIndex() = Pstream::myProcNo();
        }
    }

    label preReinsertionSize(number_of_vertices());

    Map<label> oldToNewIndices =
        this->DelaunayMesh<Delaunay>::insertPoints(vertices, reIndex);

    const label nReinserted = returnReduce
    (
        label(number_of_vertices()) - preReinsertionSize,
        sumOp<label>()
    );

    Info<< "    Reinserted " << nReinserted << " vertices out of "
        << returnReduce(vertices.size(), sumOp<label>())
        << endl;

    return oldToNewIndices;
}


void Foam::conformalVoronoiMesh::insertSurfacePointPairs
(
    const pointIndexHitAndFeatureList& surfaceHits,
    const fileName fName,
    DynamicList<Vb>& pts
)
{
    forAll(surfaceHits, i)
    {
        vectorField norm(1);

        const pointIndexHit surfaceHit = surfaceHits[i].first();
        const label featureIndex = surfaceHits[i].second();

        allGeometry_[featureIndex].getNormal
        (
            List<pointIndexHit>(1, surfaceHit),
            norm
        );

        const vector& normal = norm[0];

        const Foam::point& surfacePt(surfaceHit.hitPoint());

        extendedFeatureEdgeMesh::sideVolumeType meshableSide =
            geometryToConformTo_.meshableSide(featureIndex, surfaceHit);

        if (meshableSide == extendedFeatureEdgeMesh::BOTH)
        {
            createBafflePointPair
            (
                pointPairDistance(surfacePt),
                surfacePt,
                normal,
                true,
                pts
            );
        }
        else if (meshableSide == extendedFeatureEdgeMesh::INSIDE)
        {
            createPointPair
            (
                pointPairDistance(surfacePt),
                surfacePt,
                normal,
                true,
                pts
            );
        }
        else if (meshableSide == extendedFeatureEdgeMesh::OUTSIDE)
        {
            createPointPair
            (
                pointPairDistance(surfacePt),
                surfacePt,
                -normal,
                true,
                pts
            );
        }
        else
        {
            WarningInFunction
                << meshableSide << ", bad"
                << endl;
        }
    }

    if (foamyHexMeshControls().objOutput() && fName != fileName::null)
    {
        DelaunayMeshTools::writeOBJ(time().path()/fName, pts);
    }
}


void Foam::conformalVoronoiMesh::insertEdgePointGroups
(
    const pointIndexHitAndFeatureList& edgeHits,
    const fileName fName,
    DynamicList<Vb>& pts
)
{
    forAll(edgeHits, i)
    {
        if (edgeHits[i].first().hit())
        {
            const extendedFeatureEdgeMesh& feMesh
            (
                geometryToConformTo_.features()[edgeHits[i].second()]
            );

//            const bool isBaffle =
//                geometryToConformTo_.isBaffleFeature(edgeHits[i].second());

            createEdgePointGroup
            (
                feMesh,
                edgeHits[i].first(),
                pts
            );
        }
    }

    if (foamyHexMeshControls().objOutput() && fName != fileName::null)
    {
        DelaunayMeshTools::writeOBJ(time().path()/fName, pts);
    }
}


bool Foam::conformalVoronoiMesh::nearFeaturePt(const Foam::point& pt) const
{
    scalar exclusionRangeSqr = featurePointExclusionDistanceSqr(pt);

    pointIndexHit info;
    label featureHit;

    geometryToConformTo_.findFeaturePointNearest
    (
        pt,
        exclusionRangeSqr,
        info,
        featureHit
    );

    return info.hit();
}


bool Foam::conformalVoronoiMesh::surfacePtNearFeatureEdge
(
    const Foam::point& pt
) const
{
    scalar exclusionRangeSqr = surfacePtExclusionDistanceSqr(pt);

    pointIndexHit info;
    label featureHit;

    geometryToConformTo_.findEdgeNearest
    (
        pt,
        exclusionRangeSqr,
        info,
        featureHit
    );

    return info.hit();
}


void Foam::conformalVoronoiMesh::insertInitialPoints()
{
    Info<< nl << "Inserting initial points" << endl;

    timeCheck("Before initial points call");

    List<Point> initPts = initialPointsMethod_->initialPoints();

    timeCheck("After initial points call");

    // Assume that the initial points method made the correct decision for
    // which processor each point should be on, so give distribute = false
    insertInternalPoints(initPts);

    if (initialPointsMethod_->fixInitialPoints())
    {
        for
        (
            Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            ++vit
        )
        {
            vit->fixed() = true;
        }
    }

    if (foamyHexMeshControls().objOutput())
    {
        DelaunayMeshTools::writeOBJ
        (
            time().path()/"initialPoints.obj",
            *this,
            Foam::indexedVertexEnum::vtInternal
        );
    }
}


void Foam::conformalVoronoiMesh::distribute()
{
    if (!Pstream::parRun())
    {
        return;
    }

    DynamicList<Foam::point> points(number_of_vertices());
    DynamicList<Foam::indexedVertexEnum::vertexType> types
    (
        number_of_vertices()
    );
    DynamicList<scalar> sizes(number_of_vertices());
    DynamicList<tensor> alignments(number_of_vertices());

    for
    (
        Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->real())
        {
            points.append(topoint(vit->point()));
            types.append(vit->type());
            sizes.append(vit->targetCellSize());
            alignments.append(vit->alignment());
        }
    }

    autoPtr<mapDistribute> mapDist =
        DistributedDelaunayMesh<Delaunay>::distribute(decomposition_(), points);

    mapDist().distribute(types);
    mapDist().distribute(sizes);
    mapDist().distribute(alignments);

    // Reset the entire tessellation
    DelaunayMesh<Delaunay>::reset();

    Info<< nl << "    Inserting distributed tessellation" << endl;

    // Internal points have to be inserted first

    DynamicList<Vb> verticesToInsert(points.size());

    forAll(points, pI)
    {
        verticesToInsert.append
        (
            Vb
            (
                toPoint(points[pI]),
                -1,
                types[pI],
                Pstream::myProcNo()
            )
        );

        verticesToInsert.last().targetCellSize() = sizes[pI];
        verticesToInsert.last().alignment() = alignments[pI];
    }

    this->rangeInsertWithInfo
    (
        verticesToInsert.begin(),
        verticesToInsert.end(),
        true
    );

    Info<< "    Total number of vertices after redistribution "
        << returnReduce
           (
               label(number_of_vertices()), sumOp<label>()
           )
        << endl;
}


void Foam::conformalVoronoiMesh::buildCellSizeAndAlignmentMesh()
{
    controlMeshRefinement meshRefinement
    (
        cellShapeControl_
    );

    smoothAlignmentSolver meshAlignmentSmoother
    (
        cellShapeControl_.shapeControlMesh()
    );

    meshRefinement.initialMeshPopulation(decomposition_);

    cellShapeControlMesh& cellSizeMesh = cellShapeControl_.shapeControlMesh();

    if (Pstream::parRun())
    {
        if (!distributeBackground(cellSizeMesh))
        {
            // Synchronise the cell size mesh if it has not been distributed
            cellSizeMesh.distribute(decomposition_);
        }
    }

    const dictionary& motionControlDict
        = foamyHexMeshControls().foamyHexMeshDict().subDict("motionControl");

    label nMaxIter = readLabel
    (
        motionControlDict.lookup("maxRefinementIterations")
    );

    Info<< "Maximum number of refinement iterations : " << nMaxIter << endl;

    for (label i = 0; i < nMaxIter; ++i)
    {
        label nAdded = meshRefinement.refineMesh(decomposition_);
        // cellShapeControl_.refineMesh(decomposition_);
        reduce(nAdded, sumOp<label>());

        if (Pstream::parRun())
        {
            cellSizeMesh.distribute(decomposition_);
        }

        Info<< "    Iteration " << i
            << " Added = " << nAdded << " points"
            << endl;

        if (nAdded == 0)
        {
            break;
        }
    }

    if (Pstream::parRun())
    {
        // Need to distribute the cell size mesh to cover the background mesh
        if (!distributeBackground(cellSizeMesh))
        {
            cellSizeMesh.distribute(decomposition_);
        }
    }

    label maxSmoothingIterations = readLabel
    (
        motionControlDict.lookup("maxSmoothingIterations")
    );
    meshAlignmentSmoother.smoothAlignments(maxSmoothingIterations);

    Info<< "Background cell size and alignment mesh:" << endl;
    cellSizeMesh.printInfo(Info);

    Info<< "Triangulation is "
        << (cellSizeMesh.is_valid() ? "valid" : "not valid!" )
        << endl;

    if (foamyHexMeshControls().writeCellShapeControlMesh())
    {
        // cellSizeMesh.writeTriangulation();
        cellSizeMesh.write();
    }

    if (foamyHexMeshControls().printVertexInfo())
    {
        cellSizeMesh.printVertexInfo(Info);
    }

//    Info<< "Estimated number of cells in final mesh = "
//        << returnReduce
//           (
//               cellSizeMesh.estimateCellCount(decomposition_),
//               sumOp<label>()
//           )
//        << endl;
}


void Foam::conformalVoronoiMesh::setVertexSizeAndAlignment()
{
    Info<< nl << "Calculating target cell alignment and size" << endl;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        vit++
    )
    {
        if (vit->internalOrBoundaryPoint())
        {
            pointFromPoint pt = topoint(vit->point());

            cellShapeControls().cellSizeAndAlignment
            (
                pt,
                vit->targetCellSize(),
                vit->alignment()
            );
        }
    }
}


Foam::face Foam::conformalVoronoiMesh::buildDualFace
(
    const Delaunay::Finite_edges_iterator& eit
) const
{
    Cell_circulator ccStart = incident_cells(*eit);
    Cell_circulator cc1 = ccStart;
    Cell_circulator cc2 = cc1;

    // Advance the second circulator so that it always stays on the next
    // cell around the edge;
    cc2++;

    DynamicList<label> verticesOnFace;

    label nUniqueVertices = 0;

    do
    {
        if
        (
            cc1->hasFarPoint() || cc2->hasFarPoint()
         || is_infinite(cc1) || is_infinite(cc2)
        )
        {
            Cell_handle c = eit->first;
            Vertex_handle vA = c->vertex(eit->second);
            Vertex_handle vB = c->vertex(eit->third);

//            DelaunayMeshTools::drawDelaunayCell(Pout, cc1);
//            DelaunayMeshTools::drawDelaunayCell(Pout, cc2);

            WarningInFunction
                << "Dual face uses circumcenter defined by a "
                << "Delaunay tetrahedron with no internal "
                << "or boundary points.  Defining Delaunay edge ends: "
                << vA->info() << " "
                << vB->info() << nl
                <<endl;//<< exit(FatalError);
        }
        else
        {
            label cc1I = cc1->cellIndex();
            label cc2I = cc2->cellIndex();

            if (cc1I != cc2I)
            {
                if (findIndex(verticesOnFace, cc1I) == -1)
                {
                    nUniqueVertices++;
                }

                verticesOnFace.append(cc1I);
            }
        }

        cc1++;
        cc2++;

    } while (cc1 != ccStart);

    verticesOnFace.shrink();

    if (verticesOnFace.size() >= 3 && nUniqueVertices < 3)
    {
        // There are not enough unique vertices on this face to
        // justify its size, it may have a form like:

        // Vertices:
        // A                                  B
        // A                                  B

        // Face:
        // ABAB

        // Setting the size to be below 3, so that it will not be
        // created

        verticesOnFace.setSize(nUniqueVertices);
    }

    return face(verticesOnFace);
}


Foam::label Foam::conformalVoronoiMesh::maxFilterCount
(
    const Delaunay::Finite_edges_iterator& eit
) const
{
    Cell_circulator ccStart = incident_cells(*eit);
    Cell_circulator cc = ccStart;

    label maxFC = 0;

    do
    {
        if (cc->hasFarPoint())
        {
            Cell_handle c = eit->first;
            Vertex_handle vA = c->vertex(eit->second);
            Vertex_handle vB = c->vertex(eit->third);

            FatalErrorInFunction
                << "Dual face uses circumcenter defined by a "
                << "Delaunay tetrahedron with no internal "
                << "or boundary points.  Defining Delaunay edge ends: "
                << topoint(vA->point()) << " "
                << topoint(vB->point()) << nl
                << exit(FatalError);
        }

        if (cc->filterCount() > maxFC)
        {
            maxFC = cc->filterCount();
        }

        cc++;

    } while (cc != ccStart);

    return maxFC;
}


bool Foam::conformalVoronoiMesh::ownerAndNeighbour
(
    Vertex_handle vA,
    Vertex_handle vB,
    label& owner,
    label& neighbour
) const
{
    bool reverse = false;

    owner = -1;
    neighbour = -1;

    label dualCellIndexA = vA->index();

    if (!vA->internalOrBoundaryPoint() || vA->referred())
    {
        if (!vA->constrained())
        {
            dualCellIndexA = -1;
        }
    }

    label dualCellIndexB = vB->index();

    if (!vB->internalOrBoundaryPoint() || vB->referred())
    {
        if (!vB->constrained())
        {
            dualCellIndexB = -1;
        }
    }

    if (dualCellIndexA == -1 && dualCellIndexB == -1)
    {
        FatalErrorInFunction
            << "Attempting to create a face joining "
            << "two unindexed dual cells "
            << exit(FatalError);
    }
    else if (dualCellIndexA == -1 || dualCellIndexB == -1)
    {
        // boundary face, find which is the owner

        if (dualCellIndexA == -1)
        {
            owner = dualCellIndexB;

            reverse = true;
        }
        else
        {
            owner = dualCellIndexA;
        }
    }
    else
    {
        // internal face, find the lower cell to be the owner

        if (dualCellIndexB > dualCellIndexA)
        {
            owner = dualCellIndexA;
            neighbour = dualCellIndexB;
        }
        else
        {
            owner = dualCellIndexB;
            neighbour = dualCellIndexA;

            // reverse face order to correctly orientate normal
            reverse = true;
        }
    }

    return reverse;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conformalVoronoiMesh::conformalVoronoiMesh
(
    const Time& runTime,
    const dictionary& foamyHexMeshDict
)
:
    DistributedDelaunayMesh<Delaunay>(runTime),
    runTime_(runTime),
    rndGen_(64293*Pstream::myProcNo()),
    foamyHexMeshControls_(foamyHexMeshDict),
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
        foamyHexMeshDict.subDict("geometry"),
        foamyHexMeshDict.lookupOrDefault("singleRegionName", true)
    ),
    geometryToConformTo_
    (
        runTime_,
        rndGen_,
        allGeometry_,
        foamyHexMeshDict.subDict("surfaceConformation")
    ),
    decomposition_
    (
        Pstream::parRun()
      ? new backgroundMeshDecomposition
        (
            runTime_,
            rndGen_,
            geometryToConformTo_,
            foamyHexMeshControls().foamyHexMeshDict().subDict
            (
                "backgroundMeshDecomposition"
            )
        )
      : nullptr
    ),
    cellShapeControl_
    (
        runTime_,
        foamyHexMeshControls_,
        allGeometry_,
        geometryToConformTo_
    ),
    limitBounds_(),
    ptPairs_(*this),
    ftPtConformer_(*this),
    edgeLocationTreePtr_(),
    surfacePtLocationTreePtr_(),
    surfaceConformationVertices_(),
    initialPointsMethod_
    (
        initialPointsMethod::New
        (
            foamyHexMeshDict.subDict("initialPoints"),
            runTime_,
            rndGen_,
            geometryToConformTo_,
            cellShapeControl_,
            decomposition_
        )
    ),
    relaxationModel_
    (
        relaxationModel::New
        (
            foamyHexMeshDict.subDict("motionControl"),
            runTime_
        )
    ),
    faceAreaWeightModel_
    (
        faceAreaWeightModel::New
        (
            foamyHexMeshDict.subDict("motionControl")
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conformalVoronoiMesh::~conformalVoronoiMesh()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::conformalVoronoiMesh::initialiseForMotion()
{
    if (foamyHexMeshControls().objOutput())
    {
        geometryToConformTo_.writeFeatureObj("foamyHexMesh");
    }

    buildCellSizeAndAlignmentMesh();

    insertInitialPoints();

    insertFeaturePoints(true);

    setVertexSizeAndAlignment();

    cellSizeMeshOverlapsBackground();

    // Improve the guess that the backgroundMeshDecomposition makes with the
    // initial positions.  Use before building the surface conformation to
    // better balance the surface conformation load.
    distributeBackground(*this);

    buildSurfaceConformation();

    // The introduction of the surface conformation may have distorted the
    // balance of vertices, distribute if necessary.
    distributeBackground(*this);

    if (Pstream::parRun())
    {
        sync(decomposition_().procBounds());
    }

    // Do not store the surface conformation until after it has been
    // (potentially) redistributed.
    storeSurfaceConformation();

    // Report any Delaunay vertices that do not think that they are in the
    // domain the processor they are on.
    // reportProcessorOccupancy();

    cellSizeMeshOverlapsBackground();

    if (foamyHexMeshControls().printVertexInfo())
    {
        printVertexInfo(Info);
    }

    if (foamyHexMeshControls().objOutput())
    {
        DelaunayMeshTools::writeOBJ
        (
            time().path()/"internalPoints_" + time().timeName() + ".obj",
            *this,
            Foam::indexedVertexEnum::vtUnassigned,
            Foam::indexedVertexEnum::vtExternalFeaturePoint
        );
    }
}


void Foam::conformalVoronoiMesh::initialiseForConformation()
{
    if (Pstream::parRun())
    {
        decomposition_.reset
        (
            new backgroundMeshDecomposition
            (
                runTime_,
                rndGen_,
                geometryToConformTo_,
                foamyHexMeshControls().foamyHexMeshDict().subDict
                (
                    "backgroundMeshDecomposition"
                )
            )
        );
    }

    insertInitialPoints();

    insertFeaturePoints();

    // Improve the guess that the backgroundMeshDecomposition makes with the
    // initial positions.  Use before building the surface conformation to
    // better balance the surface conformation load.
    distributeBackground(*this);

    buildSurfaceConformation();

    // The introduction of the surface conformation may have distorted the
    // balance of vertices, distribute if necessary.
    distributeBackground(*this);

    if (Pstream::parRun())
    {
        sync(decomposition_().procBounds());
    }

    cellSizeMeshOverlapsBackground();

    if (foamyHexMeshControls().printVertexInfo())
    {
        printVertexInfo(Info);
    }
}


void Foam::conformalVoronoiMesh::move()
{
    timeCheck("Start of move");

    scalar relaxation = relaxationModel_->relaxation();

    Info<< nl << "Relaxation = " << relaxation << endl;

    pointField dualVertices(number_of_finite_cells());

    this->resetCellCount();

    // Find the dual point of each tetrahedron and assign it an index.
    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        cit->cellIndex() = Cb::ctUnassigned;

        if (cit->anyInternalOrBoundaryDualVertex())
        {
            cit->cellIndex() = getNewCellIndex();

            dualVertices[cit->cellIndex()] = cit->dual();
        }

        if (cit->hasFarPoint())
        {
            cit->cellIndex() = Cb::ctFar;
        }
    }

    dualVertices.setSize(cellCount());

    setVertexSizeAndAlignment();

    timeCheck("Determined sizes and alignments");

    Info<< nl << "Determining vertex displacements" << endl;

    vectorField cartesianDirections(3);

    cartesianDirections[0] = vector(1, 0, 0);
    cartesianDirections[1] = vector(0, 1, 0);
    cartesianDirections[2] = vector(0, 0, 1);

    vectorField displacementAccumulator
    (
        number_of_vertices(),
        Zero
    );

    PackedBoolList pointToBeRetained
    (
        number_of_vertices(),
        true
    );

    DynamicList<Point> pointsToInsert(number_of_vertices());

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

        if
        (
            (
                vA->internalPoint() && !vA->referred()
             && vB->internalOrBoundaryPoint()
            )
         || (
                vB->internalPoint() && !vB->referred()
             && vA->internalOrBoundaryPoint()
            )
        )
        {
            pointFromPoint dVA = topoint(vA->point());
            pointFromPoint dVB = topoint(vB->point());

            Field<vector> alignmentDirsA
            (
                vA->alignment().T() & cartesianDirections
            );
            Field<vector> alignmentDirsB
            (
                vB->alignment().T() & cartesianDirections
            );

            Field<vector> alignmentDirs(alignmentDirsA);

            forAll(alignmentDirsA, aA)
            {
                const vector& a = alignmentDirsA[aA];

                scalar maxDotProduct = 0.0;

                forAll(alignmentDirsB, aB)
                {
                    const vector& b = alignmentDirsB[aB];

                    const scalar dotProduct = a & b;

                    if (mag(dotProduct) > maxDotProduct)
                    {
                        maxDotProduct = mag(dotProduct);

                        alignmentDirs[aA] = a + sign(dotProduct)*b;

                        alignmentDirs[aA] /= mag(alignmentDirs[aA]);
                    }
                }
            }

            vector rAB = dVA - dVB;

            scalar rABMag = mag(rAB);

            if (rABMag < small)
            {
                // Removal of close points

                if
                (
                    vA->internalPoint() && !vA->referred() && !vA->fixed()
                 && vB->internalPoint() && !vB->referred() && !vB->fixed()
                )
                {
                    // Only insert a point at the midpoint of
                    // the short edge if neither attached
                    // point has already been identified to be
                    // removed.

                    if
                    (
                        pointToBeRetained[vA->index()] == true
                     && pointToBeRetained[vB->index()] == true
                    )
                    {
                        const Foam::point pt(0.5*(dVA + dVB));

                        if (internalPointIsInside(pt))
                        {
                            pointsToInsert.append(toPoint(pt));
                        }
                    }
                }

                if (vA->internalPoint() && !vA->referred() && !vA->fixed())
                {
                    pointToBeRetained[vA->index()] = false;
                }

                if (vB->internalPoint() && !vB->referred() && !vB->fixed())
                {
                    pointToBeRetained[vB->index()] = false;
                }

                // Do not consider this Delaunay edge any further

                continue;
            }

            forAll(alignmentDirs, aD)
            {
                vector& alignmentDir = alignmentDirs[aD];

                scalar dotProd = rAB & alignmentDir;

                if (dotProd < 0)
                {
                    // swap the direction of the alignment so that has the
                    // same sense as rAB
                    alignmentDir *= -1;
                    dotProd *= -1;
                }

                const scalar alignmentDotProd = dotProd/rABMag;

                if
                (
                    alignmentDotProd
                  > foamyHexMeshControls().cosAlignmentAcceptanceAngle()
                )
                {
                    scalar targetCellSize =
                        CGAL::indexedVertexOps::averageCellSize(vA, vB);

                    scalar targetFaceArea = sqr(targetCellSize);

                    const vector originalAlignmentDir = alignmentDir;

                    // Update cell size and face area
                    cellShapeControls().aspectRatio().updateCellSizeAndFaceArea
                    (
                        alignmentDir,
                        targetFaceArea,
                        targetCellSize
                    );

                    // Vector to move end points around middle of vector
                    // to align edge (i.e. dual face normal) with alignment
                    // directions.
                    vector delta = alignmentDir - 0.5*rAB;

                    face dualFace = buildDualFace(eit);

//                    Pout<< dualFace << endl;
//                    Pout<< "    " << vA->info() << endl;
//                    Pout<< "    " << vB->info() << endl;

                    const scalar faceArea = dualFace.mag(dualVertices);

                    // Update delta vector
                    cellShapeControls().aspectRatio().updateDeltaVector
                    (
                        originalAlignmentDir,
                        targetCellSize,
                        rABMag,
                        delta
                    );

//                    if (targetFaceArea == 0)
//                    {
//                        Pout<< vA->info() << vB->info();
//
//                        Cell_handle ch = locate(vA->point());
//                        if (is_infinite(ch))
//                        {
//                            Pout<< "vA " << vA->targetCellSize() << endl;
//                        }
//
//                        ch = locate(vB->point());
//                        if (is_infinite(ch))
//                        {
//                            Pout<< "vB " << vB->targetCellSize() << endl;
//                        }
//                    }

                    delta *= faceAreaWeightModel_->faceAreaWeight
                    (
                        faceArea/targetFaceArea
                    );

                    if
                    (
                        (
                            (vA->internalPoint() && vB->internalPoint())
                         && (!vA->referred() || !vB->referred())
//                         ||
//                            (
//                                vA->referredInternalPoint()
//                             && vB->referredInternalPoint()
//                            )
                        )
                     && rABMag
                      > foamyHexMeshControls().insertionDistCoeff()
                       *targetCellSize
                     && faceArea
                      > foamyHexMeshControls().faceAreaRatioCoeff()
                       *targetFaceArea
                     && alignmentDotProd
                      > foamyHexMeshControls().cosInsertionAcceptanceAngle()
                    )
                    {
                        // Point insertion
                        if
                        (
                            !geometryToConformTo_.findSurfaceAnyIntersection
                            (
                                dVA,
                                dVB
                            )
                        )
                        {
                            const Foam::point newPt(0.5*(dVA + dVB));

                            // Prevent insertions spanning surfaces
                            if (internalPointIsInside(newPt))
                            {
                                if (Pstream::parRun())
                                {
                                    if
                                    (
                                        decomposition().positionOnThisProcessor
                                        (
                                            newPt
                                        )
                                    )
                                    {
                                        pointsToInsert.append(toPoint(newPt));
                                    }
                                }
                                else
                                {
                                    pointsToInsert.append(toPoint(newPt));
                                }
                            }
                        }
                    }
                    else if
                    (
                        (
                            (vA->internalPoint() && !vA->referred())
                         || (vB->internalPoint() && !vB->referred())
                        )
                     && rABMag
                      < foamyHexMeshControls().removalDistCoeff()
                       *targetCellSize
                    )
                    {
                        // Point removal
                        if
                        (
                            (
                                vA->internalPoint()
                             && !vA->referred()
                             && !vA->fixed()
                            )
                         &&
                            (
                                vB->internalPoint()
                             && !vB->referred()
                             && !vB->fixed()
                            )
                        )
                        {
                            // Only insert a point at the midpoint of
                            // the short edge if neither attached
                            // point has already been identified to be
                            // removed.
                            if
                            (
                                pointToBeRetained[vA->index()] == true
                             && pointToBeRetained[vB->index()] == true
                            )
                            {
                                const Foam::point pt(0.5*(dVA + dVB));

                                if (internalPointIsInside(pt))
                                {
                                    pointsToInsert.append(toPoint(pt));
                                }
                            }
                        }

                        if
                        (
                            vA->internalPoint()
                         && !vA->referred()
                         && !vA->fixed()
                        )
                        {
                            pointToBeRetained[vA->index()] = false;
                        }

                        if
                        (
                            vB->internalPoint()
                         && !vB->referred()
                         && !vB->fixed()
                        )
                        {
                            pointToBeRetained[vB->index()] = false;
                        }
                    }
                    else
                    {
                        if
                        (
                            vA->internalPoint()
                         && !vA->referred()
                         && !vA->fixed()
                        )
                        {
                            if (vB->fixed())
                            {
                                displacementAccumulator[vA->index()] += 2*delta;
                            }
                            else
                            {
                                displacementAccumulator[vA->index()] += delta;
                            }
                        }

                        if
                        (
                            vB->internalPoint()
                         && !vB->referred()
                         && !vB->fixed()
                        )
                        {
                            if (vA->fixed())
                            {
                                displacementAccumulator[vB->index()] -= 2*delta;
                            }
                            else
                            {
                                displacementAccumulator[vB->index()] -= delta;
                            }
                        }
                    }
                }
            }
        }
    }

    Info<< "Limit displacements" << endl;

    // Limit displacements that pierce, or get too close to the surface
    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint() && !vit->referred() && !vit->fixed())
        {
            if (pointToBeRetained[vit->index()] == true)
            {
                limitDisplacement
                (
                    vit,
                    displacementAccumulator[vit->index()]
                );
            }
        }
    }

    vector totalDisp = gSum(displacementAccumulator);
    scalar totalDist = gSum(mag(displacementAccumulator));

    displacementAccumulator *= relaxation;

    Info<< "Sum displacements" << endl;

    label nPointsToRetain = 0;
    label nPointsToRemove = 0;

    for
    (
        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
        vit != finite_vertices_end();
        ++vit
    )
    {
        if (vit->internalPoint() && !vit->referred() && !vit->fixed())
        {
            if (pointToBeRetained[vit->index()] == true)
            {
                // Convert vit->point() to FOAM vector (double) to do addition,
                // avoids memory increase because a record of the constructions
                // would be kept otherwise.
                // See cgal-discuss@lists-sop.inria.fr:
                // "Memory issue with openSUSE 11.3, exact kernel, adding
                //  points/vectors"
                // 14/1/2011.
                // Only necessary if using an exact constructions kernel
                // (extended precision)
                Foam::point pt
                (
                    topoint(vit->point())
                  + displacementAccumulator[vit->index()]
                );

                if (internalPointIsInside(pt))
                {
                    pointsToInsert.append(toPoint(pt));
                    nPointsToRemove++;
                }

                nPointsToRetain++;
            }
        }
    }

    pointsToInsert.shrink();

    Info<< returnReduce
           (
               nPointsToRetain - nPointsToRemove,
               sumOp<label>()
           )
        << " internal points are outside the domain. "
        << "They will not be inserted." << endl;

    // Save displacements to file.
    if (foamyHexMeshControls().objOutput() && time().writeTime())
    {
        Info<< "Writing point displacement vectors to file." << endl;
        OFstream str
        (
            time().path()/"displacements_" + runTime_.timeName() + ".obj"
        );

        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            ++vit
        )
        {
            if (vit->internalPoint() && !vit->referred())
            {
                if (pointToBeRetained[vit->index()] == true)
                {
                    meshTools::writeOBJ(str, topoint(vit->point()));

                    str << "vn "
                        << displacementAccumulator[vit->index()][0] << " "
                        << displacementAccumulator[vit->index()][1] << " "
                        << displacementAccumulator[vit->index()][2] << " "
                        << endl;
                }
            }
        }
    }

    // Remove the entire tessellation
    DelaunayMesh<Delaunay>::reset();

    timeCheck("Displacement calculated");

    Info<< nl<< "Inserting displaced tessellation" << endl;

    insertFeaturePoints(true);

    insertInternalPoints(pointsToInsert, true);

    // Fix points that have not been significantly displaced
//    for
//    (
//        Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
//        vit != finite_vertices_end();
//        ++vit
//    )
//    {
//        if (vit->internalPoint())
//        {
//            if
//            (
//                mag(displacementAccumulator[vit->index()])
//              < 0.1*targetCellSize(topoint(vit->point()))
//            )
//            {
//                vit->setVertexFixed();
//            }
//        }
//    }

    timeCheck("Internal points inserted");

    {
        // Check that no index is shared between any of the local points
        labelHashSet usedIndices;
        for
        (
            Delaunay::Finite_vertices_iterator vit = finite_vertices_begin();
            vit != finite_vertices_end();
            ++vit
        )
        {
            if (!vit->referred() && !usedIndices.insert(vit->index()))
            {
                FatalErrorInFunction
                    << "Index already used! Could not insert: " << nl
                    << vit->info()
                    << abort(FatalError);
            }
        }
    }

    conformToSurface();

    if (foamyHexMeshControls().objOutput())
    {
        DelaunayMeshTools::writeOBJ
        (
            time().path()/"internalPoints_" + time().timeName() + ".obj",
            *this,
            Foam::indexedVertexEnum::vtInternal
        );

        if (reconformToSurface())
        {
            DelaunayMeshTools::writeBoundaryPoints
            (
                time().path()/"boundaryPoints_" + time().timeName() + ".obj",
                *this
            );

            DelaunayMeshTools::writeOBJ
            (
                time().path()/"internalBoundaryPoints_" + time().timeName()
              + ".obj",
                *this,
                Foam::indexedVertexEnum::vtInternalSurface,
                Foam::indexedVertexEnum::vtInternalFeaturePoint
            );

            DelaunayMeshTools::writeOBJ
            (
                time().path()/"externalBoundaryPoints_" + time().timeName()
              + ".obj",
                *this,
                Foam::indexedVertexEnum::vtExternalSurface,
                Foam::indexedVertexEnum::vtExternalFeaturePoint
            );

            OBJstream multipleIntersections
            (
                "multipleIntersections_"
              + time().timeName()
              + ".obj"
            );

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

                Foam::point ptA(topoint(vA->point()));
                Foam::point ptB(topoint(vB->point()));

                List<pointIndexHit> surfHits;
                labelList hitSurfaces;

                geometryToConformTo().findSurfaceAllIntersections
                (
                    ptA,
                    ptB,
                    surfHits,
                    hitSurfaces
                );

                if
                (
                    surfHits.size() >= 2
                 || (
                     surfHits.size() == 0
                  && (
                          (vA->externalBoundaryPoint()
                       && vB->internalBoundaryPoint())
                       || (vB->externalBoundaryPoint()
                       && vA->internalBoundaryPoint())
                     )
                    )
                )
                {
                    multipleIntersections.write(linePointRef(ptA, ptB));
                }
            }
        }
    }

    timeCheck("After conformToSurface");

    if (foamyHexMeshControls().printVertexInfo())
    {
        printVertexInfo(Info);
    }

    if (time().writeTime())
    {
        writeMesh(time().timeName());
    }

    Info<< nl
        << "Total displacement = " << totalDisp << nl
        << "Total distance = " << totalDist << nl
        << endl;
}


void Foam::conformalVoronoiMesh::checkCoPlanarCells() const
{
    typedef CGAL::Exact_predicates_exact_constructions_kernel   Kexact;
    typedef CGAL::Point_3<Kexact>                               PointExact;

    if (!is_valid())
    {
        Pout<< "Triangulation is invalid!" << endl;
    }

    OFstream str("badCells.obj");

    label badCells = 0;

    for
    (
        Delaunay::Finite_cells_iterator cit = finite_cells_begin();
        cit != finite_cells_end();
        ++cit
    )
    {
        const scalar quality = foamyHexMeshChecks::coplanarTet(cit, 1e-16);

        if (quality == 0)
        {
            Pout<< "COPLANAR: " << cit->info() << nl
                << "    quality = " << quality << nl
                << "    dual    = " << topoint(cit->dual()) << endl;

            DelaunayMeshTools::drawDelaunayCell(str, cit, badCells++);

            FixedList<PointExact, 4> cellVerticesExact(PointExact(0,0,0));
            forAll(cellVerticesExact, vI)
            {
                cellVerticesExact[vI] = PointExact
                (
                    cit->vertex(vI)->point().x(),
                    cit->vertex(vI)->point().y(),
                    cit->vertex(vI)->point().z()
                );
            }

            PointExact synchronisedDual = CGAL::circumcenter<Kexact>
            (
                cellVerticesExact[0],
                cellVerticesExact[1],
                cellVerticesExact[2],
                cellVerticesExact[3]
            );

            Foam::point exactPt
            (
                CGAL::to_double(synchronisedDual.x()),
                CGAL::to_double(synchronisedDual.y()),
                CGAL::to_double(synchronisedDual.z())
            );

            Info<< "inexact = " << cit->dual() << nl
                << "exact   = " << exactPt << endl;
        }
    }

    Pout<< "There are " << badCells << " bad cells out of "
        << number_of_finite_cells() << endl;


    label nNonGabriel = 0;
    for
    (
        Delaunay::Finite_facets_iterator fit = finite_facets_begin();
        fit != finite_facets_end();
        ++fit
    )
    {
        if (!is_Gabriel(*fit))
        {
            nNonGabriel++;//Pout<< "Non-gabriel face" << endl;
        }
    }

    Pout<< "There are " << nNonGabriel << " non-Gabriel faces out of "
        << number_of_finite_facets() << endl;
}


// ************************************************************************* //
