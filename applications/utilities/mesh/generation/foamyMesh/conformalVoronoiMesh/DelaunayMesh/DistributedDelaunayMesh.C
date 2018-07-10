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

#include "DistributedDelaunayMesh.H"
#include "meshSearch.H"
#include "mapDistribute.H"
#include "zeroGradientFvPatchFields.H"
#include "pointConversion.H"
#include "indexedVertexEnum.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

template<class Triangulation>
Foam::autoPtr<Foam::mapDistribute>
Foam::DistributedDelaunayMesh<Triangulation>::buildMap
(
    const List<label>& toProc
)
{
    // Determine send map
    // ~~~~~~~~~~~~~~~~~~

    // 1. Count
    labelList nSend(Pstream::nProcs(), 0);

    forAll(toProc, i)
    {
        label proci = toProc[i];

        nSend[proci]++;
    }


    // 2. Size sendMap
    labelListList sendMap(Pstream::nProcs());

    forAll(nSend, proci)
    {
        sendMap[proci].setSize(nSend[proci]);

        nSend[proci] = 0;
    }

    // 3. Fill sendMap
    forAll(toProc, i)
    {
        label proci = toProc[i];

        sendMap[proci][nSend[proci]++] = i;
    }

    // 4. Send over how many I need to receive
    labelList recvSizes;
    Pstream::exchangeSizes(sendMap, recvSizes);


    // Determine receive map
    // ~~~~~~~~~~~~~~~~~~~~~

    labelListList constructMap(Pstream::nProcs());

    // Local transfers first
    constructMap[Pstream::myProcNo()] = identity
    (
        sendMap[Pstream::myProcNo()].size()
    );

    label constructSize = constructMap[Pstream::myProcNo()].size();

    forAll(constructMap, proci)
    {
        if (proci != Pstream::myProcNo())
        {
            label nRecv = recvSizes[proci];

            constructMap[proci].setSize(nRecv);

            for (label i = 0; i < nRecv; i++)
            {
                constructMap[proci][i] = constructSize++;
            }
        }
    }

    return autoPtr<mapDistribute>
    (
        new mapDistribute
        (
            constructSize,
            sendMap.xfer(),
            constructMap.xfer()
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Triangulation>
Foam::DistributedDelaunayMesh<Triangulation>::DistributedDelaunayMesh
(
    const Time& runTime
)
:
    DelaunayMesh<Triangulation>(runTime),
    allBackgroundMeshBounds_()
{}


template<class Triangulation>
Foam::DistributedDelaunayMesh<Triangulation>::DistributedDelaunayMesh
(
    const Time& runTime,
    const word& meshName
)
:
    DelaunayMesh<Triangulation>(runTime, meshName),
    allBackgroundMeshBounds_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Triangulation>
Foam::DistributedDelaunayMesh<Triangulation>::~DistributedDelaunayMesh()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Triangulation>
bool Foam::DistributedDelaunayMesh<Triangulation>::distributeBoundBoxes
(
    const boundBox& bb
)
{
    allBackgroundMeshBounds_.reset(new List<boundBox>(Pstream::nProcs()));

    // Give the bounds of every processor to every other processor
    allBackgroundMeshBounds_()[Pstream::myProcNo()] = bb;

    Pstream::gatherList(allBackgroundMeshBounds_());
    Pstream::scatterList(allBackgroundMeshBounds_());

    return true;
}


template<class Triangulation>
bool Foam::DistributedDelaunayMesh<Triangulation>::isLocal
(
    const Vertex_handle& v
) const
{
    return isLocal(v->procIndex());
}


template<class Triangulation>
bool Foam::DistributedDelaunayMesh<Triangulation>::isLocal
(
    const label localProcIndex
) const
{
    return localProcIndex == Pstream::myProcNo();
}


template<class Triangulation>
Foam::labelList Foam::DistributedDelaunayMesh<Triangulation>::overlapProcessors
(
    const point& centre,
    const scalar radiusSqr
) const
{
    DynamicList<label> toProc(Pstream::nProcs());

    forAll(allBackgroundMeshBounds_(), proci)
    {
        // Test against the bounding box of the processor
        if
        (
            !isLocal(proci)
         && allBackgroundMeshBounds_()[proci].overlaps(centre, radiusSqr)
        )
        {
            toProc.append(proci);
        }
    }

    return toProc;
}


template<class Triangulation>
bool Foam::DistributedDelaunayMesh<Triangulation>::checkProcBoundaryCell
(
    const Cell_handle& cit,
    Map<labelList>& circumsphereOverlaps
) const
{
    const Foam::point& cc = cit->dual();

    const scalar crSqr = magSqr
    (
        cc - topoint(cit->vertex(0)->point())
    );

    labelList circumsphereOverlap = overlapProcessors
    (
        cc,
        sqr(1.01)*crSqr
    );

    cit->cellIndex() = this->getNewCellIndex();

    if (!circumsphereOverlap.empty())
    {
        circumsphereOverlaps.insert(cit->cellIndex(), circumsphereOverlap);

        return true;
    }

    return false;
}


template<class Triangulation>
void Foam::DistributedDelaunayMesh<Triangulation>::findProcessorBoundaryCells
(
    Map<labelList>& circumsphereOverlaps
) const
{
    // Start by assuming that all the cells have no index
    // If they do, they have already been visited so ignore them

    labelHashSet cellToCheck
    (
        Triangulation::number_of_finite_cells()
       /Pstream::nProcs()
    );

//    std::list<Cell_handle> infinite_cells;
//    Triangulation::incident_cells
//    (
//        Triangulation::infinite_vertex(),
//        std::back_inserter(infinite_cells)
//    );
//
//    for
//    (
//        typename std::list<Cell_handle>::iterator vcit
//            = infinite_cells.begin();
//        vcit != infinite_cells.end();
//        ++vcit
//    )
//    {
//        Cell_handle cit = *vcit;
//
//        // Index of infinite vertex in this cell.
//        label i = cit->index(Triangulation::infinite_vertex());
//
//        Cell_handle c = cit->neighbor(i);
//
//        if (c->unassigned())
//        {
//            c->cellIndex() = this->getNewCellIndex();
//
//            if (checkProcBoundaryCell(c, circumsphereOverlaps))
//            {
//                cellToCheck.insert(c->cellIndex());
//            }
//        }
//    }
//
//
//    for
//    (
//        Finite_cells_iterator cit = Triangulation::finite_cells_begin();
//        cit != Triangulation::finite_cells_end();
//        ++cit
//    )
//    {
//        if (cit->parallelDualVertex())
//        {
//            if (cit->unassigned())
//            {
//                if (checkProcBoundaryCell(cit, circumsphereOverlaps))
//                {
//                    cellToCheck.insert(cit->cellIndex());
//                }
//            }
//        }
//    }


    for
    (
        All_cells_iterator cit = Triangulation::all_cells_begin();
        cit != Triangulation::all_cells_end();
        ++cit
    )
    {
        if (Triangulation::is_infinite(cit))
        {
            // Index of infinite vertex in this cell.
            label i = cit->index(Triangulation::infinite_vertex());

            Cell_handle c = cit->neighbor(i);

            if (c->unassigned())
            {
                c->cellIndex() = this->getNewCellIndex();

                if (checkProcBoundaryCell(c, circumsphereOverlaps))
                {
                    cellToCheck.insert(c->cellIndex());
                }
            }
        }
        else if (cit->parallelDualVertex())
        {
            if (cit->unassigned())
            {
                if (checkProcBoundaryCell(cit, circumsphereOverlaps))
                {
                    cellToCheck.insert(cit->cellIndex());
                }
            }
        }
    }

    for
    (
        Finite_cells_iterator cit = Triangulation::finite_cells_begin();
        cit != Triangulation::finite_cells_end();
        ++cit
    )
    {
        if (cellToCheck.found(cit->cellIndex()))
        {
            // Get the neighbours and check them
            for (label adjCelli = 0; adjCelli < 4; ++adjCelli)
            {
                Cell_handle citNeighbor = cit->neighbor(adjCelli);

                // Ignore if has far point or previously visited
                if
                (
                    !citNeighbor->unassigned()
                 || !citNeighbor->internalOrBoundaryDualVertex()
                 || Triangulation::is_infinite(citNeighbor)
                )
                {
                    continue;
                }

                if
                (
                    checkProcBoundaryCell
                    (
                        citNeighbor,
                        circumsphereOverlaps
                    )
                )
                {
                    cellToCheck.insert(citNeighbor->cellIndex());
                }
            }

            cellToCheck.unset(cit->cellIndex());
        }
    }
}


template<class Triangulation>
void Foam::DistributedDelaunayMesh<Triangulation>::markVerticesToRefer
(
    const Map<labelList>& circumsphereOverlaps,
    PtrList<labelPairHashSet>& referralVertices,
    DynamicList<label>& targetProcessor,
    DynamicList<Vb>& parallelInfluenceVertices
)
{
    // Relying on the order of iteration of cells being the same as before
    for
    (
        Finite_cells_iterator cit = Triangulation::finite_cells_begin();
        cit != Triangulation::finite_cells_end();
        ++cit
    )
    {
        if (Triangulation::is_infinite(cit))
        {
             continue;
        }

        Map<labelList>::const_iterator iter =
            circumsphereOverlaps.find(cit->cellIndex());

        // Pre-tested circumsphere potential influence
        if (iter != circumsphereOverlaps.cend())
        {
            const labelList& citOverlaps = iter();

            forAll(citOverlaps, cOI)
            {
                label proci = citOverlaps[cOI];

                for (int i = 0; i < 4; i++)
                {
                    Vertex_handle v = cit->vertex(i);

                    if (v->farPoint())
                    {
                        continue;
                    }

                    label vProcIndex = v->procIndex();
                    label vIndex = v->index();

                    const labelPair procIndexPair(vProcIndex, vIndex);

                    // Using the hashSet to ensure that each vertex is only
                    // referred once to each processor.
                    // Do not refer a vertex to its own processor.
                    if (vProcIndex != proci)
                    {
                        if (referralVertices[proci].insert(procIndexPair))
                        {
                            targetProcessor.append(proci);

                            parallelInfluenceVertices.append
                            (
                                Vb
                                (
                                    v->point(),
                                    v->index(),
                                    v->type(),
                                    v->procIndex()
                                )
                            );

                            parallelInfluenceVertices.last().targetCellSize() =
                                v->targetCellSize();
                            parallelInfluenceVertices.last().alignment() =
                                v->alignment();
                        }
                    }
                }
            }
        }
    }
}


template<class Triangulation>
Foam::label Foam::DistributedDelaunayMesh<Triangulation>::referVertices
(
    const DynamicList<label>& targetProcessor,
    DynamicList<Vb>& parallelVertices,
    PtrList<labelPairHashSet>& referralVertices,
    labelPairHashSet& receivedVertices
)
{
    DynamicList<Vb> referredVertices(targetProcessor.size());

    const label preDistributionSize = parallelVertices.size();

    mapDistribute pointMap = buildMap(targetProcessor);

    // Make a copy of the original list.
    DynamicList<Vb> originalParallelVertices(parallelVertices);

    pointMap.distribute(parallelVertices);

    for (label proci = 0; proci < Pstream::nProcs(); proci++)
    {
        const labelList& constructMap = pointMap.constructMap()[proci];

        if (constructMap.size())
        {
            forAll(constructMap, i)
            {
                const Vb& v = parallelVertices[constructMap[i]];

                if
                (
                    v.procIndex() != Pstream::myProcNo()
                 && !receivedVertices.found(labelPair(v.procIndex(), v.index()))
                )
                {
                    referredVertices.append(v);

                    receivedVertices.insert
                    (
                        labelPair(v.procIndex(), v.index())
                    );
                }
            }
        }
    }

    label preInsertionSize = Triangulation::number_of_vertices();

    labelPairHashSet pointsNotInserted = rangeInsertReferredWithInfo
    (
        referredVertices.begin(),
        referredVertices.end(),
        true
    );

    if (!pointsNotInserted.empty())
    {
        for
        (
            typename labelPairHashSet::const_iterator iter
                = pointsNotInserted.begin();
            iter != pointsNotInserted.end();
            ++iter
        )
        {
            if (receivedVertices.found(iter.key()))
            {
                receivedVertices.erase(iter.key());
            }
        }
    }

    boolList pointInserted(parallelVertices.size(), true);

    forAll(parallelVertices, vI)
    {
        const labelPair procIndexI
        (
            parallelVertices[vI].procIndex(),
            parallelVertices[vI].index()
        );

        if (pointsNotInserted.found(procIndexI))
        {
            pointInserted[vI] = false;
        }
    }

    pointMap.reverseDistribute(preDistributionSize, pointInserted);

    forAll(originalParallelVertices, vI)
    {
        const label procIndex = targetProcessor[vI];

        if (!pointInserted[vI])
        {
            if (referralVertices[procIndex].size())
            {
                if
                (
                    !referralVertices[procIndex].unset
                    (
                        labelPair
                        (
                            originalParallelVertices[vI].procIndex(),
                            originalParallelVertices[vI].index()
                        )
                    )
                )
                {
                    Pout<< "*** not found "
                        << originalParallelVertices[vI].procIndex()
                        << " " << originalParallelVertices[vI].index() << endl;
                }
            }
        }
    }

    label postInsertionSize = Triangulation::number_of_vertices();

    reduce(preInsertionSize, sumOp<label>());
    reduce(postInsertionSize, sumOp<label>());

    label nTotalToInsert = referredVertices.size();

    reduce(nTotalToInsert, sumOp<label>());

    if (preInsertionSize + nTotalToInsert != postInsertionSize)
    {
        label nNotInserted =
            returnReduce(pointsNotInserted.size(), sumOp<label>());

        Info<< " Inserted = "
            << setw(name(label(Triangulation::number_of_finite_cells())).size())
            << nTotalToInsert - nNotInserted
            << " / " << nTotalToInsert << endl;

        nTotalToInsert -= nNotInserted;
    }
    else
    {
        Info<< " Inserted = " << nTotalToInsert << endl;
    }

    return nTotalToInsert;
}


template<class Triangulation>
void Foam::DistributedDelaunayMesh<Triangulation>::sync
(
    const boundBox& bb,
    PtrList<labelPairHashSet>& referralVertices,
    labelPairHashSet& receivedVertices,
    bool iterateReferral
)
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (allBackgroundMeshBounds_.empty())
    {
        distributeBoundBoxes(bb);
    }

    label nVerts = Triangulation::number_of_vertices();
    label nCells = Triangulation::number_of_finite_cells();

    DynamicList<Vb> parallelInfluenceVertices(0.1*nVerts);
    DynamicList<label> targetProcessor(0.1*nVerts);

    // Some of these values will not be used, i.e. for non-real cells
    DynamicList<Foam::point> circumcentre(0.1*nVerts);
    DynamicList<scalar> circumradiusSqr(0.1*nVerts);

    Map<labelList> circumsphereOverlaps(nCells);

    findProcessorBoundaryCells(circumsphereOverlaps);

    Info<< "    Influences = "
        << setw(name(nCells).size())
        << returnReduce(circumsphereOverlaps.size(), sumOp<label>()) << " / "
        << returnReduce(nCells, sumOp<label>());

    markVerticesToRefer
    (
        circumsphereOverlaps,
        referralVertices,
        targetProcessor,
        parallelInfluenceVertices
    );

    referVertices
    (
        targetProcessor,
        parallelInfluenceVertices,
        referralVertices,
        receivedVertices
    );

    if (iterateReferral)
    {
        label oldNReferred = 0;
        label nIterations = 1;

        Info<< incrIndent << indent
            << "Iteratively referring referred vertices..."
            << endl;
        do
        {
            Info<< indent << "Iteration " << nIterations++ << ":";

            circumsphereOverlaps.clear();
            targetProcessor.clear();
            parallelInfluenceVertices.clear();

            findProcessorBoundaryCells(circumsphereOverlaps);

            nCells = Triangulation::number_of_finite_cells();

            Info<< " Influences = "
                << setw(name(nCells).size())
                << returnReduce(circumsphereOverlaps.size(), sumOp<label>())
                << " / "
                << returnReduce(nCells, sumOp<label>());

            markVerticesToRefer
            (
                circumsphereOverlaps,
                referralVertices,
                targetProcessor,
                parallelInfluenceVertices
            );

            label nReferred = referVertices
            (
                targetProcessor,
                parallelInfluenceVertices,
                referralVertices,
                receivedVertices
            );

            if (nReferred == 0 || nReferred == oldNReferred)
            {
                break;
            }

            oldNReferred = nReferred;

        } while (true);

        Info<< decrIndent;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Triangulation>
Foam::scalar
Foam::DistributedDelaunayMesh<Triangulation>::calculateLoadUnbalance() const
{
    label nRealVertices = 0;

    for
    (
        Finite_vertices_iterator vit = Triangulation::finite_vertices_begin();
        vit != Triangulation::finite_vertices_end();
        ++vit
    )
    {
        // Only store real vertices that are not feature vertices
        if (vit->real() && !vit->featurePoint())
        {
            nRealVertices++;
        }
    }

    scalar globalNRealVertices = returnReduce
    (
        nRealVertices,
        sumOp<label>()
    );

    scalar unbalance = returnReduce
    (
        mag(1.0 - nRealVertices/(globalNRealVertices/Pstream::nProcs())),
        maxOp<scalar>()
    );

    Info<< "    Processor unbalance " << unbalance << endl;

    return unbalance;
}


template<class Triangulation>
bool Foam::DistributedDelaunayMesh<Triangulation>::distribute
(
    const boundBox& bb
)
{
    NotImplemented;

    if (!Pstream::parRun())
    {
        return false;
    }

    distributeBoundBoxes(bb);

    return true;
}


template<class Triangulation>
Foam::autoPtr<Foam::mapDistribute>
Foam::DistributedDelaunayMesh<Triangulation>::distribute
(
    const backgroundMeshDecomposition& decomposition,
    List<Foam::point>& points
)
{
    if (!Pstream::parRun())
    {
        return autoPtr<mapDistribute>();
    }

    distributeBoundBoxes(decomposition.procBounds());

    autoPtr<mapDistribute> mapDist = decomposition.distributePoints(points);

    return mapDist;
}


template<class Triangulation>
void Foam::DistributedDelaunayMesh<Triangulation>::sync(const boundBox& bb)
{
    if (!Pstream::parRun())
    {
        return;
    }

    if (allBackgroundMeshBounds_.empty())
    {
        distributeBoundBoxes(bb);
    }

    const label nApproxReferred =
        Triangulation::number_of_vertices()
       /Pstream::nProcs();

    PtrList<labelPairHashSet> referralVertices(Pstream::nProcs());
    forAll(referralVertices, proci)
    {
        if (!isLocal(proci))
        {
            referralVertices.set(proci, new labelPairHashSet(nApproxReferred));
        }
    }

    labelPairHashSet receivedVertices(nApproxReferred);

    sync
    (
        bb,
        referralVertices,
        receivedVertices,
        true
    );
}


template<class Triangulation>
template<class PointIterator>
typename Foam::DistributedDelaunayMesh<Triangulation>::labelPairHashSet
Foam::DistributedDelaunayMesh<Triangulation>::rangeInsertReferredWithInfo
(
    PointIterator begin,
    PointIterator end,
    bool printErrors
)
{
    const boundBox& bb = allBackgroundMeshBounds_()[Pstream::myProcNo()];

    typedef DynamicList
    <
        std::pair<scalar, label>
    > vectorPairPointIndex;

    vectorPairPointIndex pointsBbDistSqr;

    label count = 0;
    for (PointIterator it = begin; it != end; ++it)
    {
        const Foam::point samplePoint(topoint(it->point()));

        scalar distFromBbSqr = 0;

        if (!bb.contains(samplePoint))
        {
            const Foam::point nearestPoint = bb.nearest(samplePoint);

            distFromBbSqr = magSqr(nearestPoint - samplePoint);
        }

        pointsBbDistSqr.append
        (
            std::make_pair(distFromBbSqr, count++)
        );
    }

    std::random_shuffle(pointsBbDistSqr.begin(), pointsBbDistSqr.end());

    // Sort in ascending order by the distance of the point from the centre
    // of the processor bounding box
    sort(pointsBbDistSqr.begin(), pointsBbDistSqr.end());

    typename Triangulation::Vertex_handle hint;

    typename Triangulation::Locate_type lt;
    int li, lj;

    label nNotInserted = 0;

    labelPairHashSet uninserted
    (
        Triangulation::number_of_vertices()
       /Pstream::nProcs()
    );

    for
    (
        typename vectorPairPointIndex::const_iterator p =
            pointsBbDistSqr.begin();
        p != pointsBbDistSqr.end();
        ++p
    )
    {
        const size_t checkInsertion = Triangulation::number_of_vertices();

        const Vb& vert = *(begin + p->second);
        const Point& pointToInsert = vert.point();

        // Locate the point
        Cell_handle c = Triangulation::locate(pointToInsert, lt, li, lj, hint);

        bool inserted = false;

        if (lt == Triangulation::VERTEX)
        {
            if (printErrors)
            {
                Vertex_handle nearV =
                    Triangulation::nearest_vertex(pointToInsert);

                Pout<< "Failed insertion, point already exists" << nl
                    << "Failed insertion : " << vert.info()
                    << "         nearest : " << nearV->info();
            }
        }
        else if (lt == Triangulation::OUTSIDE_AFFINE_HULL)
        {
            WarningInFunction
                << "Point is outside affine hull! pt = " << pointToInsert
                << endl;
        }
        else if (lt == Triangulation::OUTSIDE_CONVEX_HULL)
        {
            // TODO: Can this be optimised?
            //
            // Only want to insert if a connection is formed between
            // pointToInsert and an internal or internal boundary point.
            hint = Triangulation::insert(pointToInsert, c);
            inserted = true;
        }
        else
        {
            // Get the cells that conflict with p in a vector V,
            // and a facet on the boundary of this hole in f.
            std::vector<Cell_handle> V;
            typename Triangulation::Facet f;

            Triangulation::find_conflicts
            (
                pointToInsert,
                c,
                CGAL::Oneset_iterator<typename Triangulation::Facet>(f),
                std::back_inserter(V)
            );

            for (size_t i = 0; i < V.size(); ++i)
            {
                Cell_handle conflictingCell = V[i];

                if
                (
                    Triangulation::dimension() < 3 // 2D triangulation
                 ||
                    (
                        !Triangulation::is_infinite(conflictingCell)
                     && (
                            conflictingCell->real()
                         || conflictingCell->hasFarPoint()
                        )
                    )
                )
                {
                    hint = Triangulation::insert_in_hole
                    (
                        pointToInsert,
                        V.begin(),
                        V.end(),
                        f.first,
                        f.second
                    );

                    inserted = true;

                    break;
                }
            }
        }

        if (inserted)
        {
            if (checkInsertion != Triangulation::number_of_vertices() - 1)
            {
                if (printErrors)
                {
                    Vertex_handle nearV =
                        Triangulation::nearest_vertex(pointToInsert);

                    Pout<< "Failed insertion : " << vert.info()
                        << "         nearest : " << nearV->info();
                }
            }
            else
            {
                hint->index() = vert.index();
                hint->type() = vert.type();
                hint->procIndex() = vert.procIndex();
                hint->targetCellSize() = vert.targetCellSize();
                hint->alignment() = vert.alignment();
            }
        }
        else
        {
            uninserted.insert(labelPair(vert.procIndex(), vert.index()));
            nNotInserted++;
        }
    }

    return uninserted;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
