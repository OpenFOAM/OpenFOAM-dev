/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "primitiveMesh.H"
#include "DynamicList.H"
#include "demandDrivenData.H"
#include "SortableList.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::primitiveMesh::getEdge
(
    List<DynamicList<label>>& pe,
    DynamicList<edge>& es,

    const label pointi,
    const label nextPointi
)
{
    // Find connection between pointi and nextPointi
    forAll(pe[pointi], ppI)
    {
        label eI = pe[pointi][ppI];

        const edge& e = es[eI];

        if (e.start() == nextPointi || e.end() == nextPointi)
        {
            return eI;
        }
    }

    // Make new edge.
    label edgeI = es.size();
    pe[pointi].append(edgeI);

    if (nextPointi != pointi)
    {
        // Very occasionally (e.g. blockMesh) a face can have duplicate
        // vertices. Make sure we register pointEdges only once.
        pe[nextPointi].append(edgeI);
    }

    if (pointi < nextPointi)
    {
        es.append(edge(pointi, nextPointi));
    }
    else
    {
        es.append(edge(nextPointi, pointi));
    }
    return edgeI;
}


void Foam::primitiveMesh::calcEdges(const bool doFaceEdges) const
{
    if (debug)
    {
        Pout<< "primitiveMesh::calcEdges(const bool) : "
            << "calculating edges, pointEdges and optionally faceEdges"
            << endl;
    }

    // It is an error to attempt to recalculate edges
    // if the pointer is already set
    if ((edgesPtr_ || pePtr_) || (doFaceEdges && fePtr_))
    {
        FatalErrorInFunction
            << "edges or pointEdges or faceEdges already calculated"
            << abort(FatalError);
    }
    else
    {
        // ALGORITHM:
        // Go through the faces list. Search pointEdges for existing edge.
        // If not found create edge and add to pointEdges.
        // In second loop sort edges in order of increasing neighbouring
        // point.
        // This algorithm replaces the one using pointFaces which used more
        // allocations but less memory and was on practical cases
        // quite a bit slower.

        const faceList& fcs = faces();

        // Size up lists
        // ~~~~~~~~~~~~~

        // Estimate pointEdges storage
        List<DynamicList<label>> pe(nPoints());
        forAll(pe, pointi)
        {
            pe[pointi].setCapacity(primitiveMesh::edgesPerPoint_);
        }

        // Estimate edges storage
        DynamicList<edge> es(pe.size()*primitiveMesh::edgesPerPoint_/2);

        // Estimate faceEdges storage
        if (doFaceEdges)
        {
            fePtr_ = new labelListList(fcs.size());
            labelListList& faceEdges = *fePtr_;
            forAll(fcs, facei)
            {
                faceEdges[facei].setSize(fcs[facei].size());
            }
        }


        // Find consecutive face points in edge list
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Edges on boundary faces
        label nExtEdges = 0;
        // Edges using no boundary point
        nInternal0Edges_ = 0;
        // Edges using 1 boundary point
        label nInt1Edges = 0;
        // Edges using two boundary points but not on boundary face:
        // edges.size()-nExtEdges-nInternal0Edges_-nInt1Edges

        // Ordering is different if points are ordered.
        if (nInternalPoints_ == -1)
        {
            // No ordering. No distinction between types.
            forAll(fcs, facei)
            {
                const face& f = fcs[facei];

                forAll(f, fp)
                {
                    label pointi = f[fp];
                    label nextPointi = f[f.fcIndex(fp)];

                    label edgeI = getEdge(pe, es, pointi, nextPointi);

                    if (doFaceEdges)
                    {
                        (*fePtr_)[facei][fp] = edgeI;
                    }
                }
            }
            // Assume all edges internal
            nExtEdges = 0;
            nInternal0Edges_ = es.size();
            nInt1Edges = 0;
        }
        else
        {
            // 1. Do external faces first. This creates external edges.
            for (label facei = nInternalFaces_; facei < fcs.size(); facei++)
            {
                const face& f = fcs[facei];

                forAll(f, fp)
                {
                    label pointi = f[fp];
                    label nextPointi = f[f.fcIndex(fp)];

                    label oldNEdges = es.size();
                    label edgeI = getEdge(pe, es, pointi, nextPointi);

                    if (es.size() > oldNEdges)
                    {
                        nExtEdges++;
                    }
                    if (doFaceEdges)
                    {
                        (*fePtr_)[facei][fp] = edgeI;
                    }
                }
            }

            // 2. Do internal faces. This creates internal edges.
            for (label facei = 0; facei < nInternalFaces_; facei++)
            {
                const face& f = fcs[facei];

                forAll(f, fp)
                {
                    label pointi = f[fp];
                    label nextPointi = f[f.fcIndex(fp)];

                    label oldNEdges = es.size();
                    label edgeI = getEdge(pe, es, pointi, nextPointi);

                    if (es.size() > oldNEdges)
                    {
                        if (pointi < nInternalPoints_)
                        {
                            if (nextPointi < nInternalPoints_)
                            {
                                nInternal0Edges_++;
                            }
                            else
                            {
                                nInt1Edges++;
                            }
                        }
                        else
                        {
                            if (nextPointi < nInternalPoints_)
                            {
                                nInt1Edges++;
                            }
                            else
                            {
                                // Internal edge with two points on boundary
                            }
                        }
                    }
                    if (doFaceEdges)
                    {
                        (*fePtr_)[facei][fp] = edgeI;
                    }
                }
            }
        }


        // For unsorted meshes the edges will be in order of occurrence inside
        // the faces. For sorted meshes the first nExtEdges will be the external
        // edges.

        if (nInternalPoints_ != -1)
        {
            nInternalEdges_ = es.size()-nExtEdges;
            nInternal1Edges_ = nInternal0Edges_+nInt1Edges;

            //Pout<< "Edge overview:" << nl
            //    << "    total number of edges           : " << es.size()
            //    << nl
            //    << "    boundary edges                  : " << nExtEdges
            //    << nl
            //    << "    edges using no boundary point   : "
            //    << nInternal0Edges_
            //    << nl
            //    << "    edges using one boundary points : " << nInt1Edges
            //   << nl
            //    << "    edges using two boundary points : "
            //    << es.size()-nExtEdges-nInternal0Edges_-nInt1Edges << endl;
        }


        // Like faces sort edges in order of increasing neighbouring point.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Automatically if points are sorted into internal and external points
        // the edges will be sorted into
        // - internal point to internal point
        // - internal point to external point
        // - external point to external point
        // The problem is that the last one mixes external edges with internal
        // edges with two points on the boundary.


        // Map to sort into new upper-triangular order
        labelList oldToNew(es.size());
        if (debug)
        {
            oldToNew = -1;
        }

        // start of edges with 0 boundary points
        label internal0EdgeI = 0;

        // start of edges with 1 boundary points
        label internal1EdgeI = nInternal0Edges_;

        // start of edges with 2 boundary points
        label internal2EdgeI = nInternal1Edges_;

        // start of external edges
        label externalEdgeI = nInternalEdges_;


        // To sort neighbouring points in increasing order. Defined outside
        // for optimisation reasons: if all connectivity size same will need
        // no reallocations
        SortableList<label> nbrPoints(primitiveMesh::edgesPerPoint_);

        forAll(pe, pointi)
        {
            const DynamicList<label>& pEdges = pe[pointi];

            nbrPoints.setSize(pEdges.size());

            forAll(pEdges, i)
            {
                const edge& e = es[pEdges[i]];

                label nbrPointi = e.otherVertex(pointi);

                if (nbrPointi < pointi)
                {
                    nbrPoints[i] = -1;
                }
                else
                {
                    nbrPoints[i] = nbrPointi;
                }
            }
            nbrPoints.sort();


            if (nInternalPoints_ == -1)
            {
                // Sort all upper-triangular
                forAll(nbrPoints, i)
                {
                    if (nbrPoints[i] != -1)
                    {
                        label edgeI = pEdges[nbrPoints.indices()[i]];

                        oldToNew[edgeI] = internal0EdgeI++;
                    }
                }
            }
            else
            {
                if (pointi < nInternalPoints_)
                {
                    forAll(nbrPoints, i)
                    {
                        label nbrPointi = nbrPoints[i];

                        label edgeI = pEdges[nbrPoints.indices()[i]];

                        if (nbrPointi != -1)
                        {
                            if (edgeI < nExtEdges)
                            {
                                // External edge
                                oldToNew[edgeI] = externalEdgeI++;
                            }
                            else if (nbrPointi < nInternalPoints_)
                            {
                                // Both points inside
                                oldToNew[edgeI] = internal0EdgeI++;
                            }
                            else
                            {
                                // One points inside, one outside
                                oldToNew[edgeI] = internal1EdgeI++;
                            }
                        }
                    }
                }
                else
                {
                    forAll(nbrPoints, i)
                    {
                        label nbrPointi = nbrPoints[i];

                        label edgeI = pEdges[nbrPoints.indices()[i]];

                        if (nbrPointi != -1)
                        {
                            if (edgeI < nExtEdges)
                            {
                                // External edge
                                oldToNew[edgeI] = externalEdgeI++;
                            }
                            else if (nbrPointi < nInternalPoints_)
                            {
                                // Not possible!
                                FatalErrorInFunction
                                    << abort(FatalError);
                            }
                            else
                            {
                                // Both points outside
                                oldToNew[edgeI] = internal2EdgeI++;
                            }
                        }
                    }
                }
            }
        }


        if (debug)
        {
            label edgeI = findIndex(oldToNew, -1);

            if (edgeI != -1)
            {
                const edge& e = es[edgeI];

                FatalErrorInFunction
                    << "Did not sort edge " << edgeI << " points:" << e
                    << " coords:" << points()[e[0]] << points()[e[1]]
                    << endl
                    << "Current buckets:" << endl
                    << "    internal0EdgeI:" << internal0EdgeI << endl
                    << "    internal1EdgeI:" << internal1EdgeI << endl
                    << "    internal2EdgeI:" << internal2EdgeI << endl
                    << "    externalEdgeI:" << externalEdgeI << endl
                    << abort(FatalError);
            }
        }



        // Renumber and transfer.

        // Edges
        edgesPtr_ = new edgeList(es.size());
        edgeList& edges = *edgesPtr_;
        forAll(es, edgeI)
        {
            edges[oldToNew[edgeI]] = es[edgeI];
        }

        // pointEdges
        pePtr_ = new labelListList(nPoints());
        labelListList& pointEdges = *pePtr_;
        forAll(pe, pointi)
        {
            DynamicList<label>& pEdges = pe[pointi];
            pEdges.shrink();
            inplaceRenumber(oldToNew, pEdges);
            pointEdges[pointi].transfer(pEdges);
            Foam::sort(pointEdges[pointi]);
        }

        // faceEdges
        if (doFaceEdges)
        {
            labelListList& faceEdges = *fePtr_;
            forAll(faceEdges, facei)
            {
                inplaceRenumber(oldToNew, faceEdges[facei]);
            }
        }
    }
}


Foam::label Foam::primitiveMesh::findFirstCommonElementFromSortedLists
(
    const labelList& list1,
    const labelList& list2
)
{
    label result = -1;

    labelList::const_iterator iter1 = list1.begin();
    labelList::const_iterator iter2 = list2.begin();

    while (iter1 != list1.end() && iter2 != list2.end())
    {
        if (*iter1 < *iter2)
        {
            ++iter1;
        }
        else if (*iter1 > *iter2)
        {
            ++iter2;
        }
        else
        {
            result = *iter1;
            break;
        }
    }
    if (result == -1)
    {
        FatalErrorInFunction
            << "No common elements in lists " << list1 << " and " << list2
            << abort(FatalError);
    }
    return result;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::edgeList& Foam::primitiveMesh::edges() const
{
    if (!edgesPtr_)
    {
        //calcEdges(true);
        calcEdges(false);
    }

    return *edgesPtr_;
}

const Foam::labelListList& Foam::primitiveMesh::pointEdges() const
{
    if (!pePtr_)
    {
        //calcEdges(true);
        calcEdges(false);
    }

    return *pePtr_;
}


const Foam::labelListList& Foam::primitiveMesh::faceEdges() const
{
    if (!fePtr_)
    {
        if (debug)
        {
            Pout<< "primitiveMesh::faceEdges() : "
                << "calculating faceEdges" << endl;
        }

        //calcEdges(true);
        const faceList& fcs = faces();
        const labelListList& pe = pointEdges();
        const edgeList& es = edges();

        fePtr_ = new labelListList(fcs.size());
        labelListList& faceEdges = *fePtr_;

        forAll(fcs, facei)
        {
            const face& f = fcs[facei];

            labelList& fEdges = faceEdges[facei];
            fEdges.setSize(f.size());

            forAll(f, fp)
            {
                label pointi = f[fp];
                label nextPointi = f[f.fcIndex(fp)];

                // Find edge between pointi, nextPontI
                const labelList& pEdges = pe[pointi];

                forAll(pEdges, i)
                {
                    label edgeI = pEdges[i];

                    if (es[edgeI].otherVertex(pointi) == nextPointi)
                    {
                        fEdges[fp] = edgeI;
                        break;
                    }
                }
            }
        }
    }

    return *fePtr_;
}


void Foam::primitiveMesh::clearOutEdges()
{
    deleteDemandDrivenData(edgesPtr_);
    deleteDemandDrivenData(pePtr_);
    deleteDemandDrivenData(fePtr_);
    labels_.clear();
    labelSet_.clear();
}


const Foam::labelList& Foam::primitiveMesh::faceEdges
(
    const label facei,
    DynamicList<label>& storage
) const
{
    if (hasFaceEdges())
    {
        return faceEdges()[facei];
    }
    else
    {
        const labelListList& pointEs = pointEdges();
        const face& f = faces()[facei];

        storage.clear();
        if (f.size() > storage.capacity())
        {
            storage.setCapacity(f.size());
        }

        forAll(f, fp)
        {
            storage.append
            (
                findFirstCommonElementFromSortedLists
                (
                    pointEs[f[fp]],
                    pointEs[f.nextLabel(fp)]
                )
            );
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::faceEdges(const label facei) const
{
    return faceEdges(facei, labels_);
}


const Foam::labelList& Foam::primitiveMesh::cellEdges
(
    const label celli,
    DynamicList<label>& storage
) const
{
    if (hasCellEdges())
    {
        return cellEdges()[celli];
    }
    else
    {
        const labelList& cFaces = cells()[celli];

        labelSet_.clear();

        forAll(cFaces, i)
        {
            const labelList& fe = faceEdges(cFaces[i]);

            forAll(fe, feI)
            {
                labelSet_.insert(fe[feI]);
            }
        }

        storage.clear();

        if (labelSet_.size() > storage.capacity())
        {
            storage.setCapacity(labelSet_.size());
        }

        forAllConstIter(labelHashSet, labelSet_, iter)
        {
            storage.append(iter.key());
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::cellEdges(const label celli) const
{
    return cellEdges(celli, labels_);
}


// ************************************************************************* //
