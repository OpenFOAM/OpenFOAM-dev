/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
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

#include "edgeSurface.H"
#include "triSurface.H"
#include "surfaceIntersection.H"
#include "meshTools.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(edgeSurface, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Write whole pointField and edges to stream
void Foam::edgeSurface::writeOBJ
(
    const pointField& points,
    const edgeList& edges,
    Ostream& os
)
{
    forAll(points, pointi)
    {
        const point& pt = points[pointi];

        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }
    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        os << "l " << e.start()+1 << ' ' << e.end()+1 << endl;
    }
}


// Write whole pointField and selected edges to stream
void Foam::edgeSurface::writeOBJ
(
    const pointField& points,
    const edgeList& edges,
    const labelList& edgeLabels,
    Ostream& os
)
{
    forAll(points, pointi)
    {
        const point& pt = points[pointi];

        os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
    }
    forAll(edgeLabels, i)
    {
        const edge& e = edges[edgeLabels[i]];

        os << "l " << e.start()+1 << ' ' << e.end()+1 << endl;
    }
}


// Pointedges in edgeSurface indices only.
void Foam::edgeSurface::calcPointEdges()
{
    pointEdges_.setSize(points_.size());

    labelList pointNEdges(points_.size(), 0);

    forAll(edges_, edgeI)
    {
        const edge& e = edges_[edgeI];

        pointNEdges[e[0]]++;
        pointNEdges[e[1]]++;
    }

    forAll(pointEdges_, pointi)
    {
        pointEdges_[pointi].setSize(pointNEdges[pointi]);
    }

    pointNEdges = 0;

    forAll(edges_, edgeI)
    {
        const edge& e = edges_[edgeI];

        labelList& pEdges0 = pointEdges_[e[0]];
        pEdges0[pointNEdges[e[0]]++] = edgeI;

        labelList& pEdges1 = pointEdges_[e[1]];
        pEdges1[pointNEdges[e[1]]++] = edgeI;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from surface and intersection description
Foam::edgeSurface::edgeSurface
(
    const triSurface& surf,
    const bool isFirstSurface,
    const surfaceIntersection& inter
)
:
    points_(surf.nPoints() + inter.cutPoints().size()),
    nSurfacePoints_(surf.nPoints()),
    edges_(),
    nSurfaceEdges_(surf.nEdges()),
    parentEdges_(0),
    faceEdges_(surf.size()),
    pointEdges_(points_.size())
{
    // Copy points (surface ones first)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label pointi = 0;

    const pointField& surfPoints = surf.localPoints();

    forAll(surfPoints, i)
    {
        points_[pointi++] = surfPoints[i];
    }

    const pointField& cutPoints = inter.cutPoints();

    forAll(cutPoints, i)
    {
        points_[pointi++] = cutPoints[i];
    }


    // Copy edges (surface ones first)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<edge> allEdges(surf.nEdges() + inter.cutEdges().size());
    DynamicList<label> allParentEdges(surf.nEdges());
    List<DynamicList<label>> allFaceEdges(surf.size());


    // Copy surface edges (can be split!)

    const edgeList& surfEdges = surf.edges();

    forAll(surfEdges, edgeI)
    {
        const edge& e = surfEdges[edgeI];

        // Get additional vertices for this edge.
        const labelList& extraVerts = inter.edgeCuts(isFirstSurface)[edgeI];

        // Store current top of allEdges.
        label freeNewEdgeI = allEdges.size();

        if (extraVerts.empty())
        {
            // No cuts across this edge. Note that vertices do not need to be
            // renumbered.
            allEdges.append(e);
        }
        else
        {
            // Edge is cut. From e.start() to extraVerts[0],
            // from extraVerts[i] to i+1 and finally to e.end().
            allEdges.append
            (
                edge
                (
                    e.start(),
                    extraVerts[0] + nSurfacePoints_
                )
            );

            for (label extraI = 1; extraI < extraVerts.size(); extraI++)
            {
                allEdges.append
                (
                    edge
                    (
                        extraVerts[extraI-1] + nSurfacePoints_,
                        extraVerts[extraI] + nSurfacePoints_
                    )
                );
            }
            allEdges.append
            (
                edge
                (
                    extraVerts.last() + nSurfacePoints_,
                    e.end()
                )
            );
        }

        // Update allFaceEdges, parentEdges_ for the newly added edges.

        // Add each edge label to all face neighbours of edgeI
        const labelList& myFaces = surf.edgeFaces()[edgeI];

        for (label eI = freeNewEdgeI; eI < allEdges.size(); eI++)
        {
            allParentEdges.append(edgeI);

            forAll(myFaces, myFacei)
            {
                allFaceEdges[myFaces[myFacei]].append(eI);
            }
        }
    }

    // Done all (possibly split) surface edges by now.
    nSurfaceEdges_ = allEdges.size();


    // Copy intersection edges
    // (note no parentEdges)
    const edgeList& cutEdges = inter.cutEdges();

    forAll(cutEdges, i)
    {
        const edge& e = cutEdges[i];

        allEdges.append(edge(e[0] + nSurfacePoints_, e[1] + nSurfacePoints_));
    }




    // Add intersection edges to faceEdges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAllConstIter(labelPairLookup, inter.facePairToEdge(), iter)
    {
        // Edge label in intersection
        const label edgeI = iter();

        // Get the face from the correct surface
        const FixedList<label, 2>& twoFaces = iter.key();

        label facei;

        if (isFirstSurface)
        {
            facei = twoFaces[0];
        }
        else
        {
            facei = twoFaces[1];
        }

        // Store on face-edge addressing. (note: offset edge)
        allFaceEdges[facei].append(edgeI + nSurfaceEdges_);
    }

    // Transfer.
    edges_.transfer(allEdges);
    parentEdges_.transfer(allParentEdges);

    forAll(allFaceEdges, facei)
    {
        faceEdges_[facei].transfer(allFaceEdges[facei]);
    }


    // Additional addressing
    // ~~~~~~~~~~~~~~~~~~~~~

    calcPointEdges();


    if (debug & 4)
    {
        Pout<< "edgeSurface : Dumping faceEdges to files" << endl;

        forAll(faceEdges_, facei)
        {
            const labelList& fEdges = faceEdges_[facei];

            if (fEdges.size() != 3)
            {
                fileName faceFName("face_" + name(facei) + ".obj");
                Pout<< "edgeSurface : Dumping faceEdges for face " << facei
                    << " to " << faceFName << endl;

                OFstream fStream(faceFName);
                writeOBJ(points_, edges_, fEdges, fStream);
            }
        }

        Pout<< "edgeSurface : Dumping edges to edges.obj" << endl;
        OFstream eStream("edges.obj");
        writeOBJ(points_, edges_, eStream);

        Pout<< "edgeSurface : Dumping intersectionEdges to"
            << " intersectionEdges.obj" << endl;
        OFstream intEdgesStream("intersectionEdges.obj");

        labelList edgeLabels(edges_.size() - nSurfaceEdges_);

        label i = 0;
        for (label edgeI = nSurfaceEdges_; edgeI < edges_.size(); edgeI++)
        {
            edgeLabels[i++] = edgeI;
        }

        writeOBJ(points_, edges_, edgeLabels, intEdgesStream);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::edgeSurface::addIntersectionEdges
(
    const label facei,
    const edgeList& additionalEdges
)
{
    if (debug & 2)
    {
        Pout<< "Old face consisted of edges:" << endl;

        const labelList& fEdges = faceEdges_[facei];
        forAll(fEdges, i)
        {
            const edge& e = edges_[fEdges[i]];

            Pout<< "    " << fEdges[i] << ' ' << e
                << points_[e.start()] << ' ' << points_[e.end()] << endl;
        }
    }

    // Make space for additional intersection edges (copies old ones)
    const label oldNEdges = edges_.size();

    edges_.setSize(oldNEdges + additionalEdges.size());

    // Append new intersection edges
    label newEdgeI = oldNEdges;

    forAll(additionalEdges, i)
    {
        edges_[newEdgeI] = additionalEdges[i];  // Vertices already in eSurf
                                                // indices.
        newEdgeI++;
    }

    // Append to faceEdges.
    labelList& fEdges = faceEdges_[facei];

    label nFEdges = fEdges.size();

    fEdges.setSize(nFEdges + additionalEdges.size());

    forAll(additionalEdges, i)
    {
        fEdges[nFEdges++] = oldNEdges + i;
    }


    // Update pointEdge addressing
    calcPointEdges();


    if (debug & 2)
    {
        const labelList& fEdges = faceEdges_[facei];

        Pout<< "New face consists of edges:" << endl;
        forAll(fEdges, i)
        {
            const edge& e = edges_[fEdges[i]];

            Pout<< "    " << fEdges[i] << ' ' << e
                << points_[e.start()] << ' ' << points_[e.end()] << endl;
        }
    }
}


// ************************************************************************* //
