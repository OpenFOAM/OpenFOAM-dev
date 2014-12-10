/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "geomCellLooper.H"
#include "polyMesh.H"
#include "DynamicList.H"
#include "plane.H"
#include "meshTools.H"
#include "SortableList.H"
#include "triSurfaceTools.H"
#include "HashSet.H"
#include "ListOps.H"
#include "transform.H"

#include "addToRunTimeSelectionTable.H"


// Extension factor of edges to make sure we catch intersections through
// edge endpoints
const Foam::scalar Foam::geomCellLooper::pointEqualTol_ = 1e-3;


// Snap cuts through edges onto edge endpoints. Fraction of edge length.
Foam::scalar Foam::geomCellLooper::snapTol_ = 0.1;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(geomCellLooper, 0);
addToRunTimeSelectionTable(cellLooper, geomCellLooper, word);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::geomCellLooper::minEdgeLen(const label vertI) const
{
    scalar minLen = GREAT;

    const labelList& pEdges = mesh().pointEdges()[vertI];

    forAll(pEdges, pEdgeI)
    {
        const edge& e = mesh().edges()[pEdges[pEdgeI]];

        minLen = min(minLen, e.mag(mesh().points()));
    }
    return minLen;
}


// Cut edge with plane. Return true and set weight to fraction between
// edge-start and edge-end
bool Foam::geomCellLooper::cutEdge
(
    const plane& cutPlane,
    const label edgeI,
    scalar& weight
) const
{
    const pointField& pts = mesh().points();

    const edge& e = mesh().edges()[edgeI];

    scalar s = cutPlane.normalIntersect(pts[e.start()], e.vec(pts));

    if ((s > -pointEqualTol_) && (s < 1 + pointEqualTol_))
    {
        weight = s;

        return true;
    }
    else
    {
        // Make sure we don't use this value
        weight = -GREAT;

        return false;
    }
}


// If edge close enough to endpoint snap to endpoint.
Foam::label Foam::geomCellLooper::snapToVert
(
    const scalar tol,
    const label edgeI,
    const scalar weight
) const
{
    const edge& e = mesh().edges()[edgeI];

    if (weight < tol)
    {
        return e.start();
    }
    else if (weight > (1-tol))
    {
        return e.end();
    }
    else
    {
        return -1;
    }
}


void Foam::geomCellLooper::getBase(const vector& n, vector& e0, vector& e1)
 const
{
    // Guess for vector normal to n.
    vector base(1,0,0);

    scalar nComp = n & base;

    if (mag(nComp) > 0.8)
    {
        // Was bad guess. Try with different vector.

        base.x() = 0;
        base.y() = 1;

        nComp = n & base;

        if (mag(nComp) > 0.8)
        {
            base.y() = 0;
            base.z() = 1;

            nComp = n & base;
        }
    }


    // Use component normal to n as base vector.
    e0 = base - nComp*n;

    e0 /= mag(e0) + VSMALL;

    e1 = n ^ e0;

    //Pout<< "Coord system:" << endl
    //    << "    n  : " << n << ' ' << mag(n) << endl
    //    << "    e0 : " << e0 << ' ' << mag(e0) << endl
    //    << "    e1 : " << e1 << ' ' << mag(e1) << endl
    //    << endl;
}


// Return true if the cut edge at loop[index] is preceded by cuts through
// the edge end points.
bool Foam::geomCellLooper::edgeEndsCut
(
    const labelList& loop,
    const label index
) const
{
    label edgeI = getEdge(loop[index]);

    const edge& e = mesh().edges()[edgeI];

    const label prevCut = loop[loop.rcIndex(index)];
    const label nextCut = loop[loop.fcIndex(index)];

    if (!isEdge(prevCut) && !isEdge(nextCut))
    {
        // Cuts before and after are both vertices. Check if both
        // the edge endpoints
        label v0 = getVertex(prevCut);
        label v1 = getVertex(nextCut);

        if
        (
            (v0 == e[0] && v1 == e[1])
         || (v0 == e[1] && v1 == e[0])
        )
        {
            return true;
        }
    }
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::geomCellLooper::geomCellLooper(const polyMesh& mesh)
:
    cellLooper(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::geomCellLooper::~geomCellLooper()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::geomCellLooper::cut
(
    const vector& refDir,
    const label cellI,
    const boolList& vertIsCut,
    const boolList& edgeIsCut,
    const scalarField& edgeWeight,

    labelList& loop,
    scalarField& loopWeights
) const
{
    // Cut through cell centre normal to refDir.
    return cut
    (
        plane(mesh().cellCentres()[cellI], refDir),
        cellI,
        vertIsCut,
        edgeIsCut,
        edgeWeight,
        loop,
        loopWeights
    );
}


bool Foam::geomCellLooper::cut
(
    const plane& cutPlane,
    const label cellI,
    const boolList&,
    const boolList&,
    const scalarField&,

    labelList& loop,
    scalarField& loopWeights
) const
{
    const pointField& points = mesh().points();
    const edgeList& edges = mesh().edges();

    // Find all cuts through edges.
    // Special cases:
    // - edge cut close to endpoint. Snap to endpoint.
    // - edge endpoint close to plane (but edge not cut). Snap to endpoint.
    // - both endpoints close to plane. Edge parallel to plane. No need to
    //   cut to edge.
    // Note: any snap-to-point check must be based on all edges using a point
    // since otherwise cut through close to point but on neighbouring edge
    // might not be snapped.

    // Size overly big.
    label nEstCuts = 2*mesh().cells()[cellI].size();

    DynamicList<label> localLoop(nEstCuts);
    DynamicList<scalar> localLoopWeights(nEstCuts);

    // Points checked. Used to make sure we don't cut edge and edge endpoints
    // at the same time.
    labelHashSet checkedPoints(nEstCuts);

    const labelList& cellEdges = mesh().cellEdges()[cellI];

    forAll(cellEdges, i)
    {
        label edgeI = cellEdges[i];

        const edge& e = edges[edgeI];

        bool useStart = false;

        bool useEnd = false;

        //
        // Check distance of endpoints to cutPlane
        //

        if (!checkedPoints.found(e.start()))
        {
            checkedPoints.insert(e.start());

            scalar typStartLen = pointEqualTol_ * minEdgeLen(e.start());

            // Check distance of startPt to plane.
            if (cutPlane.distance(points[e.start()]) < typStartLen)
            {
                // Use point.
                localLoop.append(vertToEVert(e.start()));
                localLoopWeights.append(-GREAT);

                useStart = true;
            }
        }
        if (!checkedPoints.found(e.end()))
        {
            checkedPoints.insert(e.end());

            scalar typEndLen = pointEqualTol_ * minEdgeLen(e.end());

            // Check distance of endPt to plane.
            if (cutPlane.distance(points[e.end()]) < typEndLen)
            {
                // Use point.
                localLoop.append(vertToEVert(e.end()));
                localLoopWeights.append(-GREAT);

                useEnd = true;
            }
        }

        //
        // Check cut of edge
        //

        if (!useEnd && !useStart)
        {
            // Edge end points not close to plane so edge not close to
            // plane. Cut edge.
            scalar cutWeight;

            if (cutEdge(cutPlane, edgeI, cutWeight))
            {
                // Snap edge cut to vertex.
                label cutVertI = snapToVert(snapTol_, edgeI, cutWeight);

                if (cutVertI == -1)
                {
                    // Proper cut of edge
                    localLoop.append(edgeToEVert(edgeI));
                    localLoopWeights.append(cutWeight);
                }
                else
                {
                    // Cut through edge end point. Might be duplicate
                    // since connected edge might also be snapped to same
                    // endpoint so only insert if unique.
                    label cut = vertToEVert(cutVertI);

                    if (findIndex(localLoop, cut) == -1)
                    {
                        localLoop.append(vertToEVert(cutVertI));
                        localLoopWeights.append(-GREAT);
                    }
                }
            }
        }
    }

    if (localLoop.size() <= 2)
    {
        return false;
    }

    localLoop.shrink();
    localLoopWeights.shrink();


    // Get points on loop and centre of loop
    pointField loopPoints(localLoop.size());
    point ctr(vector::zero);

    forAll(localLoop, i)
    {
        loopPoints[i] = coord(localLoop[i], localLoopWeights[i]);
        ctr += loopPoints[i];
    }
    ctr /= localLoop.size();


    // Get base vectors of coordinate system normal to refDir
    vector e0, e1;
    getBase(cutPlane.normal(), e0, e1);


    // Get sorted angles from point on loop to centre of loop.
    SortableList<scalar> sortedAngles(localLoop.size());

    forAll(sortedAngles, i)
    {
        vector toCtr(loopPoints[i] - ctr);
        toCtr /= mag(toCtr);

        sortedAngles[i] = pseudoAngle(e0, e1, toCtr);
    }
    sortedAngles.sort();

    loop.setSize(localLoop.size());
    loopWeights.setSize(loop.size());


    // Fill loop and loopweights according to sorted angles

    const labelList& indices = sortedAngles.indices();

    forAll(indices, i)
    {
        loop[i] = localLoop[indices[i]];
        loopWeights[i] = localLoopWeights[indices[i]];
    }


    // Check for cut edges along already cut vertices.
    bool filterLoop = false;

    forAll(loop, i)
    {
        label cut = loop[i];

        if (isEdge(cut) && edgeEndsCut(loop, i))
        {
            filterLoop = true;
        }
    }

    if (filterLoop)
    {
        // Found edges in loop where both end points are also in loop. This
        // is superfluous information so we can remove it.

        labelList filteredLoop(loop.size());
        scalarField filteredLoopWeights(loopWeights.size());
        label filterI = 0;

        forAll(loop, i)
        {
            label cut = loop[i];

            if (isEdge(cut) && edgeEndsCut(loop, i))
            {
                // Superfluous cut. Edge end points cut and edge itself as well.
            }
            else
            {
                filteredLoop[filterI] = loop[i];
                filteredLoopWeights[filterI] = loopWeights[i];
                filterI++;
            }
        }
        filteredLoop.setSize(filterI);
        filteredLoopWeights.setSize(filterI);

        loop.transfer(filteredLoop);
        loopWeights.transfer(filteredLoopWeights);
    }

    if (debug&2)
    {
        Pout<< "cell:" << cellI << endl;
        forAll(loop, i)
        {
            Pout<< "At angle:" << sortedAngles[i] << endl
                << "    cut:";

            writeCut(Pout, loop[i], loopWeights[i]);

            Pout<< "  coord:" << coord(loop[i], loopWeights[i]);

            Pout<< endl;
        }
    }

    return true;
}


// ************************************************************************* //
