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

Description

\*---------------------------------------------------------------------------*/

#include "router.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::router::count(const label weight) const
{
    label cnt = 0;

    forAll(weights_, nodeI)
    {
        cnt += weights_[nodeI];
    }

    return cnt;
}


// Given connections between nodes set minimum distance from nodeI
void Foam::router::setWeights
(
    const label weight,
    const label nodeI
)
{
    // Set weight at current node
    weights_[nodeI] = weight;

    const labelList& myNeighbours = connections_[nodeI];

    forAll(myNeighbours, neighbourI)
    {
        if (weights_[myNeighbours[neighbourI]] > weight + 1)
        {
            // Distribute weight+1 to neighbours
            setWeights
            (
                weight+1,
                myNeighbours[neighbourI]
            );
        }
    }
}


// Mark shortest path from endNode to startNode by setting the weights
// to 0.
void Foam::router::fixWeights
(
    const label startNodeI,
    const label endNodeI,

    const label nodeI,
    const label prevNodeI
)
{
    // Mark this node
    weights_[nodeI] = 0;

    label minNodeI = -1;
    label minDist = labelMax;
    label nMinNodes = 0;

    const labelList& myNeighbours = connections_[nodeI];

    forAll(myNeighbours, neighbourI)
    {
        label nbrNodeI = myNeighbours[neighbourI];

        if (nbrNodeI != prevNodeI)
        {
            if (weights_[nbrNodeI] == 0)
            {
                // Reached end
                minDist = 0;
                break;
            }
            else if (weights_[nbrNodeI] > 0)
            {
                if (weights_[nbrNodeI] < minDist)
                {
                    minDist = weights_[nbrNodeI];
                    minNodeI = nbrNodeI;
                    nMinNodes = 1;
                }
                else if (weights_[nbrNodeI] == minDist)
                {
                    nMinNodes++;
                }
            }
        }
    }

    if (minDist == 0)
    {
        // Reached starting point.
        return;
    }

    if (minNodeI == -1)
    {
        WarningInFunction
            << "Cannot route from node " << nodeI
            << " since all neighbours of node "
            << "already allocated:" << endl;

        forAll(myNeighbours, neighbourI)
        {
            label nbrNodeI = myNeighbours[neighbourI];

            WarningInFunction
                << "  weight:" << weights_[nbrNodeI] << endl;
        }
        return;
    }


    if (nMinNodes > 1)
    {
        // Multiple paths, all with same weight. Use some heuristic
        // to choose one. Here: smallest angle to vector end-start
        vector n(coords_[endNodeI] - coords_[startNodeI]);

        scalar maxCosAngle = -great;

        forAll(myNeighbours, neighbourI)
        {
            label nbrNodeI = myNeighbours[neighbourI];

            if (weights_[nbrNodeI] == minDist)
            {
                vector n2(coords_[nbrNodeI] - coords_[startNodeI]);

                scalar magN2 = mag(n2);

                if (magN2 > small)
                {
                    n2 /= mag(n2);

                    scalar cosAngle = n & n2;

                    if (cosAngle > maxCosAngle)
                    {
                        maxCosAngle = cosAngle;
                        minNodeI = nbrNodeI;
                    }
                }
            }
        }
    }


    // Recursively go mark the path at minNode
    fixWeights
    (
        startNodeI,
        endNodeI,
        minNodeI,
        nodeI
    );
}


Foam::label Foam::router::getValue(const label pathValue) const
{
    forAll(weights_, nodeI)
    {
        if (weights_[nodeI] == pathValue)
        {
            return nodeI;
        }
    }
    return -1;
}


// Find node which has no neighbours with pathValue
Foam::label Foam::router::findEndNode
(
    const label startNodeI,
    const label prevNodeI,
    const label pathValue
) const
{
    const labelList& myNeighbours = connections_[startNodeI];

    forAll(myNeighbours, neighbourI)
    {
        label nodeI = myNeighbours[neighbourI];

        if (nodeI != prevNodeI)
        {
            if (weights_[nodeI] == pathValue)
            {
                return findEndNode(nodeI, startNodeI, pathValue);
            }
        }
    }

    // No neighbours with pathValue found. Return this node
    return startNodeI;
}


// Append all pathValue weights to route.
void Foam::router::storeRoute
(
    const label startNodeI,
    const label prevNodeI,
    const label pathValue,
    DynamicList<label>& route
) const
{
    const labelList& myNeighbours = connections_[startNodeI];

    forAll(myNeighbours, neighbourI)
    {
        label nodeI = myNeighbours[neighbourI];

        if (nodeI != prevNodeI)
        {
            if (weights_[nodeI] == pathValue)
            {
                route.append(nodeI);

                storeRoute
                (
                    nodeI,
                    startNodeI,
                    pathValue,
                    route
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from connections, route later
Foam::router::router
(
    const labelListList& connections,
    const List<point>& coords
)
:
    connections_(connections),
    coords_(coords),
    weights_(coords.size(), 0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::router::route(const labelList& path, const label pathValue)
{
    if (pathValue >= 0)
    {
        FatalErrorInFunction
            << "Illegal pathValue " << pathValue << exit(FatalError);
    }

    // Reset all non-allocated weights to maximum distance
    forAll(weights_, nodeI)
    {
        if (weights_[nodeI] >= 0)
        {
            weights_[nodeI] = labelMax;
        }
    }

    if (weights_[path[0]] < 0)
    {
        // Already used
        return false;
    }
    // Get weights according to distance to starting node
    setWeights(0, path[0]);

    // Check if all endPoints can be reached
    for (label leafI = 1; leafI < path.size(); leafI++)
    {
        if (weights_[path[leafI]] == labelMax)
        {
            //Info<< "Cannot route leaf from " << path[0]
            //    << " to " << path[leafI] << " of path " << path
            //    << " since there is no valid route between them" << endl;

            // Do not fix any paths but return
            return false;
        }
    }

    // Search back from all endpoints to start and fix weights
    for (label leafI = 1; leafI < path.size(); leafI++)
    {
        fixWeights
        (
            path[0],
            path[leafI],
            path[leafI],
            -1
        );

        if (leafI < path.size() - 1)
        {
            // Update distance to take new connections into account
            forAll(weights_, nodeI)
            {
                if (weights_[nodeI]== 0)
                {
                    // Include these nodes in distance calculation
                    setWeights(0, nodeI);
                }
            }
        }
    }

    // All nodes on the path will now have value 0.
    // Mark these nodes with the (negative) pathvalue
    forAll(weights_, nodeI)
    {
        if (weights_[nodeI] == 0)
        {
            weights_[nodeI] = pathValue;
        }
    }

    // Reset unallocated weights to 0
    forAll(weights_, nodeI)
    {
        if (weights_[nodeI] > 0)
        {
            weights_[nodeI] = 0;
        }
    }

    return true;
}


Foam::labelList Foam::router::getRoute(const label pathValue) const
{
    // Find a starting point
    label pathNodeI = getValue(pathValue);

    if (pathNodeI == -1)
    {
        FatalErrorInFunction
            << "No route with value " << pathValue << endl;
    }

    // Find end or start by walking
    label startNodeI = findEndNode(pathNodeI, -1, pathValue);

    // Walk to other end and store

    DynamicList<label> route(weights_.size());

    route.append(startNodeI);

    storeRoute(startNodeI, -1, pathValue, route);

    route.shrink();

    return route;
}

// ************************************************************************* //
