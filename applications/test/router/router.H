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

Class
    Foam::router

Description
    Lee's PCB routing algorithm. Construct with list of connections between
    nodes (i.e. topology) and list of coordinates of nodes (can be vector::zero)

    Use e.g.

        // Enter topology/geometry
        router cellRouter
        (
            mesh().cellCells(),
            mesh().cellCentres()
        );

        // Try to route connections one by one. Specify unique value (<0) to
        // mark path with.
        forAll(wantedConnections, i)
        {
            bool success = cellRouter.route(wantedConnections[i], -(i+1));
        }


    The coordinates are only used at the moment for diagonal preference of
    routing:

    So not:

    +A
    |
    |
    |
    |
    ------+B

    But:

    + A
    |_
      |_
        |_
          |_
            |
            + B


    Lee algo: take array with same dimensions as grid of nodes. Initialize to
    large number. Put 0 at starting point. Now recursively assign neighbours
    as current value plus one. Stop if you hit node which has smaller number.
    Phase two is where you search path with lowest value. These are assigned
    negative number so they for next route are not overwritten.

SourceFiles
    router.C

\*---------------------------------------------------------------------------*/

#ifndef router_H
#define router_H

#include "labelList.H"
#include "pointField.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class router Declaration
\*---------------------------------------------------------------------------*/

class router
{
    // Private data

        //- Connections
        const labelListList connections_;

        //- Coordinates of nodes
        const pointField coords_;

        //- Routing table
        labelList weights_;

    // Private Member Functions

        //- Return number of weights. Utility function
        label count(const label weight) const;

        //- Set distance from nodeI
        void setWeights
        (
            const label weight,
            const label nodeI
        );

        //- Finds path from nodeI to startNodeI by travelling in direction
        //  of lowest weight
        void fixWeights
        (
            const label startNodeI,
            const label endNodeI,
            const label nodeI,
            const label prevNodeI
        );

        //- Routes single path
        bool shortestPath
        (
            const labelList& path,
            const label pathValue
        );

        //- Linear search for element with weight
        label getValue(const label) const;

        //- Find node which has no neighbours with pathValue
        label findEndNode
        (
            const label startNodeI,
            const label prevNodeI,
            const label pathValue
        ) const;

        //- Append all pathValue weights to route.
        void storeRoute
        (
            const label startNodeI,
            const label prevNodeI,
            const label pathValue,
            DynamicList<label>& route
        ) const;

        //- Disallow default bitwise copy construct
        router(const router&);

        //- Disallow default bitwise assignment
        void operator=(const router&);


public:

    // Constructors

        //- Construct from connections, route later.
        router
        (
            const labelListList& connections,
            const List<point>& coords
        );


    // Member Functions

        // Access

            const labelList& weights() const
            {
                return weights_;
            }

        // Edit

            //- Find path from first element in path to all other elements
            //  Mark resulting path in weights with (negative) pathValue.
            //  Returns false and does not mark any elements if cannot route.
            bool route(const labelList& path, const label pathValue);

            //- Extract labels of route with given value.
            labelList getRoute(const label pathValue) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
