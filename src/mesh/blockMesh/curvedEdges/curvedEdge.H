/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::curvedEdge

Description
    Define a curved edge that is parameterized for 0<lambda<1
    between the start and end point.

SourceFiles
    curvedEdge.C

\*---------------------------------------------------------------------------*/

#ifndef curvedEdges_H
#define curvedEdges_H

#include "edge.H"
#include "pointField.H"
#include "typeInfo.H"
#include "HashTable.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class curvedEdge Declaration
\*---------------------------------------------------------------------------*/

class curvedEdge
{
protected:

    // Protected data

        const pointField& points_;
        const label start_;
        const label end_;

    // Protected Member Functions

        //- Return a complete point field by appending the start/end points
        //  to the given list
        static pointField appendEndPoints
        (
            const pointField&,
            const label start,
            const label end,
            const pointField& otherKnots
        );


public:

    //- Runtime type information
    TypeName("curvedEdge");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            curvedEdge,
            Istream,
            (
                const pointField& points,
                Istream& is
            ),
            (points, is)
        );


    // Constructors

        //- Construct from components
        curvedEdge
        (
            const pointField& points,
            const label start,
            const label end
        );

        //- Construct from Istream setting pointsList
        curvedEdge(const pointField&, Istream&);

        //- Copy construct
        curvedEdge(const curvedEdge&);

        //- Clone function
        virtual autoPtr<curvedEdge> clone() const;

        //- New function which constructs and returns pointer to a curvedEdge
        static autoPtr<curvedEdge> New(const pointField&, Istream&);


    //- Destructor
    virtual ~curvedEdge(){}


    // Member Functions

        //- Return label of start point
        inline label start() const;

        //- Return label of end point
        inline label end() const;

        //- Compare the given start and end points with this curve
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const curvedEdge&) const;

        //- Compare the given start and end points with this curve
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const edge&) const;

        //- Compare the given start and end points with this curve
        //  Return:
        //  -  0: different
        //  - +1: identical
        //  - -1: same edge, but different orientation
        inline int compare(const label start, const label end) const;

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        virtual point position(const scalar) const = 0;

        //- Return the length of the curve
        virtual scalar length() const = 0;


    // Member operators

        void operator=(const curvedEdge&);


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const curvedEdge&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "curvedEdgeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
