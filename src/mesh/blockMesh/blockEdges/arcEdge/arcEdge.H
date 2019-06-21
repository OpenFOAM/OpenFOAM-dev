/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::blockEdges::arcEdge

Description
    Defines the arcEdge of a circle in terms of 3 points on its circumference

SourceFiles
    arcEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_arcEdge_H
#define blockEdges_arcEdge_H

#include "blockEdge.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                           Class arcEdge Declaration
\*---------------------------------------------------------------------------*/

class arcEdge
:
    public blockEdge
{
    // Private Data

        point p1_, p2_, p3_;
        cylindricalCS cs_;
        scalar angle_;
        scalar radius_;

    // Private Member Functions

        //- Calculate the coordinate system, angle and radius
        cylindricalCS calcAngle();


public:

    //- Runtime type information
    TypeName("arc");


    // Constructors

        //- Construct from components
        arcEdge
        (
            const pointField& points,
            const label start, const label end,
            const point& pMid
        );

        //- Construct from Istream setting pointsList
        arcEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            const pointField& points,
            Istream&
        );

        //- Disallow default bitwise copy construction
        arcEdge(const arcEdge&) = delete;


    //- Destructor
    virtual ~arcEdge()
    {}


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar) const;

        //- Return the length of the curve
        scalar length() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const arcEdge&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End of namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
