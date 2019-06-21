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
    Foam::blockEdges::polyLineEdge

Description
    A blockEdge defined in terms of a series of straight line segments.

SourceFiles
    polyLineEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_polyLineEdge_H
#define blockEdges_polyLineEdge_H

#include "blockEdge.H"
#include "polyLine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                        Class polyLineEdge Declaration
\*---------------------------------------------------------------------------*/

class polyLineEdge
:
    public blockEdge,
    public polyLine
{
public:

    // Static Data Members
    TypeName("polyLine");


    // Constructors

        //- Construct from components
        polyLineEdge
        (
            const pointField&,
            const label start,
            const label end,
            const pointField& otherPoints
        );

        //- Construct from Istream
        polyLineEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            const pointField&,
            Istream&
        );

        //- Disallow default bitwise copy construction
        polyLineEdge(const polyLineEdge&) = delete;


        //- Destructor
        virtual ~polyLineEdge();


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar lambda) const;

        //- Return the length of the curve
        scalar length() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const polyLineEdge&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End of namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
