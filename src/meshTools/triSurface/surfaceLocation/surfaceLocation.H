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
    Foam::surfaceLocation

Description
    Contains information about location on a triSurface:
    - pointIndexHit:
        - location
        - bool: hit/miss
        - index (of triangle/point/edge)
    - elementType():
        - what index above relates to. In triangle::proxType
    - triangle():
        - last known triangle

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef surfaceLocation_H
#define surfaceLocation_H

#include "pointIndexHit.H"
#include "triPointRef.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurface;

/*---------------------------------------------------------------------------*\
                           Class surfaceLocation Declaration
\*---------------------------------------------------------------------------*/

class surfaceLocation
:
    public pointIndexHit
{
    // Private data

        triPointRef::proxType elementType_;

        label triangle_;

public:


    // Constructors

        //- Construct null
        surfaceLocation()
        :
            pointIndexHit(),
            elementType_(triPointRef::NONE),
            triangle_(-1)
        {}

        //- Construct from components
        surfaceLocation
        (
            const pointIndexHit& pih,
            const triPointRef::proxType elementType,
            const label triangle
        )
        :
            pointIndexHit(pih),
            elementType_(elementType),
            triangle_(triangle)
        {}

        //- Construct from Istream
        surfaceLocation(Istream& is)
        :
            pointIndexHit(is),
            elementType_(triPointRef::proxType(readLabel(is))),
            triangle_(readLabel(is))
        {}


    // Member Functions

        triPointRef::proxType& elementType()
        {
            return elementType_;
        }

        triPointRef::proxType elementType() const
        {
            return elementType_;
        }

        label& triangle()
        {
            return triangle_;
        }

        label triangle() const
        {
            return triangle_;
        }

        //- Normal. Approximate for points.
        vector normal(const triSurface& s) const;

        //- Return info proxy.
        //  Used to print token information to a stream
        InfoProxy<surfaceLocation> info() const
        {
            return *this;
        }

        //- Write info to os
        void write(Ostream& os, const triSurface& s) const;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, surfaceLocation& sl);

        friend Ostream& operator<<(Ostream& os, const surfaceLocation& sl);

        friend Ostream& operator<<
        (
            Ostream&,
            const InfoProxy<surfaceLocation>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
