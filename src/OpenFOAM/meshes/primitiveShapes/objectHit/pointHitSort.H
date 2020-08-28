/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::pointHitSort

Description
    Container for sorting intersections

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef pointHitSort_H
#define pointHitSort_H

#include "pointHit.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                        Class pointHitSort Declaration
\*---------------------------------------------------------------------------*/

class pointHitSort
{
    // Private Data

        //- intersection
        pointHit inter_;

        //- Original index
        label index_;

public:

    // Constructors

        //- Construct null
        pointHitSort()
        :
            inter_(false, Zero, great, false),
            index_(-1)
        {}

        //- Construct from intersection, index
        pointHitSort(const pointHit& inter, const label index)
        :
            inter_(inter),
            index_(index)
        {}


    // Member Functions

        const pointHit& inter() const
        {
            return inter_;
        }

        label index() const
        {
            return index_;
        }

    // Member Operators

        bool operator==(const pointHitSort& rhs) const
        {
            return inter_.distance() == rhs.inter().distance();
        }

        bool operator<(const pointHitSort& rhs) const
        {
            return inter_.distance() < rhs.inter().distance();
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
