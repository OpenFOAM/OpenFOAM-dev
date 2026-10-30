/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::FilteredLimiter

Description
    Foam::FilteredLimiter

\*---------------------------------------------------------------------------*/

#ifndef Filtered_H
#define Filtered_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class FilteredLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimitedScheme>
class FilteredLimiter
:
    public LimitedScheme
{
    //- Lower limit of the variable below which the limiter is not applied
    scalar lowerLimit_;


public:

    FilteredLimiter
    (
        const scalar lowerLimit,
        Istream& is
    )
    :
        LimitedScheme(is),
        lowerLimit_(lowerLimit)
    {}

    FilteredLimiter(Istream& is)
    :
        LimitedScheme(is),
        lowerLimit_(readScalar(is))
    {}


    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const scalar phiP,
        const scalar phiN,
        const vector& gradcP,
        const vector& gradcN,
        const vector& d
    ) const
    {
        if
        (
            mag(phiP) > lowerLimit_
         && mag(phiN) > lowerLimit_
         && mag(d & gradcP) > lowerLimit_
         && mag(d & gradcN) > lowerLimit_
        )
        {
            return LimitedScheme::limiter
            (
                cdWeight,
                faceFlux,
                phiP,
                phiN,
                gradcP,
                gradcN,
                d
            );
        }
        else
        {
            return LimitedScheme::maxLimit;
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
