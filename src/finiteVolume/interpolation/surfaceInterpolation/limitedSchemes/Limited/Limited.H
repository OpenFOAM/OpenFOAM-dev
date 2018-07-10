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
    Foam::LimitedLimiter

Description
    Foam::LimitedLimiter

\*---------------------------------------------------------------------------*/

#ifndef Limited_H
#define Limited_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class LimitedLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimitedScheme>
class LimitedLimiter
:
    public LimitedScheme
{
    //- Lower and upper bound of the variable
    scalar lowerBound_, upperBound_;

    void checkParameters(Istream& is)
    {
        if (lowerBound_ > upperBound_)
        {
            FatalIOErrorInFunction(is)
                << "Invalid bounds.  Lower = " << lowerBound_
                << "  Upper = " << upperBound_
                << ".  Lower bound is higher than the upper bound."
                << exit(FatalIOError);
        }
    }


public:

    LimitedLimiter
    (
        const scalar lowerBound,
        const scalar upperBound,
        Istream& is
    )
    :
        LimitedScheme(is),
        lowerBound_(lowerBound),
        upperBound_(upperBound)
    {
        checkParameters(is);
    }

    LimitedLimiter(Istream& is)
    :
        LimitedScheme(is),
        lowerBound_(readScalar(is)),
        upperBound_(readScalar(is))
    {
        checkParameters(is);
    }


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
        // If not between the lower and upper bounds use upwind
        if
        (
            (faceFlux > 0 && (phiP < lowerBound_ || phiN > upperBound_))
         || (faceFlux < 0 && (phiN < lowerBound_ || phiP > upperBound_))
         /*
            phiP < lowerBound_
         || phiP > upperBound_
         || phiN < lowerBound_
         || phiN > upperBound_
         */
        )
        {
            return 0;
        }
        else
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
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
