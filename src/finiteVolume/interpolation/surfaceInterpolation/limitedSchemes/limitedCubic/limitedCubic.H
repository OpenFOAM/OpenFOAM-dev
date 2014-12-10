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
    Foam::limitedCubicLimiter

Description
    Class with limiter function which returns the limiter for the
    TVD limited centred-cubic differencing scheme based on r obtained from
    the LimiterFunc class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    limitedCubic.C

\*---------------------------------------------------------------------------*/

#ifndef limitedCubic_H
#define limitedCubic_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class limitedCubicLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class limitedCubicLimiter
:
    public LimiterFunc
{
    scalar k_;
    scalar twoByk_;

public:

    limitedCubicLimiter(Istream& is)
    :
        k_(readScalar(is))
    {
        if (k_ < 0 || k_ > 1)
        {
            FatalIOErrorIn("limitedCubicLimiter(Istream& is)", is)
                << "coefficient = " << k_
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }

        // Avoid the /0 when k_ = 0
        twoByk_ = 2.0/max(k_, SMALL);
    }

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        scalar twor = twoByk_*LimiterFunc::r
        (
            faceFlux, phiP, phiN, gradcP, gradcN, d
        );

        scalar phiU;

        if (faceFlux > 0)
        {
            phiU = phiP;
        }
        else
        {
            phiU = phiN;
        }

        // Calculate the face value using cubic interpolation
        scalar phif =
            cdWeight*(phiP - 0.25*(d & gradcN))
          + (1 - cdWeight)*(phiN + 0.25*(d & gradcP));

        scalar phiCD = cdWeight*phiP + (1 - cdWeight)*phiN;

        // Calculate the effective limiter for the cubic interpolation
        scalar cubicLimiter = (phif - phiU)/stabilise(phiCD - phiU, SMALL);

        /*
        if (twor < 0.05)
        {
            cubicLimiter = twor;
        }

        return max(min(cubicLimiter, 2), 0);
        */

        // Limit the limiter to obey the TVD constraint
        return max(min(min(twor, cubicLimiter), 2), 0);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
