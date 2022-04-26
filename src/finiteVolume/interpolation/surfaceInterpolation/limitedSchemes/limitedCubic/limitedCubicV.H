/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
    Foam::limitedCubicVLimiter

Description
    Class with limiter function which returns the limiter for the
    limitedCubicV interpolation scheme based on r obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    limitedCubicV.C

\*---------------------------------------------------------------------------*/

#ifndef limitedCubicV_H
#define limitedCubicV_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class limitedCubicVLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class limitedCubicVLimiter
:
    public LimiterFunc
{
    scalar k_;
    scalar twoByk_;

public:

    limitedCubicVLimiter(Istream& is)
    :
        k_(readScalar(is))
    {
        if (k_ < 0 || k_ > 1)
        {
            FatalIOErrorInFunction(is)
                << "coefficient = " << k_
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }

        // Avoid the /0 when k_ = 0
        twoByk_ = 2.0/max(k_, small);
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

        vector fV = cdWeight*phiP + (1.0 - cdWeight)*phiN;

        scalar fVphiP = fV & phiP;
        scalar fVphiN = fV & phiN;

        scalar fVphiU;

        if (faceFlux > 0)
        {
            fVphiU = fVphiP;
        }
        else
        {
            fVphiU = fVphiN;
        }

        // Calculate the face value using cubic interpolation
        scalar fVphif =
            cdWeight*(fVphiP - 0.25*(fV & (d & gradcN)))
          + (1 - cdWeight)*(fVphiN + 0.25*(fV & (d & gradcP)));

        scalar fVphiCD = cdWeight*fVphiP + (1 - cdWeight)*fVphiN;

        // Calculate the effective limiter for the cubic interpolation
        scalar cubicLimiter =
            (fVphif - fVphiU)/stabilise(fVphiCD - fVphiU, small);

        // Limit the limiter to obey the TVD constraint
        return max(min(min(twor, cubicLimiter), 2), 0);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
