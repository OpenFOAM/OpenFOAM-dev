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
    Foam::filteredLinear2VLimiter

Description
    Class to generate weighting factors for the filteredLinear2V
    interpolation scheme.

    The aim is to remove high-frequency modes with "staggering"
    characteristics from vector fields by comparing the face gradient in
    the direction of maximum gradient with both neighbouring cell gradients
    and introduce small amounts of upwind in order to damp these modes.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    filteredLinear2V.C

\*---------------------------------------------------------------------------*/

#ifndef filteredLinear2V_H
#define filteredLinear2V_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class filteredLinear2VLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class filteredLinear2VLimiter
:
    public LimiterFunc
{
    // Private Data

        // Scaling coefficient for the gradient ratio,
        // 0 = linear
        // 1 = fully limited
        scalar k_;

        // Maximum allowed overshoot/undershoot relative to the difference
        // across the face.
        // On input:
        //     0 = no overshoot/undershoot
        //     1 = overshoot/undershoot equal to the difference across the face
        // Note: After input 1 is added to l_
        scalar l_;


public:

    filteredLinear2VLimiter(Istream& is)
    :
        k_(readScalar(is)),
        l_(readScalar(is))
    {
        if (k_ < 0 || k_ > 1)
        {
            FatalIOErrorInFunction(is)
                << "coefficient = " << k_
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }

        if (l_ < 0 || l_ > 1)
        {
            FatalIOErrorInFunction(is)
                << "coefficient = " << l_
                << " should be >= 0 and <= 1"
                << exit(FatalIOError);
        }

        l_ += 1.0;
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
        // Difference across face
        vector dfV = phiN - phiP;

        // Scalar difference across the face
        // in the direction in which the difference is largest
        scalar df = dfV & dfV;

        // Twice differences across face-neighbour cells
        // in the direction in which the face-difference is largest
        scalar tdcP = 2*(dfV & (d & gradcP));
        scalar tdcN = 2*(dfV & (d & gradcN));

        // Calculate the limiter according to the sign of the face difference
        scalar limiter;

        if (df > 0)
        {
            limiter = l_
                - k_*min(max(df - tdcP, 0), max(df - tdcN, 0))
                /(max(mag(df), max(mag(tdcP), mag(tdcN))) + vSmall);
        }
        else
        {
            limiter = l_
                - k_*min(max(tdcP - df, 0), max(tdcN - df, 0))
                /(max(mag(df), max(mag(tdcP), mag(tdcN))) + vSmall);
        }

        // Limit the limiter between linear and upwind
        return max(min(limiter, 1), 0);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
