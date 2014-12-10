/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::NVDVTVDV

Description
    Foam::NVDVTVDV

\*---------------------------------------------------------------------------*/

#ifndef NVDVTVDV_H
#define NVDVTVDV_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class LimitedScheme Declaration
\*---------------------------------------------------------------------------*/

class NVDVTVDV
{

public:

    typedef vector phiType;
    typedef tensor gradPhiType;

    // Null Constructor

        NVDVTVDV()
        {}


    // Member Functions

        scalar phict
        (
            const scalar faceFlux,
            const vector& phiP,
            const vector& phiN,
            const tensor& gradcP,
            const tensor& gradcN,
            const vector& d
        ) const
        {
            vector gradfV = phiN - phiP;
            scalar gradf = gradfV & gradfV;

            scalar gradcf;

            if (faceFlux > 0)
            {
                gradcf = gradfV & (d & gradcP);
            }
            else
            {
                gradcf = gradfV & (d & gradcN);
            }

            if (mag(gradf) >= 1000*mag(gradcf))
            {
                return 1 - 0.5*1000*sign(gradcf)*sign(gradf);
            }
            else
            {
                return 1 - 0.5*gradf/gradcf;
            }
        }


        scalar r
        (
            const scalar faceFlux,
            const vector& phiP,
            const vector& phiN,
            const tensor& gradcP,
            const tensor& gradcN,
            const vector& d
        ) const
        {
            vector gradfV = phiN - phiP;
            scalar gradf = gradfV & gradfV;

            scalar gradcf;

            if (faceFlux > 0)
            {
                gradcf = gradfV & (d & gradcP);
            }
            else
            {
                gradcf = gradfV & (d & gradcN);
            }

            if (mag(gradcf) >= 1000*mag(gradf))
            {
                return 2*1000*sign(gradcf)*sign(gradf) - 1;
            }
            else
            {
                return 2*(gradcf/gradf) - 1;
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
