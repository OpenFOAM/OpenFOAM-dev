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
    Foam::linearFitPolynomial

Description
    Linear polynomial for interpolation fitting.

    Can be used with the CentredFit scheme to create a linear surface
    interpolation scheme

\*---------------------------------------------------------------------------*/

#ifndef linearFitPolynomial_H
#define linearFitPolynomial_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class linearFitPolynomial Declaration
\*---------------------------------------------------------------------------*/

class linearFitPolynomial
{
public:

    // Member Functions

        static label nTerms(const direction dim)
        {
            return
            (
                dim == 1 ? 2 :
                dim == 2 ? 3 :
                dim == 3 ? 4 : 0
            );
        }

        static void addCoeffs
        (
            scalar* coeffs,
            const vector& d,
            const scalar weight,
            const direction dim
        )
        {
            label curIdx = 0;

            coeffs[curIdx++] = weight;
            coeffs[curIdx++] = weight*d.x();

            if (dim >= 2)
            {
                coeffs[curIdx++] = weight*d.y();
            }
            if (dim == 3)
            {
                coeffs[curIdx++] = weight*d.z();
            }
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
