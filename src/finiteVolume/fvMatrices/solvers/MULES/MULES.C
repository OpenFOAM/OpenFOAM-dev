/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "MULES.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MULES, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::MULES::control::control(const dictionary& dict)
{
    read(dict);
}


void Foam::MULES::control::read(const dictionary& dict)
{
    const dictionary& MULEScontrols = dict.optionalSubDict("MULES");

    smoothingCoeff = MULEScontrols.lookupOrDefault<scalar>("smoothingCoeff", 0);

    extremaCoeff = MULEScontrols.lookupOrDefault<scalar>("extremaCoeff", 0);

    boundaryExtremaCoeff =
        MULEScontrols.lookupOrDefault<scalar>("boundaryExtremaCoeff", 0);

    if (dict.found("MULES"))
    {
        nIter = MULEScontrols.lookupOrDefault<label>("nIter", 3);

        tol =
            nIter == 1
          ? 0
          : MULEScontrols.lookupOrDefault<scalar>("tolerance", 0);
    }
    else
    {
        nIter = MULEScontrols.lookupOrDefault<label>("nLimiterIter", 3);

        tol =
            nIter == 1
          ? 0
          : MULEScontrols.lookupOrDefault<scalar>("MULEStolerance", 0);
    }
}


void Foam::MULES::limitSum(UPtrList<scalarField>& phiPsiCorrs)
{
    forAll(phiPsiCorrs[0], facei)
    {
        scalar sumPos = 0;
        scalar sumNeg = 0;

        for (int phasei=0; phasei<phiPsiCorrs.size(); phasei++)
        {
            if (phiPsiCorrs[phasei][facei] > 0)
            {
                sumPos += phiPsiCorrs[phasei][facei];
            }
            else
            {
                sumNeg += phiPsiCorrs[phasei][facei];
            }
        }

        const scalar sum = sumPos + sumNeg;

        if (sum > 0 && sumPos > vSmall)
        {
            const scalar lambda = -sumNeg/sumPos;

            forAll (phiPsiCorrs, phasei)
            {
                if (phiPsiCorrs[phasei][facei] > 0)
                {
                    phiPsiCorrs[phasei][facei] *= lambda;
                }
            }
        }
        else if (sum < 0 && sumNeg < -vSmall)
        {
            const scalar lambda = -sumPos/sumNeg;

            forAll (phiPsiCorrs, phasei)
            {
                if (phiPsiCorrs[phasei][facei] < 0)
                {
                    phiPsiCorrs[phasei][facei] *= lambda;
                }
            }
        }
    }
}


// ************************************************************************* //
