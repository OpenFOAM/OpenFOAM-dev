/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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

\*---------------------------------------------------------------------------*/

#include "MULES.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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

        scalar sum = sumPos + sumNeg;

        if (sum > 0 && sumPos > vSmall)
        {
            scalar lambda = -sumNeg/sumPos;

            for (int phasei=0; phasei<phiPsiCorrs.size(); phasei++)
            {
                if (phiPsiCorrs[phasei][facei] > 0)
                {
                    phiPsiCorrs[phasei][facei] *= lambda;
                }
            }
        }
        else if (sum < 0 && sumNeg < -vSmall)
        {
            scalar lambda = -sumPos/sumNeg;

            for (int phasei=0; phasei<phiPsiCorrs.size(); phasei++)
            {
                if (phiPsiCorrs[phasei][facei] < 0)
                {
                    phiPsiCorrs[phasei][facei] *= lambda;
                }
            }
        }
    }
}


void Foam::MULES::limitSum
(
    const UPtrList<const scalarField>& alphas,
    UPtrList<scalarField>& phiPsiCorrs,
    const labelHashSet& fixed
)
{
    labelHashSet notFixed(identity(phiPsiCorrs.size()));
    notFixed -= fixed;

    forAll(phiPsiCorrs[0], facei)
    {
        scalar alphaNotFixed = 0, corrNotFixed = 0;
        forAllConstIter(labelHashSet, notFixed, iter)
        {
            alphaNotFixed += alphas[iter.key()][facei];
            corrNotFixed += phiPsiCorrs[iter.key()][facei];
        }

        scalar corrFixed = 0;
        forAllConstIter(labelHashSet, fixed, iter)
        {
            corrFixed += phiPsiCorrs[iter.key()][facei];
        }

        const scalar sumCorr = corrNotFixed + corrFixed;

        const scalar lambda = - sumCorr/alphaNotFixed;

        forAllConstIter(labelHashSet, notFixed, iter)
        {
            phiPsiCorrs[iter.key()][facei] += lambda*alphas[iter.key()][facei];
        }
    }
}


// ************************************************************************* //
