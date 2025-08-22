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

    globalBounds = MULEScontrols.lookupOrDefault<Switch>("globalBounds", false);

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


void Foam::MULES::limitSumCorr(UPtrList<scalarField>& phiPsiCorrs)
{
    forAll(phiPsiCorrs[0], facei)
    {
        scalar sumPos = 0;
        scalar sumNeg = 0;

        forAll(phiPsiCorrs, phasei)
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

            forAll(phiPsiCorrs, phasei)
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

            forAll(phiPsiCorrs, phasei)
            {
                if (phiPsiCorrs[phasei][facei] < 0)
                {
                    phiPsiCorrs[phasei][facei] *= lambda;
                }
            }
        }
    }
}


void Foam::MULES::limitSumCorr
(
    UPtrList<scalarField>& phiPsiCorrs,
    const UPtrList<scalarField>& phiPsiFixedCorrs
)
{
    forAll(phiPsiCorrs[0], facei)
    {
        scalar sumFixed = 0;

        forAll(phiPsiFixedCorrs, i)
        {
            sumFixed += phiPsiFixedCorrs[i][facei];
        }

        scalar sumPos = 0;
        scalar sumNeg = 0;

        forAll (phiPsiCorrs, i)
        {
            if (phiPsiCorrs[i][facei] > 0)
            {
                sumPos += phiPsiCorrs[i][facei];
            }
            else
            {
                sumNeg += phiPsiCorrs[i][facei];
            }
        }

        const scalar sum = sumPos + sumNeg;

        if (sum > 0 && sumPos > vSmall)
        {
            const scalar lambda = -(sumNeg + sumFixed)/sumPos;

            forAll (phiPsiCorrs, i)
            {
                if (phiPsiCorrs[i][facei] > 0)
                {
                    phiPsiCorrs[i][facei] *= lambda;
                }
            }
        }
        else if (sum < 0 && sumNeg < -vSmall)
        {
            const scalar lambda = -(sumPos + sumFixed)/sumNeg;

            forAll (phiPsiCorrs, i)
            {
                if (phiPsiCorrs[i][facei] < 0)
                {
                    phiPsiCorrs[i][facei] *= lambda;
                }
            }
        }
    }
}


void Foam::MULES::limitSumCorr
(
    const UPtrList<const volScalarField>& psis,
    UPtrList<surfaceScalarField>& psiPhiCorrs,
    const surfaceScalarField& phi
)
{
    {
        UPtrList<scalarField> psiPhiCorrsInternal(psiPhiCorrs.size());

        forAll(psiPhiCorrsInternal, phasei)
        {
            psiPhiCorrsInternal.set(phasei, &psiPhiCorrs[phasei]);
        }

        limitSumCorr(psiPhiCorrsInternal);
    }

    const surfaceScalarField::Boundary& phibf = phi.boundaryField();

    forAll(phibf, patchi)
    {
        label nFixed = 0;

        forAll(psis, phasei)
        {
            if
            (
               !psis[phasei].boundaryField()[patchi].assignable()
             || psiPhiCorrs[phasei].boundaryField()[patchi].fixesValue()
            )
            {
                nFixed++;
            }
        }

        if (nFixed == 0)
        {
            UPtrList<scalarField> psiPhiCorrsPatch(psiPhiCorrs.size());

            forAll(psiPhiCorrsPatch, phasei)
            {
                psiPhiCorrsPatch.set
                (
                    phasei,
                    &psiPhiCorrs[phasei].boundaryFieldRef()[patchi]
                );
            }

            limitSumCorr(psiPhiCorrsPatch);
        }
        else if (nFixed < psiPhiCorrs.size())
        {
            UPtrList<scalarField> psiPhiCorrsPatch(psiPhiCorrs.size());
            UPtrList<scalarField> psiPhiCorrsFixedPatch(psiPhiCorrs.size());

            label i = 0;
            label fixedi = 0;

            forAll(psis, phasei)
            {
                if
                (
                   !psis[phasei].boundaryField()[patchi].assignable()
                 || psiPhiCorrs[phasei].boundaryField()[patchi].fixesValue()
                )
                {
                    psiPhiCorrsFixedPatch.set
                    (
                        fixedi++,
                        &psiPhiCorrs[phasei].boundaryFieldRef()[patchi]
                    );
                }
                else
                {
                    psiPhiCorrsPatch.set
                    (
                        i++,
                        &psiPhiCorrs[phasei].boundaryFieldRef()[patchi]
                    );
                }
            }

            psiPhiCorrsPatch.setSize(i);
            psiPhiCorrsFixedPatch.setSize(fixedi);

            limitSumCorr(psiPhiCorrsPatch, psiPhiCorrsFixedPatch);
        }
    }
}


void Foam::MULES::limitSum
(
    const UPtrList<const volScalarField>& psis,
    const PtrList<surfaceScalarField>& alphaPhiBDs,
    UPtrList<surfaceScalarField>& psiPhis,
    const surfaceScalarField& phi
)
{
    forAll(psiPhis, phasei)
    {
        psiPhis[phasei] == psiPhis[phasei] - alphaPhiBDs[phasei];
    }

    limitSumCorr(psis, psiPhis, phi);

    forAll(psiPhis, phasei)
    {
        psiPhis[phasei] == psiPhis[phasei] + alphaPhiBDs[phasei];
    }
}


void Foam::MULES::limitSum
(
    const UPtrList<const volScalarField>& psis,
    UPtrList<surfaceScalarField>& psiPhis,
    const surfaceScalarField& phi
)
{
    PtrList<surfaceScalarField> alphaPhiBDs(psiPhis.size());

    forAll(psiPhis, phasei)
    {
        alphaPhiBDs.set
        (
            phasei,
            upwind<scalar>(phi.mesh(), phi).flux(psis[phasei])
        );
    }

    limitSum(psis, alphaPhiBDs, psiPhis, phi);
}


// ************************************************************************* //
