/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "compressibleVoF.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleVoF::alphaSuSp
(
    tmp<volScalarField::Internal>& tSu,
    tmp<volScalarField::Internal>& tSp,
    const dictionary& alphaControls
)
{
    const scalar vDotResidualAlpha
    (
        alphaControls.lookupOrDefault("vDotResidualAlpha", 1e-4)
    );

    const dimensionedScalar Szero(vDot.dimensions(), 0);

    tSp = volScalarField::Internal::New("Sp", mesh, Szero);
    tSu = volScalarField::Internal::New("Su", mesh, Szero);

    volScalarField::Internal& Sp = tSp.ref();
    volScalarField::Internal& Su = tSu.ref();

    if (fvModels().addsSupToField(mixture.rho1().name()))
    {
        const volScalarField::Internal alpha2ByRho1(alpha2()/mixture.rho1()());
        const fvScalarMatrix alphaRho1Sup
        (
            fvModels().sourceProxy(alpha1, mixture.rho1(), alpha1)
        );

        Su += alpha2ByRho1*alphaRho1Sup.Su();
        Sp += alpha2ByRho1*alphaRho1Sup.Sp();
    }

    if (fvModels().addsSupToField(mixture.rho2().name()))
    {
        const volScalarField::Internal alpha1ByRho2(alpha1()/mixture.rho2()());
        const fvScalarMatrix alphaRho2Sup
        (
            fvModels().sourceProxy(alpha2, mixture.rho2(), alpha2)
        );

        Su -= alpha1ByRho2*(alphaRho2Sup.Su() + alphaRho2Sup.Sp());
        Sp += alpha1ByRho2*alphaRho2Sup.Sp();
    }

    forAll(vDot, celli)
    {
        if (vDot[celli] > 0.0)
        {
            Sp[celli] -=
                vDot[celli]/max(1.0 - alpha1[celli], vDotResidualAlpha);
            Su[celli] +=
                vDot[celli]/max(1.0 - alpha1[celli], vDotResidualAlpha);
        }
        else if (vDot[celli] < 0.0)
        {
            Sp[celli] += vDot[celli]/max(alpha1[celli], vDotResidualAlpha);
        }
    }
}


// ************************************************************************* //
