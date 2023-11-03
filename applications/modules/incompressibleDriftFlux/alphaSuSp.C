/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "incompressibleDriftFlux.H"
#include "fvcFlux.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::incompressibleDriftFlux::alphaPhi
(
    const surfaceScalarField& phi,
    const volScalarField& alpha,
    const dictionary& alphaControls
)
{
    return fvc::flux
    (
        phi + fvc::flux(relativeVelocity->Udm()),
        alpha,
        divAlphaName
    );
}


void Foam::solvers::incompressibleDriftFlux::alphaSuSp
(
    tmp<volScalarField::Internal>& tSu,
    tmp<volScalarField::Internal>& tSp,
    const dictionary& alphaControls
)
{
    if (!divergent()) return;

    const dimensionedScalar Szero(dimless/dimTime, 0);

    tSp = volScalarField::Internal::New("Sp", mesh, Szero);
    tSu = volScalarField::Internal::New("Su", mesh, Szero);

    volScalarField::Internal& Sp = tSp.ref();
    volScalarField::Internal& Su = tSu.ref();

    if (fvModels().addsSupToField(alpha1.name()))
    {
        const fvScalarMatrix alpha1Sup(fvModels().source(alpha1));

        Su += alpha2()*alpha1Sup.Su();
        Sp += alpha2()*alpha1Sup.Sp();
    }

    if (fvModels().addsSupToField(alpha2.name()))
    {
        const fvScalarMatrix alpha2Sup(fvModels().source(alpha2));

        Su -= alpha1()*(alpha2Sup.Su() + alpha2Sup.Sp());
        Sp += alpha1()*alpha2Sup.Sp();
    }
}


// ************************************************************************* //
