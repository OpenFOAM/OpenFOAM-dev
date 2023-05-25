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
    tmp<volScalarField::Internal>& Su,
    tmp<volScalarField::Internal>& Sp
)
{
    if (divergent())
    {
        // Phase change alpha1 source
        const fvScalarMatrix alphaSup(fvModels().source(alpha1));

        Su = alphaSup.Su();
        Sp = alphaSup.Sp();
    }
}


// ************************************************************************* //
