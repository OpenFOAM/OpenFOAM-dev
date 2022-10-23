/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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

#include "incompressibleMomentumTransportModel.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleMomentumTransportModel::incompressibleMomentumTransportModel
(
    const word& type,
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    momentumTransportModel(U, alphaRhoPhi, phi, viscosity),
    alpha_(alpha),
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::incompressibleMomentumTransportModel>
Foam::incompressibleMomentumTransportModel::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
{
    return momentumTransportModel::New<incompressibleMomentumTransportModel>
    (
        geometricOneField(),
        geometricOneField(),
        U,
        phi,
        phi,
        viscosity
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField>
Foam::incompressibleMomentumTransportModel::devSigma() const
{
    return devTau();
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::incompressibleMomentumTransportModel::divDevSigma(volVectorField& U) const
{
    return divDevTau(U);
}


Foam::tmp<Foam::volSymmTensorField>
Foam::incompressibleMomentumTransportModel::devTau() const
{
    NotImplemented;
    return devSigma();
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::incompressibleMomentumTransportModel::divDevTau
(
    volVectorField& U
) const
{
    NotImplemented;
    return divDevSigma(U);
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::incompressibleMomentumTransportModel::divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    NotImplemented;
    return divDevSigma(U);
}


// ************************************************************************* //
