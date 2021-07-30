/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

#include "phaseIncompressibleMomentumTransportModel.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseIncompressibleMomentumTransportModel, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseIncompressibleMomentumTransportModel::
phaseIncompressibleMomentumTransportModel
(
    const word& type,
    const volScalarField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    incompressibleMomentumTransportModel
    (
        type,
        geometricOneField(),
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
    alpha_(alpha)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseIncompressibleMomentumTransportModel>
Foam::phaseIncompressibleMomentumTransportModel::New
(
    const volScalarField& alpha,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
{
    return momentumTransportModel::New
    <
        phaseIncompressibleMomentumTransportModel
    >
    (
        alpha,
        geometricOneField(),
        U,
        alphaRhoPhi,
        phi,
        viscosity
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::phaseIncompressibleMomentumTransportModel::pPrime() const
{
    return volScalarField::New
    (
        IOobject::groupName("pPrime", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(dimPressure, 0)
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::phaseIncompressibleMomentumTransportModel::pPrimef() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName("pPrimef", this->alphaRhoPhi_.group()),
        this->mesh_,
        dimensionedScalar(dimPressure, 0)
    );
}


Foam::tmp<Foam::volSymmTensorField>
Foam::phaseIncompressibleMomentumTransportModel::devSigma() const
{
    return devTau();
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::phaseIncompressibleMomentumTransportModel::divDevSigma
(
    volVectorField& U
) const
{
    return divDevTau(U);
}


Foam::tmp<Foam::volSymmTensorField>
Foam::phaseIncompressibleMomentumTransportModel::
devTau() const
{
    NotImplemented;

    return devSigma();
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::phaseIncompressibleMomentumTransportModel::divDevTau
(
    volVectorField& U
) const
{
    NotImplemented;

    return divDevSigma(U);
}


// ************************************************************************* //
