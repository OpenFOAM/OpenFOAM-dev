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

#include "compressibleMomentumTransportModel.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleMomentumTransportModel::compressibleMomentumTransportModel
(
    const word& type,
    const geometricOneField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
:
    momentumTransportModel
    (
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
    alpha_(alpha),
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::compressibleMomentumTransportModel>
Foam::compressibleMomentumTransportModel::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const viscosity& viscosity
)
{
    return momentumTransportModel::New<compressibleMomentumTransportModel>
    (
        geometricOneField(),
        rho,
        U,
        phi,
        phi,
        viscosity
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::compressibleMomentumTransportModel::phi() const
{
    if (phi_.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        return phi_;
    }
    else
    {
        return phi_/fvc::interpolate(rho_);
    }
}


// ************************************************************************* //
