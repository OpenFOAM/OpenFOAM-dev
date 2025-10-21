/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2025 OpenFOAM Foundation
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

#include "ubMomentumTransportModel.H"
#include "fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

defineTypeNameAndDebug(ubMomentumTransportModel, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ubMomentumTransportModel::ubMomentumTransportModel
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const compressible::momentumTransportModel& mixtureMomentumTransport
)
:
    RASModel<phaseCompressibleMomentumTransportModel>
    (
        RASModels::ubMomentumTransportModel::typeName,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
    mixtureMomentumTransport_(mixtureMomentumTransport)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool ubMomentumTransportModel::read()
{
    return true;
}


tmp<fvVectorMatrix> ubMomentumTransportModel::divDevTau(volVectorField& U) const
{
    NotImplemented;
    return tmp<fvVectorMatrix>(nullptr);
}


tmp<fvVectorMatrix> ubMomentumTransportModel::divDevTau
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    NotImplemented;
    return tmp<fvVectorMatrix>(nullptr);
}


void ubMomentumTransportModel::correct()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
