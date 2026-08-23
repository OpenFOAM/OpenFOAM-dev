/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025-2026 OpenFOAM Foundation
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

#include "coupledToFluid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(coupledToFluid, 0);
}
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::clouds::coupledToFluid::updateCarrier()
{
    coupled::updateCarrier();

    if (autoPtr<carrierDynamicViscosity>::valid())
    {
        autoPtr<carrierDynamicViscosity>::operator()().correct();
    }

    if (rhocPhasePtr_.valid())
    {
        rhocPhasePtr_.operator()().correct();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::coupledToFluid::coupledToFluid
(
    const cloud& c,
    const carried& carriedCloud,
    const carrierDynamicViscosity& mucOrNull,
    const fluidSurfaceThermo& thermocsOrNull,
    const carrierDensity& rhocPhaseOrNull
)
:
    autoPtr<carrierDynamicViscosity>
    (
        notNull(mucOrNull)
      ? nullptr
      : carrierDynamicViscosity::New
        (
            c,
            carriedCloud,
            carriedCloud.carrierPhaseName()
        ).ptr()
    ),
    coupled
    (
        c,
        carriedCloud,
        notNull(mucOrNull)
      ? mucOrNull
      : autoPtr<carrierDynamicViscosity>::operator()(),
        thermocsOrNull
    ),
    muc_
    (
        notNull(mucOrNull)
      ? mucOrNull
      : autoPtr<carrierDynamicViscosity>::operator()()
    ),
    thermocsOrNull_(thermocsOrNull),
    rhocPhasePtr_
    (
        !carriedCloud.hasPhase()
      ? nullptr
      : notNull(rhocPhaseOrNull)
      ? nullptr
      : carrierDensity::New
        (
            c,
            carriedCloud,
            carriedCloud.phaseName()
        ).ptr()
    ),
    rhoc(muc_.rho()),
    rhocPhase
    (
        !carriedCloud.hasPhase()
      ? carriedCloud.noCarrierField<scalar>("rho", "density", true)
      : notNull(rhocPhaseOrNull)
      ? rhocPhaseOrNull.rho()
      : rhocPhasePtr_->rho()
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupledToFluid::~coupledToFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::LagrangianSubScalarSubField> Foam::clouds::coupledToFluid::muc
(
    const LagrangianSubMesh& subMesh
) const
{
    return
        notNull(thermocsOrNull_)
      ? thermocsOrNull_.muNear(subMesh)
      : tmp<LagrangianSubScalarSubField>(muc_.mu().ref(subMesh));
}


// ************************************************************************* //
