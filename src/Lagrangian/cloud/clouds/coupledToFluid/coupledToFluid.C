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
#include "fluidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(coupledToFluid, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::clouds::coupledToFluid::getRhocVf(const word& phaseName) const
{
    const word rhocName = IOobject::groupName("rho", phaseName);

    if (mesh_.poly().foundObject<volScalarField>(rhocName))
    {
        return mesh_.poly().lookupObject<volScalarField>(rhocName);
    }

    const word thermocName =
        IOobject::groupName(physicalProperties::typeName, phaseName);

    if (mesh_.poly().foundObject<basicThermo>(thermocName))
    {
        return mesh_.poly().lookupObject<basicThermo>(thermocName).rho();
    }

    FatalErrorInFunction
        << "Could not determine the carrier density"
        << exit(FatalError);

    return tmp<volScalarField>(nullptr);
}


const Foam::volScalarField&
Foam::clouds::coupledToFluid::getMucVf(const word& phaseName) const
{
    const word mucName = IOobject::groupName("mu", phaseName);

    if (mesh_.poly().foundObject<volScalarField>(mucName))
    {
        return mesh_.poly().lookupObject<volScalarField>(mucName);
    }

    const word thermocName =
        IOobject::groupName(physicalProperties::typeName, phaseName);

    if (mesh_.poly().foundObject<fluidThermo>(thermocName))
    {
        return mesh_.poly().lookupObject<fluidThermo>(thermocName).mu();
    }

    return NullObjectRef<volScalarField>();
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::clouds::coupledToFluid::calcNuc
(
    const LagrangianModelRef& model,
    const LagrangianSubMesh& subMesh
) const
{
    return muc(model, subMesh)/rhoc(model, subMesh);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::clouds::coupledToFluid::updateCarrier()
{
    coupled::updateCarrier();

    if (trhocVf_.isTmp())
    {
        trhocVf_.ref() = getRhocVf(carriedCloud_.carrierPhaseName());
    }

    if (trhocPhaseVf_.isTmp())
    {
        trhocPhaseVf_.ref() = getRhocVf(carriedCloud_.phaseName());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::coupledToFluid::coupledToFluid
(
    const cloud& c,
    const carried& carriedCloud
)
:
    coupled(c, carriedCloud),
    mesh_(c.mesh()),
    carriedCloud_(carriedCloud),
    trhocVf_(getRhocVf(carriedCloud.carrierPhaseName())),
    trhocPhaseVf_
    (
        carriedCloud.hasPhase()
      ? getRhocVf(carriedCloud.phaseName())
      : tmp<volScalarField>(NullObjectRef<volScalarField>())
    ),
    mucVf_(getMucVf(carriedCloud.carrierPhaseName())),
    rhoc(carriedCloud.carrierField<scalar>(trhocVf_())),
    rhocPhase
    (
        carriedCloud.hasPhase()
      ? carriedCloud.carrierField<scalar>(trhocPhaseVf_())
      : carriedCloud.noCarrierField<scalar>("rho", "density", true)
    ),
    muc
    (
        isNull(mucVf_)
      ? c.derivedField<scalar>
        (
            [&]
            (
                const LagrangianModelRef& model,
                const LagrangianSubMesh& subMesh
            )
            {
                return rhoc(model, subMesh)*nuc(model, subMesh);
            }
        )
      : carriedCloud.carrierField<scalar>(mucVf_)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupledToFluid::~coupledToFluid()
{}


// ************************************************************************* //
