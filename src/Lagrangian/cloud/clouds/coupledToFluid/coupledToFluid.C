/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

    if (mesh_.mesh().foundObject<volScalarField>(rhocName))
    {
        return mesh_.mesh().lookupObject<volScalarField>(rhocName);
    }

    const word thermocName =
        IOobject::groupName(physicalProperties::typeName, phaseName);

    if (mesh_.mesh().foundObject<basicThermo>(thermocName))
    {
        return mesh_.mesh().lookupObject<basicThermo>(thermocName).rho();
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

    if (mesh_.mesh().foundObject<volScalarField>(mucName))
    {
        return mesh_.mesh().lookupObject<volScalarField>(mucName);
    }

    const word thermocName =
        IOobject::groupName(physicalProperties::typeName, phaseName);

    if (mesh_.mesh().foundObject<fluidThermo>(thermocName))
    {
        return mesh_.mesh().lookupObject<fluidThermo>(thermocName).mu();
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
        trhocVf_.ref() = getRhocVf(carrierPhaseName());
    }

    if (trhocPhaseVf_.isTmp())
    {
        trhocPhaseVf_.ref() = getRhocVf(phaseName());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::coupledToFluid::coupledToFluid
(
    const cloud& c,
    const dictionary& dict
)

:
    coupled(c, dict),
    mesh_(c.mesh()),
    trhocVf_(getRhocVf(carrierPhaseName())),
    trhocPhaseVf_
    (
        hasPhase()
      ? getRhocVf(phaseName())
      : tmp<volScalarField>(NullObjectRef<volScalarField>())
    ),
    mucVf_(getMucVf(carrierPhaseName())),
    rhoc(carrierField<scalar>(trhocVf_())),
    rhocPhase
    (
        hasPhase()
      ? carrierField<scalar>(trhocPhaseVf_())
      : carrierField<scalar>
        (
            "rhocPhase",
            [&]()
            {
                FatalErrorInFunction
                    << "Cloud " << c.name() << " does not have a corresponding "
                    << "Eulerian phase density" << exit(FatalError);
                return tmp<volScalarField>(nullptr);
            }
        )
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
      : carrierField<scalar>(mucVf_)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupledToFluid::~coupledToFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
