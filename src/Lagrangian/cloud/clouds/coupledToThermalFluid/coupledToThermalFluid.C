/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "coupledToThermalFluid.H"
#include "basicThermo.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(coupledToThermalFluid, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::coupledToThermalFluid::coupledToThermalFluid
(
    const cloud& c,
    const dictionary& dict
)
:
    coupledToFluid(c, dict),
    mesh_(c.mesh()),
    thermoc_
    (
        mesh_.lookupObject<fluidThermo>
        (
            IOobject::groupName
            (
                physicalProperties::typeName,
                carrierPhaseName()
            )
        )
    ),
    thermocPhase_
    (
        hasPhase()
      ? mesh_.lookupObject<basicThermo>
        (
            IOobject::groupName
            (
                physicalProperties::typeName,
                phaseName()
            )
        )
      : NullObjectRef<basicThermo>()
    ),
    pc(carrierField<scalar>(thermoc_.p())),
    Tc(carrierField<scalar>(thermoc_.T())),
    hec(carrierField<scalar>(thermoc_.he())),
    hecPhase
    (
        hasPhase()
      ? carrierField<scalar>(thermocPhase_.he())
      : carrierField<scalar>
        (
            IOobject::groupName("hec", phaseName()),
            [&]()
            {
                FatalErrorInFunction
                    << "Cloud " << c.name() << " does not have a corresponding "
                    << "Eulerian phase enthalpy/energy" << exit(FatalError);
                return tmp<volScalarField>(nullptr);
            }
        )
    ),
    Cpvc(carrierField<scalar>(thermoc_.Cpv())),
    Cpc
    (
        &thermoc_.Cp() == &thermoc_.Cpv()
      ? Cpvc
      : carrierField<scalar>(thermoc_.Cp())
    ),
    Cvc
    (
        &thermoc_.Cv() == &thermoc_.Cpv()
      ? Cpvc
      : carrierField<scalar>(thermoc_.Cv())
    ),
    kappac(carrierField<scalar>(thermoc_.kappa())),
    Prc
    (
        c.derivedField<scalar>
        (
            [&]
            (
                const LagrangianModelRef& model,
                const LagrangianSubMesh& subMesh
            )
            {
                return
                    Cpc(model, subMesh)
                   *muc(model, subMesh)
                   /kappac(model, subMesh);
            }
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupledToThermalFluid::~coupledToThermalFluid()
{}


// ************************************************************************* //
