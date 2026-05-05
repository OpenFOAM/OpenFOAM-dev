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
#include "multicomponentThermo.H"
#include "fluidThermo.H"
#include "fluidMulticomponentThermo.H"
#include "multicomponentLagrangianThermo.H"

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
    const carried& carriedCloud,
    const thermal& thermalCloud
)
:
    coupledToFluid(c, carriedCloud),
    mesh_(c.mesh()),
    thermoc_
    (
        mesh_.poly().lookupObject<fluidThermo>
        (
            IOobject::groupName
            (
                physicalProperties::typeName,
                carriedCloud.carrierPhaseName()
            )
        )
    ),
    multicomponentThermoc_
    (
        refCastNull<const fluidMulticomponentThermo>(thermoc_)
    ),
    thermocPhase_
    (
        carriedCloud.hasPhase()
      ? mesh_.poly().lookupObject<basicThermo>
        (
            IOobject::groupName
            (
                physicalProperties::typeName,
                carriedCloud.phaseName()
            )
        )
      : NullObjectRef<basicThermo>()
    ),
    multicomponentThermocPhase_
    (
        refCastNull<const multicomponentThermo>(thermocPhase_)
    ),
    pc(carriedCloud.carrierField<scalar>(thermoc_.p())),
    Tc(carriedCloud.carrierField<scalar>(thermoc_.T())),
    hec(carriedCloud.carrierField<scalar>(thermoc_.he())),
    hecPhase
    (
        carriedCloud.hasPhase()
      ? carriedCloud.carrierField<scalar>(thermocPhase_.he())
      : carriedCloud.carrierField<scalar>
        (
            IOobject::groupName("hec", carriedCloud.phaseName()),
            [&]()
            {
                FatalErrorInFunction
                    << "Cloud " << c.name() << " does not have a corresponding "
                    << "Eulerian phase enthalpy/energy" << exit(FatalError);
                return tmp<volScalarField>(nullptr);
            }
        )
    ),
    Cpvc(carriedCloud.carrierField<scalar>(thermoc_.Cpv())),
    Cpc
    (
        &thermoc_.Cp() == &thermoc_.Cpv()
      ? Cpvc
      : carriedCloud.carrierField<scalar>(thermoc_.Cp())
    ),
    Cvc
    (
        &thermoc_.Cv() == &thermoc_.Cpv()
      ? Cpvc
      : carriedCloud.carrierField<scalar>(thermoc_.Cv())
    ),
    kappac(carriedCloud.carrierField<scalar>(thermoc_.kappa())),
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
    ),
    iToic(),
    iToicPhase(),
    Yc(),
    YcPhase()
{
    // Create maps from cloud to Eulerian specie indices
    if (thermalCloud.isThermo<multicomponentLagrangianThermo>())
    {
        const multicomponentLagrangianThermo& thermo =
            thermalCloud.thermo<multicomponentLagrangianThermo>();

        iToic.resize(thermo.species().size(), -1);
        iToicPhase.resize(thermo.species().size(), -1);

        if (notNull(multicomponentThermoc_))
        {
            forAll(iToic, i)
            {
                const word& specieName = thermo.species()[i];

                iToic[i] =
                    multicomponentThermoc_.containsSpecie(specieName)
                  ? multicomponentThermoc_.species()[specieName]
                  : -1;
            }
        }
        if (notNull(multicomponentThermocPhase_))
        {
            forAll(iToicPhase, i)
            {
                const word& specieName = thermo.species()[i];

                iToicPhase[i] =
                    multicomponentThermocPhase_.containsSpecie(specieName)
                  ? multicomponentThermocPhase_.species()[specieName]
                  : -1;
            }
        }
    }

    // Create carrier fields for the mass fractions
    if (notNull(multicomponentThermoc_))
    {
        const PtrList<volScalarField>& YcVf =
            multicomponentThermoc_.Y();

        Yc.resize(YcVf.size());

        forAll(Yc, ic)
        {
            Yc.set
            (
                ic,
                &carriedCloud.carrierField<scalar>(YcVf[ic])
            );
        }
    }
    if (notNull(multicomponentThermocPhase_))
    {
        const PtrList<volScalarField>& YcPhaseVf =
            multicomponentThermocPhase_.Y();

        YcPhase.resize(YcPhaseVf.size());

        forAll(YcPhase, icPhase)
        {
            YcPhase.set
            (
                icPhase,
                &carriedCloud.carrierField<scalar>(YcPhaseVf[icPhase])
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupledToThermalFluid::~coupledToThermalFluid()
{}


// ************************************************************************* //
