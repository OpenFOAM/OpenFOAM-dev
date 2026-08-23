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

#include "LagrangianSubFieldsFwd.H"
#include "coupledToThermalFluid.H"
#include "multicomponentLagrangianThermo.H"
#include "multicomponentCarrierThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace clouds
{
    defineTypeNameAndDebug(coupledToThermalFluid, 0);
}
}


const Foam::UPtrList<const Foam::CarrierField<Foam::scalar>>
    Foam::clouds::coupledToThermalFluid::nullY_;


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::clouds::coupledToThermalFluid::updateCarrier()
{
    coupledToFluid::updateCarrier();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::clouds::coupledToThermalFluid::coupledToThermalFluid
(
    const cloud& c,
    const carried& carriedCloud,
    const thermal& thermalCloud
)
:
    autoPtr<fluidCarrierThermo>
    (
        fluidCarrierThermo::New
        (
            c,
            carriedCloud,
            carriedCloud.carrierPhaseName()
        ).ptr()
    ),
    autoPtr<fluidSurfaceThermo>
    (
        autoPtr<fluidCarrierThermo>::valid()
      ? fluidSurfaceThermo::New(autoPtr<fluidCarrierThermo>::operator()()).ptr()
      : nullptr
    ),
    autoPtr<carrierThermo>
    (
        carriedCloud.hasPhase()
      ? carrierThermo::New
        (
            c,
            carriedCloud,
            carriedCloud.phaseName()
        ).ptr()
      : nullptr
    ),
    coupledToFluid
    (
        c,
        carriedCloud,
        autoPtr<fluidCarrierThermo>::valid()
      ? autoPtr<fluidCarrierThermo>::operator()()
      : NullObjectRef<carrierDynamicViscosity>(),
        autoPtr<fluidSurfaceThermo>::valid()
      ? autoPtr<fluidSurfaceThermo>::operator()()
      : NullObjectRef<fluidSurfaceThermo>(),
        autoPtr<carrierThermo>::valid()
      ? autoPtr<carrierThermo>::operator()()
      : NullObjectRef<carrierDensity>()
    ),
    cloud_(c),
    pc
    (
        hasThermoc()
      ? thermoc().p()
      : carriedCloud.carrierField<scalar>
        (
            cloud_.mesh().poly().lookupObject<volScalarField>("p")
        )
    ),
    Tc
    (
        hasThermoc()
      ? thermoc().T()
      : carriedCloud.carrierField<scalar>
        (
            cloud_.mesh().poly().lookupObject<volScalarField>("T")
        )
    ),
    hec
    (
        hasThermoc()
      ? thermoc().he()
      : carriedCloud.noCarrierField<scalar>("he", "enthalpy/energy", false)
    ),
    hecPhase
    (
        hasThermocPhase()
      ? thermocPhase().he()
      : carriedCloud.noCarrierField<scalar>("he", "enthalpy/energy", true)
    ),
    Yc
    (
        isThermoc<multicomponentCarrierThermo>()
      ? thermoc<multicomponentCarrierThermo>().Y()
      : nullY_
    ),
    YcPhase
    (
        isThermocPhase<multicomponentCarrierThermo>()
      ? thermocPhase<multicomponentCarrierThermo>().Y()
      : nullY_
    ),
    iToic(),
    iToicPhase()
{
    // Create maps from cloud to Eulerian specie indices
    if (thermalCloud.isThermo<multicomponentLagrangianThermo>())
    {
        const multicomponentLagrangianThermo& thermo =
            thermalCloud.thermo<multicomponentLagrangianThermo>();

        iToic.resize(thermo.species().size(), -1);
        iToicPhase.resize(thermo.species().size(), -1);

        if (isThermoc<multicomponentCarrierThermo>())
        {
            const multicomponentCarrierThermo& thermoc =
                this->thermoc<multicomponentCarrierThermo>();

            forAll(iToic, i)
            {
                const word& specieName = thermo.species()[i];

                iToic[i] =
                    thermoc.containsSpecie(specieName)
                  ? thermoc.species()[specieName]
                  : -1;
            }
        }

        if (isThermocPhase<multicomponentCarrierThermo>())
        {
            const multicomponentCarrierThermo& thermocPhase =
                this->thermocPhase<multicomponentCarrierThermo>();

            forAll(iToicPhase, i)
            {
                const word& specieName = thermo.species()[i];

                iToicPhase[i] =
                    thermocPhase.containsSpecie(specieName)
                  ? thermocPhase.species()[specieName]
                  : -1;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::clouds::coupledToThermalFluid::~coupledToThermalFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::clouds::coupledToThermalFluid::hasThermoc() const
{
    return autoPtr<fluidCarrierThermo>::valid();
}


bool Foam::clouds::coupledToThermalFluid::hasThermocs() const
{
    return autoPtr<fluidSurfaceThermo>::valid();
}


bool Foam::clouds::coupledToThermalFluid::hasThermocPhase() const
{
    return autoPtr<carrierThermo>::valid();
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::clouds::coupledToThermalFluid::kappac
(
    const LagrangianSubMesh& subMesh
) const
{
    return thermocs().kappaNear(subMesh);
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::clouds::coupledToThermalFluid::Prc
(
    const LagrangianSubMesh& subMesh
) const
{
    return thermocs().PrNear(subMesh);
}


// ************************************************************************* //
