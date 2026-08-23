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

#include "IOdictionary.H"
#include "LagrangianMesh.H"
#include "carrierThermo.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(carrierThermo, 0);
    defineRunTimeSelectionTable(carrierThermo, cloud);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::carrierThermo::implementation::implementation
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const basicThermo& thermo
)
:
    thermo_(thermo),
    mesh_(c.mesh()),
    T_(carriedCloud.carrierField<scalar>(thermo.T())),
    trhoVf_
    (
        thermo.mesh().foundObject<volScalarField>
        (
            thermo.phasePropertyName("rho")
        )
      ? tmp<volScalarField>
        (
            thermo.mesh().lookupObject<volScalarField>
            (
                thermo.phasePropertyName("rho")
            )
        )
      : thermo.rho()
    ),
    rho_(carriedCloud.carrierField<scalar>(trhoVf_())),
    he_(carriedCloud.carrierField<scalar>(thermo.he())),
    Cpv_(carriedCloud.carrierField<scalar>(thermo.Cpv())),
    Cp_
    (
        &thermo.Cp() == &thermo.Cpv()
      ? Cpv_
      : carriedCloud.carrierField<scalar>(thermo.Cp())
    ),
    Cv_
    (
        &thermo.Cv() == &thermo.Cpv()
      ? Cpv_
      : carriedCloud.carrierField<scalar>(thermo.Cv())
    ),
    kappa_(carriedCloud.carrierField<scalar>(thermo.kappa()))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::carrierThermo>
Foam::carrierThermo::New
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const word& phaseName
)
{
    return New<carrierThermo, basicThermo>(c, carriedCloud, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::carrierThermo::~carrierThermo()
{}


Foam::carrierThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::basicThermo& Foam::carrierThermo::implementation::thermo() const
{
    return thermo_;
}


const Foam::IOdictionary&
Foam::carrierThermo::implementation::properties() const
{
    return thermo_.properties();
}


const Foam::LagrangianMesh&
Foam::carrierThermo::implementation::mesh() const
{
    return mesh_;
}


const Foam::word& Foam::carrierThermo::implementation::phaseName() const
{
    return thermo_.phaseName();
}


void Foam::carrierThermo::implementation::correct()
{
    if (trhoVf_.isTmp())
    {
        trhoVf_.ref() == thermo_.rho();
    }
}


const Foam::CarrierField<Foam::scalar>&
Foam::carrierThermo::implementation::T() const
{
    return T_;
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::carrierThermo::implementation::T
(
    const LagrangianSubMesh& subMesh
) const
{
    return T_.ref(subMesh);
}


const Foam::CarrierField<Foam::scalar>&
Foam::carrierThermo::implementation::rho() const
{
    return rho_;
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::carrierThermo::implementation::rho
(
    const LagrangianSubMesh& subMesh
) const
{
    return rho_.ref(subMesh);
}


const Foam::CarrierField<Foam::scalar>&
Foam::carrierThermo::implementation::he() const
{
    return he_;
}


const Foam::CarrierField<Foam::scalar>&
Foam::carrierThermo::implementation::Cpv() const
{
    return Cpv_;
}


const Foam::CarrierField<Foam::scalar>&
Foam::carrierThermo::implementation::Cp() const
{
    return Cp_;
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::carrierThermo::implementation::Cp
(
    const LagrangianSubMesh& subMesh
) const
{
    return Cp_(subMesh);
}


const Foam::CarrierField<Foam::scalar>&
Foam::carrierThermo::implementation::Cv() const
{
    return Cv_;
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::carrierThermo::implementation::Cv
(
    const LagrangianSubMesh& subMesh
) const
{
    return Cv_.ref(subMesh);
}


const Foam::CarrierField<Foam::scalar>&
Foam::carrierThermo::implementation::kappa() const
{
    return kappa_;
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::carrierThermo::implementation::kappa
(
    const LagrangianSubMesh& subMesh
) const
{
    return kappa_.ref(subMesh);
}


// ************************************************************************* //
