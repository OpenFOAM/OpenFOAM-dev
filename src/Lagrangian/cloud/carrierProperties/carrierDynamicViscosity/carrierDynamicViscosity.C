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

#include "carrierDynamicViscosity.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(carrierDynamicViscosity, 0);
}


// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

const Foam::volScalarField& Foam::carrierDynamicViscosity::getMuVf
(
    const cloud& c,
    const word& phaseName,
    const bool allowNull
)
{
    const word muName = IOobject::groupName("mu", phaseName);

    if (c.mesh().poly().foundObject<volScalarField>(muName))
    {
        return c.mesh().poly().lookupObject<volScalarField>(muName);
    }

    const word thermoName =
        IOobject::groupName(physicalProperties::typeName, phaseName);

    if (c.mesh().poly().foundObject<fluidThermo>(thermoName))
    {
        return c.mesh().poly().lookupObject<fluidThermo>(thermoName).mu();
    }

    if (!allowNull)
    {
        FatalErrorInFunction
            << "Could not determine the"
            << (phaseName == word::null ? " carrier" : "")
            << " dynamic viscosity"
            << (phaseName == word::null ? "" : " of phase " + phaseName)
            << " for cloud " << c.name()
            << exit(FatalError);
    }

    return NullObjectRef<volScalarField>();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::carrierDynamicViscosity::implementation::implementation
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const word& phaseName
)
:
    carrierDensity::implementation(c, carriedCloud, phaseName),
    cloud_(c),
    phaseName_(phaseName),
    tnuVf_(getNuVf(cloud_, phaseName, true)),
    muVf_(getMuVf(cloud_, phaseName, true)),
    nu_
    (
        tnuVf_.valid()
      ? carriedCloud.carrierField<scalar>(tnuVf_())
      : c.derivedField<scalar>
        (
            [&]
            (
                const LagrangianModelRef& model,
                const LagrangianSubMesh& subMesh
            )
            {
                return mu_(model, subMesh)/rho()(model, subMesh);
            }
        )
    ),
    mu_
    (
        notNull(muVf_)
      ? carriedCloud.carrierField<scalar>(muVf_)
      : c.derivedField<scalar>
        (
            [&]
            (
                const LagrangianModelRef& model,
                const LagrangianSubMesh& subMesh
            )
            {
                return nu_(model, subMesh)*rho()(model, subMesh);
            }
        )
    )
{
    if (!tnuVf_.valid() && isNull(muVf_))
    {
        FatalErrorInFunction
            << "Could not determine the"
            << (phaseName == word::null ? " carrier" : "")
            << " kinematic or dynamic viscosity"
            << (phaseName == word::null ? "" : " of phase " + phaseName)
            << " for cloud " << c.name()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::carrierDynamicViscosity>
Foam::carrierDynamicViscosity::New
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const word& phaseName
)
{
    return
        autoPtr<carrierDynamicViscosity>
        (
            new carrierDynamicViscosity::implementation
            (
                c,
                carriedCloud,
                phaseName
            )
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::carrierDynamicViscosity::~carrierDynamicViscosity()
{}


Foam::carrierDynamicViscosity::implementation::~implementation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::carrierDynamicViscosity::implementation::correct()
{
    if (tnuVf_.isTmp())
    {
        tnuVf_.ref() == getNuVf(cloud_, phaseName_, false);
    }
}


const Foam::CloudDerivedField<Foam::scalar>&
Foam::carrierDynamicViscosity::implementation::nu() const
{
    return nu_;
}


const Foam::CloudDerivedField<Foam::scalar>&
Foam::carrierDynamicViscosity::implementation::mu() const
{
    return mu_;
}


// ************************************************************************* //
