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

#include "carrierViscosity.H"
#include "physicalProperties.H"
#include "viscosity.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(carrierViscosity, 0);
}


// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::carrierViscosity::getNuVf
(
    const cloud& c,
    const word& phaseName,
    const bool allowEmpty
)
{
    const word nuName = IOobject::groupName("nu", phaseName);

    if (c.mesh().poly().foundObject<volScalarField>(nuName))
    {
        return c.mesh().poly().lookupObject<volScalarField>(nuName);
    }

    const word viscosityName =
        IOobject::groupName(physicalProperties::typeName, phaseName);

    if (c.mesh().poly().foundObject<viscosity>(viscosityName))
    {
        return c.mesh().poly().lookupObject<viscosity>(viscosityName).nu();
    }

    if (!allowEmpty)
    {
        FatalErrorInFunction
            << "Could not determine the"
            << (phaseName == word::null ? " carrier" : "")
            << " kinematic viscosity"
            << (phaseName == word::null ? "" : " of phase " + phaseName)
            << " for cloud " << c.name()
            << exit(FatalError);
    }

    return tmp<volScalarField>(nullptr);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::carrierViscosity::implementation::implementation
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const word& phaseName
)
:
    cloud_(c),
    phaseName_(phaseName),
    tnuVf_(getNuVf(cloud_, phaseName_, false)),
    nu_(carriedCloud.carrierField<scalar>(tnuVf_()))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::carrierViscosity> Foam::carrierViscosity::New
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const word& phaseName
)
{
    return
        autoPtr<carrierViscosity>
        (
            new carrierViscosity::implementation(c, carriedCloud, phaseName)
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::carrierViscosity::~carrierViscosity()
{}


Foam::carrierViscosity::implementation::~implementation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::carrierViscosity::implementation::correct()
{
    if (tnuVf_.isTmp())
    {
        tnuVf_.ref() == getNuVf(cloud_, phaseName_, false);
    }
}


const Foam::CloudDerivedField<Foam::scalar>&
Foam::carrierViscosity::implementation::nu() const
{
    return nu_;
}


// ************************************************************************* //
