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

#include "carrierDensity.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(carrierDensity, 0);
}


// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::carrierDensity::getRhoVf
(
    const cloud& c,
    const word& phaseName,
    const bool allowEmpty
)
{
    const word rhoName = IOobject::groupName("rho", phaseName);

    if (c.mesh().poly().foundObject<volScalarField>(rhoName))
    {
        return c.mesh().poly().lookupObject<volScalarField>(rhoName);
    }

    const word thermoName =
        IOobject::groupName(physicalProperties::typeName, phaseName);

    if (c.mesh().poly().foundObject<basicThermo>(thermoName))
    {
        return c.mesh().poly().lookupObject<basicThermo>(thermoName).rho();
    }

    if (!allowEmpty)
    {
        FatalErrorInFunction
            << "Could not determine the"
            << (phaseName == word::null ? " carrier" : "")
            << " density"
            << (phaseName == word::null ? "" : " of phase " + phaseName)
            << " for cloud " << c.name()
            << exit(FatalError);
    }

    return tmp<volScalarField>(nullptr);

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::carrierDensity::implementation::implementation
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const word& phaseName
)
:
    cloud_(c),
    phaseName_(phaseName),
    trhoVf_(getRhoVf(cloud_, phaseName_, false)),
    rho_(carriedCloud.carrierField<scalar>(trhoVf_()))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::carrierDensity> Foam::carrierDensity::New
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const word& phaseName
)
{
    return
        autoPtr<carrierDensity>
        (
            new carrierDensity::implementation(c, carriedCloud, phaseName)
        );
}



Foam::carrierDensity::~carrierDensity()
{}


Foam::carrierDensity::implementation::~implementation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::carrierDensity::implementation::correct()
{
    if (trhoVf_.isTmp())
    {
        trhoVf_.ref() == getRhoVf(cloud_, phaseName_, false);
    }
}


const Foam::CarrierField<Foam::scalar>&
Foam::carrierDensity::implementation::rho() const
{
    return rho_;
}


// ************************************************************************* //
