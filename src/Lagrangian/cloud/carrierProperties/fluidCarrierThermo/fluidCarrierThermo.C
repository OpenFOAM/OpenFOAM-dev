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

#include "fluidCarrierThermo.H"
#include "fluidThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidCarrierThermo, 0);
    defineRunTimeSelectionTable(fluidCarrierThermo, cloud);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidCarrierThermo::implementation::implementation
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const basicThermo& thermo
)
:
    thermo_(refCast<const fluidThermo>(thermo)),
    p_(carriedCloud.carrierField<scalar>(thermo_.p())),
    mu_(carriedCloud.carrierField<scalar>(thermo_.mu())),
    nu_
    (
        c.derivedField<scalar>
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
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidCarrierThermo>
Foam::fluidCarrierThermo::New
(
    const cloud& c,
    const clouds::carried& carriedCloud,
    const word& phaseName
)
{
    return
        carrierThermo::New<fluidCarrierThermo, fluidThermo>
        (
            c,
            carriedCloud,
            phaseName
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidCarrierThermo::~fluidCarrierThermo()
{}


Foam::fluidCarrierThermo::implementation::~implementation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::CarrierField<Foam::scalar>&
Foam::fluidCarrierThermo::implementation::p() const
{
    return p_;
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::fluidCarrierThermo::implementation::p
(
    const LagrangianSubMesh& subMesh
) const
{
    return p_.ref(subMesh);
}


const Foam::CloudDerivedField<Foam::scalar>&
Foam::fluidCarrierThermo::implementation::mu() const
{
    return mu_;
}


Foam::tmp<Foam::LagrangianSubScalarSubField>
Foam::fluidCarrierThermo::implementation::mu
(
    const LagrangianSubMesh& subMesh
) const
{
    return mu_.ref(subMesh);
}


const Foam::CloudDerivedField<Foam::scalar>&
Foam::fluidCarrierThermo::implementation::nu() const
{
    return nu_;
}


Foam::tmp<Foam::LagrangianSubScalarField>
Foam::fluidCarrierThermo::implementation::nu
(
    const LagrangianSubMesh& subMesh
) const
{
    return nu_(subMesh);
}


// ************************************************************************* //
