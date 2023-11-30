/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "solidProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidProperties, 0);
    defineRunTimeSelectionTable(solidProperties,);
    defineRunTimeSelectionTable(solidProperties, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidProperties::solidProperties
(
    scalar rho,
    scalar Cp,
    scalar kappa,
    scalar hf,
    scalar emissivity
)
:
    rho_(rho),
    Cp_(Cp),
    kappa_(kappa),
    hf_(hf),
    emissivity_(emissivity)
{}


Foam::solidProperties::solidProperties(const dictionary& dict)
:
    rho_(dict.lookup<scalar>("rho")),
    Cp_(dict.lookup<scalar>("Cp")),
    kappa_(dict.lookupBackwardsCompatible<scalar>({"kappa", "K"})),
    hf_(dict.lookupBackwardsCompatible<scalar>({"hf", "Hf"})),
    emissivity_(dict.lookup<scalar>("emissivity"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidProperties::readIfPresent(const dictionary& dict)
{
    dict.readIfPresent("rho", rho_);
    dict.readIfPresent("Cp", Cp_);
    dict.readIfPresent("K", kappa_);
    dict.readIfPresent("kappa", kappa_);
    dict.readIfPresent("hf", hf_);
    dict.readIfPresent("emissivity", emissivity_);
}


void Foam::solidProperties::write(Ostream& os) const
{
    os  << rho_ << token::SPACE
        << Cp_ << token::SPACE
        << kappa_ << token::SPACE
        << hf_ << token::SPACE
        << emissivity_;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const solidProperties& s)
{
    s.write(os);
    return os;
}


// ************************************************************************* //
