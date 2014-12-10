/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    defineRunTimeSelectionTable(solidProperties, Istream);
    defineRunTimeSelectionTable(solidProperties, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidProperties::solidProperties
(
    scalar rho,
    scalar Cp,
    scalar K,
    scalar Hf,
    scalar emissivity
)
:
    rho_(rho),
    Cp_(Cp),
    K_(K),
    Hf_(Hf),
    emissivity_(emissivity)
{}


Foam::solidProperties::solidProperties(Istream& is)
:
    rho_(readScalar(is)),
    Cp_(readScalar(is)),
    K_(readScalar(is)),
    Hf_(readScalar(is)),
    emissivity_(readScalar(is))
{}


Foam::solidProperties::solidProperties(const dictionary& dict)
:
    rho_(readScalar(dict.lookup("rho"))),
    Cp_(readScalar(dict.lookup("Cp"))),
    K_(readScalar(dict.lookup("K"))),
    Hf_(readScalar(dict.lookup("Hf"))),
    emissivity_(readScalar(dict.lookup("emissivity")))
{}


Foam::solidProperties::solidProperties(const solidProperties& s)
:
    rho_(s.rho_),
    Cp_(s.Cp_),
    K_(s.K_),
    Hf_(s.Hf_),
    emissivity_(s.emissivity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidProperties::writeData(Ostream& os) const
{
    os  << rho_ << token::SPACE
        << Cp_ << token::SPACE
        << K_ << token::SPACE
        << Hf_ << token::SPACE
        << emissivity_;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const solidProperties& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
