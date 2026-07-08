/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2026 OpenFOAM Foundation
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

#include "ArdenBuck.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{
    addStreamConstructableScalarFunction1(ArdenBuck);
}
}


const Foam::scalar Foam::Function1s::ArdenBuck::zeroC_ =
    Foam::units::lookup("K").toStandard(273.15);

const Foam::scalar Foam::Function1s::ArdenBuck::A_ =
    Foam::units::lookup("Pa").toStandard(611.21);

const Foam::scalar Foam::Function1s::ArdenBuck::B_ = 18.678;

const Foam::scalar Foam::Function1s::ArdenBuck::C_ =
    Foam::units::lookup("K").toStandard(234.5);

const Foam::scalar Foam::Function1s::ArdenBuck::D_ =
    Foam::units::lookup("K").toStandard(257.14);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::ArdenBuck::ArdenBuck
(
    const word& name,
    const unitSets& units,
    const dictionary& dict
)
:
    FieldFunction1<scalar, ArdenBuck>(name)
{
    assertNoConvertUnits(typeName, units, dict);
}


Foam::Function1s::ArdenBuck::ArdenBuck
(
    const word& name,
    const unitSets& units,
    Istream& is
)
:
    FieldFunction1<scalar, ArdenBuck>(name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1s::ArdenBuck::~ArdenBuck()
{}


void Foam::Function1s::ArdenBuck::write
(
    Ostream& os,
    const unitSets& units
) const
{}


// ************************************************************************* //
