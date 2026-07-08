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

#include "AntoineExtended.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{
    addScalarFunction1(AntoineExtended);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::AntoineExtended::AntoineExtended
(
    const word& name,
    const unitSets& units,
    const dictionary& dict
)
:
    FieldFunction1<scalar, AntoineExtended>(name),
    A_(dict.lookup<scalar>("A", dimless)),
    B_(dict.lookup<scalar>("B", units.x)),
    C_(dict.lookup<scalar>("C", units.x)),
    D_(dict.lookup<scalar>("D", dimless/units.x)),
    E_(dict.lookup<scalar>("E", dimless)),
    G_(dict.lookup<scalar>("G", dimless)),
    F_(dict.lookup<scalar>("F", dimless/pow(units.x, G_)))
{
    assertNoConvertUnits(typeName, {units::any, units.value}, dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1s::AntoineExtended::~AntoineExtended()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Function1s::AntoineExtended::write
(
    Ostream& os,
    const unitSets& units
) const
{
    writeEntry(os, "A", units::unitless, A_);
    writeEntry(os, "B", units.x, B_);
    writeEntry(os, "C", units.x, C_);
    writeEntry(os, "D", dimless/units.x, D_);
    writeEntry(os, "E", dimless, E_);
    writeEntry(os, "G", dimless, G_);
    writeEntry(os, "F", dimless/pow(units.x, G_), F_);
}


// ************************************************************************* //
