/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "APIdiffCoef.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function2s
{
    addScalarFunction2(APIdiffCoef);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function2s::APIdiffCoef::APIdiffCoef
(
    const word& name,
    const scalar a,
    const scalar b,
    const scalar wf,
    const scalar wa
)
:
    FieldFunction2<scalar, APIdiffCoef>(name),
    a_(a),
    b_(b),
    wf_(wf),
    wa_(wa),
    alpha_(sqrt(1/wf_ + 1/wa_)),
    beta_(sqr(cbrt(a_) + cbrt(b_)))
{}


Foam::Function2s::APIdiffCoef::APIdiffCoef
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction2<scalar, APIdiffCoef>(name),
    a_(dict.lookup<scalar>("a")),
    b_(dict.lookup<scalar>("b")),
    wf_(dict.lookup<scalar>("wf")),
    wa_(dict.lookup<scalar>("wa")),
    alpha_(sqrt(1/wf_ + 1/wa_)),
    beta_(sqr((cbrt(a_) + cbrt(b_))))
{
    assertNoConvertUnits(typeName, units, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Function2s::APIdiffCoef::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, "a", a_);
    writeEntry(os, "b", b_);
    writeEntry(os, "wf", wf_);
    writeEntry(os, "wa", wa_);
}


// ************************************************************************* //
