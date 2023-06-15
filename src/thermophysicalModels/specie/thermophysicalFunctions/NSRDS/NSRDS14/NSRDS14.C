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

#include "NSRDS14.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{
    makeScalarFunction1(NSRDS14);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::NSRDS14::NSRDS14
(
    const word& name,
    const scalar Tc,
    const scalar a,
    const scalar b,
    const scalar c,
    const scalar d
)
:
    FieldFunction1<scalar, NSRDS14>(name),
    Tc_(Tc),
    a_(a),
    b_(b),
    c_(c),
    d_(d)
{}


Foam::Function1s::NSRDS14::NSRDS14
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction1<scalar, NSRDS14>(name),
    Tc_(dict.lookup<scalar>("Tc")),
    a_(dict.lookup<scalar>("a")),
    b_(dict.lookup<scalar>("b")),
    c_(dict.lookup<scalar>("c")),
    d_(dict.lookup<scalar>("d"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1s::NSRDS14::integral
(
    const scalar x1,
    const scalar x2
) const
{
    NotImplemented;
    return 0;
}


void Foam::Function1s::NSRDS14::write(Ostream& os) const
{
    writeEntry(os, "Tc", Tc_);
    writeEntry(os, "a", a_);
    writeEntry(os, "b", b_);
    writeEntry(os, "c", c_);
    writeEntry(os, "d", d_);
}


// ************************************************************************* //
