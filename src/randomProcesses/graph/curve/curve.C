/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "curve.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::curve::curve
(
    const string& name,
    const curveStyle& style,
    const label l
)
:
    scalarField(l, 0.0),
    name_(name),
    style_(style)
{}


Foam::curve::curve
(
    const string& name,
    const curveStyle& style,
    const scalarField& y
)
:
    scalarField(y),
    name_(name),
    style_(style)
{}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const curve& c)
{
    os  << nl
        << c.name_ << nl
        << c.style_ << nl
        << static_cast<const scalarField&>(c);

    os.check("Ostream& operator>>(Ostream&, const curve&)");

    return os;
}


// ************************************************************************* //
