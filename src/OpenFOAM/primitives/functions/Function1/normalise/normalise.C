/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "normalise.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1s
{
    addScalarFunction1(normalise);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1s::normalise::normalise
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction1<scalar, normalise>(name),
    bounds_(dict.lookup<Pair<scalar>>("bounds", units.x)),
    value_(Function1<scalar>::New("value", {units.x, unitAny}, dict)),
    scale_(1/value_->integral(bounds_[0], bounds_[1]))
{}


Foam::Function1s::normalise::normalise(const normalise& se)
:
    FieldFunction1<scalar, normalise>(se),
    bounds_(se.bounds_),
    value_(se.value_, false),
    scale_(se.scale_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1s::normalise::~normalise()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Function1s::normalise::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, units.x, bounds_);
    writeEntry(os, units, value_());
}


// ************************************************************************* //
