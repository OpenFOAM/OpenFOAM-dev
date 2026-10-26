/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2026 OpenFOAM Foundation
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

#include "reaction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reaction, 0);
    defineRunTimeSelectionTable(reaction, dictionary);
    defineRunTimeSelectionTable(reaction, objectRegistry);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reaction::reaction(const reactionCoeffs& r)
:
    reactionCoeffs(r),
    lWork_(lhs().size()),
    rWork_(rhs().size())
{}


Foam::reaction::reaction
(
    const speciesTable& species,
    const dictionary& dict
)
:
    reactionCoeffs(species, dict),
    lWork_(lhs().size()),
    rWork_(rhs().size())
{}


Foam::reaction::state::state(const label nSpecie)
:
    c(nSpecie),
    specieWork_(scalarField(nSpecie))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::reaction::selectionType() const
{
    return type();
}


void Foam::reaction::write(Ostream& os) const
{
    reactionCoeffs::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Foam::Ostream& os, const reaction& r)
{
    r.write(os);
    return os;
}


// ************************************************************************* //
