/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

#include "phasePairKey.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phasePairKey::hash::hash()
{}


Foam::phasePairKey::phasePairKey()
{}


Foam::phasePairKey::phasePairKey
(
    const word& name1,
    const word& name2,
    const bool ordered
)
:
    Pair<word>(name1, name2),
    ordered_(ordered)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phasePairKey::~phasePairKey()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::phasePairKey::ordered() const
{
    return ordered_;
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

Foam::label Foam::phasePairKey::hash::operator()
(
    const phasePairKey& key
) const
{
    if (key.ordered_)
    {
        return
            word::hash()
            (
                key.first(),
                word::hash()(key.second())
            );
    }
    else
    {
        return
            word::hash()(key.first())
          + word::hash()(key.second());
    }
}


// * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * * //

bool Foam::operator==
(
    const phasePairKey& a,
    const phasePairKey& b
)
{
    const label c = Pair<word>::compare(a, b);

    return
        (a.ordered_ == b.ordered_)
     && (
            (a.ordered_ && (c == 1))
         || (!a.ordered_ && (c != 0))
        );
}


bool Foam::operator!=
(
    const phasePairKey& a,
    const phasePairKey& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, phasePairKey& key)
{
    const FixedList<word, 3> temp(is);

    key.first() = temp[0];

    if (temp[1] == "and")
    {
        key.ordered_ = false;
    }
    else if (temp[1] == "in")
    {
        key.ordered_ = true;
    }
    else
    {
        FatalErrorInFunction
            << "Phase pair type is not recognised. "
            << temp
            << "Use (phaseDispersed in phaseContinuous) for an ordered"
            << "pair, or (phase1 and pase2) for an unordered pair."
            << exit(FatalError);
    }

    key.second() = temp[2];

    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const phasePairKey& key)
{
    os  << token::BEGIN_LIST
        << key.first()
        << token::SPACE
        << (key.ordered_ ? "in" : "and")
        << token::SPACE
        << key.second()
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
