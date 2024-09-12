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

#include "Repeat.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::Function1s::Repeat<Type>::readPeriod
(
    const unitConversions& units,
    const dictionary& dict
)
{
    const bool havePeriod = dict.found("period");
    const bool haveFrequency = dict.found("frequency");

    if (havePeriod == haveFrequency)
    {
        FatalIOErrorInFunction(dict)
            << (havePeriod ? "both" : "neither") << " of keywords "
            << "period " << (havePeriod ? "and" : "or")
            << " frequency defined in dictionary " << dict.name()
            << exit(FatalIOError);
    }

    return
        havePeriod
      ? dict.lookupOrDefault<scalar>("period", units.x, 0)
      : 1/dict.lookupOrDefault<scalar>("frequency", unitless/units.x, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Repeat<Type>::Repeat
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction1<Type, Repeat<Type>>(name),
    period_(readPeriod(units, dict)),
    start_(dict.lookupOrDefault<scalar>("start", units.x, scalar(0))),
    value_(Function1<Type>::New("value", units, dict))
{}


template<class Type>
Foam::Function1s::Repeat<Type>::Repeat(const Repeat<Type>& se)
:
    FieldFunction1<Type, Repeat<Type>>(se),
    period_(se.period_),
    start_(se.start_),
    value_(se.value_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Repeat<Type>::~Repeat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Repeat<Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, "period", units.x, period_);
    writeEntryIfDifferent(os, "start", units.x, scalar(0), start_);
    writeEntry(os, units, value_());
}


// ************************************************************************* //
