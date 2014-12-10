/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "TimeDataEntry.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::TimeDataEntry<Type>::TimeDataEntry
(
    const Time& t,
    const word& name,
    const dictionary& dict
)
:
    time_(t),
    name_(name),
    entry_(DataEntry<Type>::New(name, dict))
{
    entry_->convertTimeBase(t);
}


template<class Type>
Foam::TimeDataEntry<Type>::TimeDataEntry(const Time& t, const word& name)
:
    time_(t),
    name_(name),
    entry_(NULL)
{}


template<class Type>
Foam::TimeDataEntry<Type>::TimeDataEntry
(
    const TimeDataEntry<Type>& tde
)
:
    time_(tde.time_),
    name_(tde.name_),
    entry_()
{
    if (tde.entry_.valid())
    {
        entry_.reset(tde.entry_->clone().ptr());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::TimeDataEntry<Type>::~TimeDataEntry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::TimeDataEntry<Type>::reset(const dictionary& dict)
{
    entry_.reset
    (
        DataEntry<Type>::New
        (
            name_,
            dict
        ).ptr()
    );

    entry_->convertTimeBase(time_);
}


template<class Type>
const Foam::word& Foam::TimeDataEntry<Type>::name() const
{
    return entry_->name();
}


template<class Type>
Type Foam::TimeDataEntry<Type>::value(const scalar x) const
{
    return entry_->value(x);
}


template<class Type>
Type Foam::TimeDataEntry<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    return entry_->integrate(x1, x2);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const TimeDataEntry<Type>& de
)
{
    return de.entry_->operator<<(os, de);
}


template<class Type>
void Foam::TimeDataEntry<Type>::writeData(Ostream& os) const
{
    entry_->writeData(os);
}


// ************************************************************************* //
