/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "TimeFunction1.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type>
Foam::TimeFunction1<Type>::TimeFunction1
(
    const Time& t,
    const word& name,
    const dictionary& dict
)
:
    time_(t),
    name_(name),
    entry_(Function1<Type>::New(name, dict))
{
    entry_->convertTimeBase(t);
}


template<class Type>
Foam::TimeFunction1<Type>::TimeFunction1(const Time& t, const word& name)
:
    time_(t),
    name_(name),
    entry_(nullptr)
{}


template<class Type>
Foam::TimeFunction1<Type>::TimeFunction1
(
    const TimeFunction1<Type>& tde
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
Foam::TimeFunction1<Type>::~TimeFunction1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::TimeFunction1<Type>::reset(const dictionary& dict)
{
    entry_.reset
    (
        Function1<Type>::New
        (
            name_,
            dict
        ).ptr()
    );

    entry_->convertTimeBase(time_);
}


template<class Type>
const Foam::word& Foam::TimeFunction1<Type>::name() const
{
    return entry_->name();
}


template<class Type>
Type Foam::TimeFunction1<Type>::value(const scalar x) const
{
    return entry_->value(x);
}


template<class Type>
Type Foam::TimeFunction1<Type>::integrate
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
    const TimeFunction1<Type>& de
)
{
    return de.entry_->operator<<(os, de);
}


template<class Type>
void Foam::TimeFunction1<Type>::writeData(Ostream& os) const
{
    entry_->writeData(os);
}


// ************************************************************************* //
