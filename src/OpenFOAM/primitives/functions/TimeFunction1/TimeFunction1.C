/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2020 OpenFOAM Foundation
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
    const Time& time,
    const word& name,
    const dictionary& dict
)
:
    time_(time),
    name_(name),
    function_(Function1<Type>::New(name, dict))
{}


template<class Type>
Foam::TimeFunction1<Type>::TimeFunction1
(
    const Time& time,
    const word& name
)
:
    time_(time),
    name_(name),
    function_(nullptr)
{}


template<class Type>
Foam::TimeFunction1<Type>::TimeFunction1
(
    const TimeFunction1<Type>& tf
)
:
    time_(tf.time_),
    name_(tf.name_),
    function_(tf.function_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::TimeFunction1<Type>::~TimeFunction1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::TimeFunction1<Type>::reset(const dictionary& dict)
{
    function_.reset(Function1<Type>::New(name_, dict).ptr());
}


template<class Type>
Type Foam::TimeFunction1<Type>::value(const scalar x) const
{
    return function_->value(time_.userTimeToTime(x));
}


template<class Type>
Type Foam::TimeFunction1<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    return
        time_.timeToUserTimeRatio()
       *function_->integral
        (
            time_.userTimeToTime(x1),
            time_.userTimeToTime(x2)
        );
}


template<class Type>
void Foam::TimeFunction1<Type>::write(Ostream& os) const
{
    writeEntry(os, function_());
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const TimeFunction1<Type>& tf
)
{
    return os << tf.function_();
}


// ************************************************************************* //
