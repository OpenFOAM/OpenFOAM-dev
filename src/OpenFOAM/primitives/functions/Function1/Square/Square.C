/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "Square.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::Function1s::Square<Type>::readDutyCycle
(
    const dictionary& dict
) const
{
    const bool haveDutyCycle = dict.found("dutyCycle");
    const bool haveMarkSpace = dict.found("markSpace");

    if (haveDutyCycle && haveMarkSpace)
    {
        FatalIOErrorInFunction(dict)
            << "both keywords dutyCycle and markSpace defined in dictionary "
            << dict.name() << exit(FatalIOError);
    }

    if (haveDutyCycle)
    {
        return dict.lookup<scalar>("dutyCycle", units::unitless);
    }

    if (haveMarkSpace)
    {
        const scalar r = dict.lookup<scalar>("markSpace", units::unitless);
        return r/(1 + r);
    }

    return 0.5;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Square<Type>::Square
(
    const word& name,
    const unitSets& units,
    const dictionary& dict
)
:
    FieldFunction1<Type, Square<Type>>(name),
    amplitude_(Function1<Type>::New("amplitude", units, dict)),
    frequency_(dict.lookup<scalar>("frequency", units::unitless/units.x)),
    start_(dict.lookupOrDefault<scalar>("start", units.x, 0)),
    level_(Function1<Type>::New("level", units, dict)),
    dutyCycle_(readDutyCycle(dict)),
    integrable_(amplitude_->constant() && level_->constant())
{}


template<class Type>
Foam::Function1s::Square<Type>::Square(const Square<Type>& se)
:
    FieldFunction1<Type, Square<Type>>(se),
    amplitude_(se.amplitude_, false),
    frequency_(se.frequency_),
    start_(se.start_),
    level_(se.level_, false),
    dutyCycle_(se.dutyCycle_),
    integrable_(se.integrable_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Square<Type>::~Square()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Square<Type>::write
(
    Ostream& os,
    const unitSets& units
) const
{
    writeEntry(os, units, amplitude_());
    writeEntry(os, "frequency", units::unitless/units.x, frequency_);
    writeEntry(os, "start", units.x, start_);
    writeEntry(os, units, level_());
    writeEntry(os, "dutyCycle", units::unitless, dutyCycle_);
}


// ************************************************************************* //
