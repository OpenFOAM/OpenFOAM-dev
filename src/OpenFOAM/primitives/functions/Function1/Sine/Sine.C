/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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

#include "Sine.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Sine<Type>::read(const dictionary& coeffs)
{
    amplitude_ = Function1<Type>::New("amplitude", coeffs);
    frequency_ = coeffs.lookup<scalar>("frequency");
    start_ = coeffs.lookupOrDefault<scalar>("start", 0);
    level_ = Function1<Type>::New("level", coeffs);

    integrable_ =
        isA<Constant<Type>>(amplitude_())
     && isA<Constant<Type>>(level_());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Sine<Type>::Sine
(
    const word& entryName,
    const dictionary& dict
)
:
    FieldFunction1<Type, Sine<Type>>(entryName)
{
    read(dict);
}


template<class Type>
Foam::Function1s::Sine<Type>::Sine(const Sine<Type>& se)
:
    FieldFunction1<Type, Sine<Type>>(se),
    amplitude_(se.amplitude_, false),
    frequency_(se.frequency_),
    start_(se.start_),
    level_(se.level_, false),
    integrable_(se.integrable_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Sine<Type>::~Sine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Sine<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os  << token::END_STATEMENT << nl;
    os  << indent << word(this->name() + "Coeffs") << nl;
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    amplitude_->writeData(os);
    writeEntry(os, "frequency", frequency_);
    writeEntry(os, "start", start_);
    level_->writeData(os);
    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// ************************************************************************* //
