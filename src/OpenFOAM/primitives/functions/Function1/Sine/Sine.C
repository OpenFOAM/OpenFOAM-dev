/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Sine<Type>::read(const dictionary& coeffs)
{
    t0_ = coeffs.lookupOrDefault<scalar>("t0", 0);
    amplitude_ = Function1<scalar>::New("amplitude", coeffs);
    frequency_ = Function1<scalar>::New("frequency", coeffs);
    scale_ = Function1<Type>::New("scale", coeffs);
    level_ = Function1<Type>::New("level", coeffs);
}


template<class Type>
Foam::Function1Types::Sine<Type>::Sine
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName)
{
    read(dict);
}


template<class Type>
Foam::Function1Types::Sine<Type>::Sine(const Sine<Type>& se)
:
    Function1<Type>(se),
    t0_(se.t0_),
    amplitude_(se.amplitude_, false),
    frequency_(se.frequency_, false),
    scale_(se.scale_, false),
    level_(se.level_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Sine<Type>::~Sine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Sine<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os  << token::END_STATEMENT << nl;
    os  << indent << word(this->name() + "Coeffs") << nl;
    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    os.writeKeyword("t0") << t0_ << token::END_STATEMENT << nl;
    amplitude_->writeData(os);
    frequency_->writeData(os);
    scale_->writeData(os);
    level_->writeData(os);
    os  << decrIndent << indent << token::END_BLOCK << endl;
}


// ************************************************************************* //
