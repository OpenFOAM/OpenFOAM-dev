/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "Constant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Constant<Type>::Constant
(
    const word& entryName,
    const Type& val
)
:
    FieldFunction1<Type, Constant<Type>>(entryName),
    value_(val)
{}


template<class Type>
Foam::Function1s::Constant<Type>::Constant
(
    const word& entryName,
    const dictionary& dict
)
:
    FieldFunction1<Type, Constant<Type>>(entryName),
    value_(Zero)
{
    if (!dict.found(entryName))
    {
        dict.lookup("value") >> value_;
    }
    else
    {
        Istream& is(dict.lookup(entryName));
        word entryType(is);
        if (is.eof())
        {
            dict.lookup("value") >> value_;
        }
        else
        {
            is  >> value_;
        }
    }
}


template<class Type>
Foam::Function1s::Constant<Type>::Constant
(
    const word& entryName,
    Istream& is
)
:
    FieldFunction1<Type, Constant<Type>>(entryName),
    value_(pTraits<Type>(is))
{}


template<class Type>
Foam::Function1s::Constant<Type>::Constant(const Constant<Type>& cnst)
:
    FieldFunction1<Type, Constant<Type>>(cnst),
    value_(cnst.value_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Constant<Type>::~Constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1s::Constant<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << token::SPACE << value_ << token::END_STATEMENT << nl;
}


// ************************************************************************* i/
