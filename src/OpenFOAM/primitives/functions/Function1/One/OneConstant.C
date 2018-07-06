/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2018 OpenFOAM Foundation
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

#include "OneConstant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::OneConstant<Type>::OneConstant(const word& entryName)
:
    Function1<Type>(entryName)
{}


template<class Type>
Foam::Function1Types::OneConstant<Type>::OneConstant
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::OneConstant<Type>::~OneConstant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1Types::OneConstant<Type>::value
(
    const scalarField& x
) const
{
    return tmp<Field<Type>>(new Field<Type>(x.size(), pTraits<Type>::one));
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::Function1Types::OneConstant<Type>::integrate
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    return (x2 - x1)*pTraits<Type>::one;
}


template<class Type>
void Foam::Function1Types::OneConstant<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << token::END_STATEMENT << nl;
}


// ************************************************************************* //
