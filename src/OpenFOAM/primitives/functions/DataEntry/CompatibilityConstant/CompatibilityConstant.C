/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "CompatibilityConstant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CompatibilityConstant<Type>::CompatibilityConstant
(
    const word& entryName,
    const dictionary& dict
)
:
    DataEntry<Type>(entryName),
    value_(pTraits<Type>::zero),
    dimensions_(dimless)
{
    Istream& is(dict.lookup(entryName));

    token firstToken(is);
    if (firstToken.isWord())
    {
        token nextToken(is);
        if (nextToken == token::BEGIN_SQR)
        {
            is.putBack(nextToken);
            is >> dimensions_;
            is >> value_;
        }
    }
    else
    {
        is.putBack(firstToken);
        is  >> value_;
    }
}


template<class Type>
Foam::CompatibilityConstant<Type>::CompatibilityConstant
(
    const CompatibilityConstant<Type>& cnst
)
:
    DataEntry<Type>(cnst),
    value_(cnst.value_),
    dimensions_(cnst.dimensions_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::CompatibilityConstant<Type>::~CompatibilityConstant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::CompatibilityConstant<Type>::value(const scalar x) const
{
    return value_;
}


template<class Type>
Type Foam::CompatibilityConstant<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    return (x2 - x1)*value_;
}


template<class Type>
Foam::dimensioned<Type> Foam::CompatibilityConstant<Type>::
dimValue(const scalar x) const
{
    return dimensioned<Type>("dimensionedValue", dimensions_, value_);
}


template<class Type>
Foam::dimensioned<Type> Foam::CompatibilityConstant<Type>::dimIntegrate
(
    const scalar x1, const scalar x2
) const
{
    return dimensioned<Type>("dimensionedValue", dimensions_, (x2-x1)*value_);
}

// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "CompatibilityConstantIO.C"


// ************************************************************************* //
