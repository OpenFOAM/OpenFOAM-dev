/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "Polynomial1.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::List<Type> Foam::Function1s::Polynomial<Type>::convertRead
(
    const unitConversions& units,
    List<Type> coeffs
)
{
    unitConversion unitsXPowI(units.x);
    forAll(coeffs, i)
    {
        coeffs[i] = units.value.toStandard(unitsXPowI.toUser(coeffs[i]));
        unitsXPowI.reset(units.x*unitsXPowI);
    }
    return coeffs;
}


template<class Type>
Foam::List<Type> Foam::Function1s::Polynomial<Type>::convertWrite
(
    const unitConversions& units,
    List<Type> coeffs
)
{
    unitConversion unitsXPowI(units.x);
    forAll(coeffs, i)
    {
        coeffs[i] = units.value.toUser(unitsXPowI.toStandard(coeffs[i]));
        unitsXPowI.reset(units.x*unitsXPowI);
    }
    return coeffs;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Polynomial<Type>::Polynomial
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction1<Type, Polynomial<Type>>(name),
    coeffs_
    (
        convertRead
        (
            dict.found("units")
          ? Function1s::unitConversions(dict.lookup("units"))
          : units,
            dict.lookup("coeffs")
        )
    )
{
    if (!coeffs_.size())
    {
        FatalIOErrorInFunction(dict)
            << typeName.capitalise() << ' ' << name
            << " must have at least one coefficient"
            << exit(FatalError);
    }
}


template<class Type>
Foam::Function1s::Polynomial<Type>::Polynomial
(
    const word& name,
    const unitConversions& units,
    Istream& is
)
:
    FieldFunction1<Type, Polynomial<Type>>(name),
    coeffs_(convertRead(units, is))
{
    if (!coeffs_.size())
    {
        FatalIOErrorInFunction(is)
            << typeName.capitalise() << ' ' << name
            << " must have at least one coefficient"
            << exit(FatalError);
    }
}


template<class Type>
Foam::Function1s::Polynomial<Type>::Polynomial(const Polynomial& poly)
:
    FieldFunction1<Type, Polynomial<Type>>(poly),
    coeffs_(poly.coeffs_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::Polynomial<Type>::~Polynomial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1s::Polynomial<Type>::value(const scalar x) const
{
    Type y = coeffs_[coeffs_.size() - 1];

    for (label i = coeffs_.size() - 2; i >= 0; i --)
    {
        y = y*x + coeffs_[i];
    }

    return y;
}


template<class Type>
Type Foam::Function1s::Polynomial<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    Type Sy1 = coeffs_[coeffs_.size() - 1]/coeffs_.size(), Sy2 = Sy1;

    for (label i = coeffs_.size() - 2; i >= 0; i --)
    {
        Sy1 = Sy1*x1 + coeffs_[i]/(i + 1);
        Sy2 = Sy2*x2 + coeffs_[i]/(i + 1);
    }

    return x2*Sy2 - x1*Sy1;
}


template<class Type>
void Foam::Function1s::Polynomial<Type>::write
(
    Ostream& os,
    const unitConversions& units
) const
{
    writeEntry(os, "coeffs", convertWrite(units, coeffs_));
}


// ************************************************************************* //
