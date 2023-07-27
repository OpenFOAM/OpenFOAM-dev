/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "DimensionedFunction1.H"
#include "dimensionedType.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::Function1s::Dimensioned<Type>::value(const scalar x) const
{
    return
        value_->value(x/xUnits_.value())
       *units_.value();
}


template<class Type>
inline Type Foam::Function1s::Dimensioned<Type>::integral
(
    const scalar x1,
    const scalar x2
) const
{
    return
        value_->integral(x1/xUnits_.value(), x2/xUnits_.value())
       *xUnits_.value()*units_.value();
}


template<class Type>
inline Foam::dimensioned<Type> Foam::Function1s::Dimensioned<Type>::value
(
    const dimensionedScalar& x
) const
{
    x.dimensions() = xUnits_.dimensions();

    return
        dimensioned<Type>
        (
            units_.dimensions(),
            value(x.value())
        );
}


template<class Type>
inline Foam::dimensioned<Type> Foam::Function1s::Dimensioned<Type>::integral
(
    const dimensionedScalar& x1,
    const dimensionedScalar& x2
) const
{
    x1.dimensions() = xUnits_.dimensions();
    x2.dimensions() = xUnits_.dimensions();

    return
        dimensioned<Type>
        (
            units_.dimensions()*xUnits_.dimensions(),
            integral(x1.value(), x2.value())
        );
}


template<class Type>
inline bool Foam::Function1s::Dimensioned<Type>::constant() const
{
    return value_->constant();
}


// ************************************************************************* //
