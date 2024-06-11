/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "distribution.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Base, class Derived>
template<class Type>
Type Foam::FieldDistribution<Base, Derived>::sample() const
{
    Type value;
    for (direction i = 0; i < pTraits<Type>::nComponents; ++ i)
    {
        setComponent(value, i) =
            static_cast<const Derived&>(*this).sample();
    }
    return value;
}


template<class Base, class Derived>
Foam::tmp<Foam::scalarField> Foam::FieldDistribution<Base, Derived>::sample
(
    const label n
) const
{
    tmp<scalarField> tResult(new scalarField(n));
    scalarField& result = tResult.ref();

    forAll(result, i)
    {
        result[i] =
            static_cast<const Derived&>(*this).sample();
    }

    return tResult;
}


// ************************************************************************* //
