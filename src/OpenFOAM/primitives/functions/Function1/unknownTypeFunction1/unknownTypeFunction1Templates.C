/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "unknownTypeFunction1.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::unknownTypeFunction1::build(const unitConversion& valueUnits) const
{
    if (!functionPtr_.autoPtr<Function1<Type>>::valid())
    {
        functionPtr_.autoPtr<Function1<Type>>::set
        (
            Function1<Type>::New
            (
                name_,
                xUnits_,
                valueUnits,
                topDict_.scopedDict(topDictKeyword_)
            ).ptr()
        );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::unknownTypeFunction1::setValueUnits
(
    const unitConversion& valueUnits
) const
{
    build<Type>(valueUnits);
}


template<class Type>
Type Foam::unknownTypeFunction1::value
(
    const scalar x
) const
{
    build<Type>(unitAny);

    return functionPtr_.autoPtr<Function1<Type>>::operator*().value(x);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::unknownTypeFunction1::value
(
    const scalarField& x
) const
{
    build<Type>(unitAny);

    return functionPtr_.autoPtr<Function1<Type>>::operator*().value(x);
}


template<class Type>
Type Foam::unknownTypeFunction1::integral
(
    const scalar x1,
    const scalar x2
) const
{
    build<Type>(unitAny);

    return functionPtr_.autoPtr<Function1<Type>>::operator*().integral(x1, x2);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::unknownTypeFunction1::integral
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    build<Type>(unitAny);

    return functionPtr_.autoPtr<Function1<Type>>::operator*().integral(x1, x2);
}


// ************************************************************************* //
