/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

#include "objectFunction1.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class PrimitiveType>
Foam::objectFunction1::objectFunction1
(
    const word& name,
    const dictionary& dict,
    const type<PrimitiveType>&
)
:
    autoPtr<Function1<PrimitiveType>>
    (
        Function1<PrimitiveType>::New(name, dict).ptr()
    )
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<template<class> class ObjectType>
Foam::autoPtr<Foam::objectFunction1> Foam::objectFunction1::New
(
    const word& name,
    const dictionary& dict,
    const word& objectName,
    const objectRegistry& db,
    const bool error
)
{
    autoPtr<objectFunction1> ptr
    (
        db.foundObject<ObjectType<scalar>>(objectName)
      ? new objectFunction1(name, dict, type<scalar>())
      : db.foundObject<ObjectType<vector>>(objectName)
      ? new objectFunction1(name, dict, type<vector>())
      : db.foundObject<ObjectType<symmTensor>>(objectName)
      ? new objectFunction1(name, dict, type<symmTensor>())
      : db.foundObject<ObjectType<sphericalTensor>>(objectName)
      ? new objectFunction1(name, dict, type<sphericalTensor>())
      : db.foundObject<ObjectType<tensor>>(objectName)
      ? new objectFunction1(name, dict, type<tensor>())
      : nullptr
    );

    if (error && !ptr.valid())
    {
        // Spit lookup error
        db.lookupObject<regIOobject>(objectName);
    }

    return ptr;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class PrimitiveType>
PrimitiveType Foam::objectFunction1::value
(
    const scalar x
) const
{
    return autoPtr<Function1<PrimitiveType>>::operator*().value(x);
}


template<class PrimitiveType>
Foam::tmp<Foam::Field<PrimitiveType>> Foam::objectFunction1::value
(
    const scalarField& x
) const
{
    return autoPtr<Function1<PrimitiveType>>::operator*().value(x);
}


template<class PrimitiveType>
PrimitiveType Foam::objectFunction1::integral
(
    const scalar x1,
    const scalar x2
) const
{
    return autoPtr<Function1<PrimitiveType>>::operator*().integral(x1, x2);
}


template<class PrimitiveType>
Foam::tmp<Foam::Field<PrimitiveType>> Foam::objectFunction1::integral
(
    const scalarField& x1,
    const scalarField& x2
) const
{
    return autoPtr<Function1<PrimitiveType>>::operator*().integral(x1, x2);
}


// ************************************************************************* //
