/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "cloud.H"
#include "CloudDerivedField.H"
#include "CloudAverageField.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Type> Foam::cloud::New
(
    const polyMesh& pMesh,
    const word& name,
    const contextType context,
    const dictionary& dict,
    const IOobject::readOption readOption,
    const IOobject::writeOption writeOption
)
{
    return
        autoPtr<Type>
        (
            new Type
            (
                mesh(pMesh, name, readOption, writeOption),
                context,
                dict
            )
        );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, class ... Args>
const Foam::CloudDerivedField<Type>& Foam::cloud::derivedField
(
    const Args& ... args
) const
{
    CloudDerivedField<Type>* ptr = new CloudDerivedField<Type>(args ...);
    derivedFields<Type>().append(ptr);
    return *ptr;
}


template<class Type, class ... Args>
const Foam::CloudAverageField<Type>& Foam::cloud::averageField
(
    const Args& ... args
) const
{
    CloudAverageField<Type>* ptr = new CloudAverageField<Type>(args ...);
    averageFields<Type>().append(ptr);
    return *ptr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
