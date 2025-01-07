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
#include "CloudStateField.H"
#include "CloudDerivedField.H"
#include "CloudAverageField.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, class ... Args>
Foam::CloudStateField<Type>& Foam::cloud::stateField
(
    const Args& ... args
) const
{
    CloudStateField<Type>* ptr = new CloudStateField<Type>(args ...);
    stateFields<Type>().append(ptr);
    return *ptr;
}


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
