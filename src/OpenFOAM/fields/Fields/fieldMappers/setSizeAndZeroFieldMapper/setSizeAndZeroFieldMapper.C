/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2023 OpenFOAM Foundation
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

#include "setSizeAndZeroFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::setSizeAndZeroFieldMapper::map
(
    Field<Type>& f,
    const Field<Type>& mapF
) const
{
    f.setSize(size_);
    f = pTraits<Type>::zero;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::setSizeAndZeroFieldMapper::map
(
    const Field<Type>& mapF
) const
{
    return tmp<Field<Type>>(new Field<Type>(size_, pTraits<Type>::zero));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FIELD_MAPPER_MAP_OPERATOR,
    setSizeAndZeroFieldMapper
)


IMPLEMENT_FIELD_MAPPER_MAP_OPERATOR(label, setSizeAndZeroFieldMapper)


// ************************************************************************* //
