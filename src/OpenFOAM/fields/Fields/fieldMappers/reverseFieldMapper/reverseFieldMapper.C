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

#include "reverseFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::reverseFieldMapper::map
(
    Field<Type>& f,
    const Field<Type>& mapF
) const
{
    f.rmap(mapF, addressing_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::reverseFieldMapper::map
(
    const Field<Type>& mapF
) const
{
    tmp<Field<Type>> tf(new Field<Type>(max(addressing_) + 1));
    map(tf.ref(), mapF);
    return tf;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

FOR_ALL_FIELD_TYPES(IMPLEMENT_FIELD_MAPPER_OPERATOR, reverseFieldMapper)


IMPLEMENT_FIELD_MAPPER_OPERATOR(label, reverseFieldMapper)


// ************************************************************************* //
