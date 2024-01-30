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

#include "forwardOrAssignFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::forwardOrAssignFieldMapper::unmappedError() const
{
    FatalErrorInFunction
        << "A field could not be mapped because a value was not provided for "
        << "unmapped elements" << exit(FatalError);
}


void Foam::forwardOrAssignPatchFieldMapper::unmappedError() const
{
    FatalErrorInFunction
        << "The " << fieldTypeName_ << " " << fieldBaseTypeName_
        << " for patch " << patchName_ << " of field " << internalFieldName_
        << " could not be mapped because a value was not provided for "
        << "unmapped elements " << exit(FatalError);
}


template<class Type>
void Foam::forwardOrAssignFieldMapper::map
(
    Field<Type>& f,
    const Field<Type>& mapF
) const
{
    if (hasUnmapped_)
    {
        unmappedError();
    }

    f.map(mapF, addressing_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::forwardOrAssignFieldMapper::map
(
    const Field<Type>& mapF
) const
{
    tmp<Field<Type>> tf(new Field<Type>(addressing_.size()));
    map(tf.ref(), mapF);
    return tf;
}


template<class Type>
void Foam::forwardOrAssignFieldMapper::mapOrAssign
(
    Field<Type>& f,
    const Field<Type>& mapF,
    const Type& unmappedVal
) const
{
    if (hasUnmapped_)
    {
        f = unmappedVal;
    }

    f.map(mapF, addressing_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::forwardOrAssignFieldMapper::mapOrAssign
(
    const Field<Type>& mapF,
    const Type& unmappedVal
) const
{
    tmp<Field<Type>> tf(new Field<Type>(addressing_.size()));
    mapOrAssign(tf.ref(), mapF, unmappedVal);
    return tf;
}


template<class Type>
void Foam::forwardOrAssignFieldMapper::mapOrAssign
(
    Field<Type>& f,
    const Field<Type>& mapF,
    const FieldFunctor<Type>& unmappedFunc
) const
{
    if (hasUnmapped_)
    {
        f = unmappedFunc();
    }

    f.map(mapF, addressing_);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::forwardOrAssignFieldMapper::mapOrAssign
(
    const Field<Type>& mapF,
    const FieldFunctor<Type>& unmappedFunc
) const
{
    tmp<Field<Type>> tf(new Field<Type>(addressing_.size()));
    mapOrAssign(tf.ref(), mapF, unmappedFunc);
    return tf;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FIELD_MAPPER_MAP_OPERATOR,
    forwardOrAssignFieldMapper
)


IMPLEMENT_FIELD_MAPPER_MAP_OPERATOR(label, forwardOrAssignFieldMapper)


FOR_ALL_FIELD_TYPES
(
    IMPLEMENT_FIELD_MAPPER_MAP_OR_ASSIGN_OPERATOR,
    forwardOrAssignFieldMapper
)


IMPLEMENT_FIELD_MAPPER_MAP_OR_ASSIGN_OPERATOR(label, forwardOrAssignFieldMapper)


// ************************************************************************* //
