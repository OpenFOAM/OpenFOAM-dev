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

#include "forwardFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::forwardFieldMapper::map
(
    Field<Type>& f,
    const Field<Type>& mapF
) const
{
    if (notNull(addressing_) && addressing_.size())
    {
        f.map(mapF, addressing_);
    }
    else
    {
        f.setSize(0);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::forwardFieldMapper::map
(
    const Field<Type>& mapF
) const
{
    if (notNull(addressing_) && addressing_.size())
    {
        tmp<Field<Type>> tf(new Field<Type>(addressing_.size()));
        map(tf.ref(), mapF);
        return tf;
    }
    else
    {
        return tmp<Field<Type>>(new Field<Type>(0));
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

FOR_ALL_FIELD_TYPES(IMPLEMENT_FIELD_MAPPER_OPERATOR, forwardFieldMapper)


IMPLEMENT_FIELD_MAPPER_OPERATOR(label, forwardFieldMapper)


// ************************************************************************* //
