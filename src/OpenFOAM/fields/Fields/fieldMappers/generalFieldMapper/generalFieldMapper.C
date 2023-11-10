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

#include "generalFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::generalFieldMapper::map
(
    Field<Type>& f,
    const Field<Type>& mapF
) const
{
    if (direct())
    {
        if (notNull(directAddressing()) && directAddressing().size())
        {
            f.map(mapF, directAddressing());
        }
        else
        {
            f.setSize(0);
        }
    }
    else
    {
        if (notNull(addressing()) && addressing().size())
        {
            f.map(mapF, addressing(), weights());
        }
        else
        {
            f.setSize(0);
        }
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::generalFieldMapper::map
(
    const Field<Type>& mapF
) const
{
    tmp<Field<Type>> tf
    (
        new Field<Type>
        (
            direct() ? directAddressing().size() : addressing().size()
        )
    );
    map(tf.ref(), mapF);
    return tf;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelUList& Foam::generalFieldMapper::directAddressing() const
{
    FatalErrorInFunction
        << "attempt to access null direct addressing"
        << abort(FatalError);

    return labelUList::null();
}


const Foam::labelListList& Foam::generalFieldMapper::addressing() const
{
    FatalErrorInFunction
        << "attempt to access null interpolation addressing"
        << abort(FatalError);

    return labelListList::null();
}


const Foam::scalarListList& Foam::generalFieldMapper::weights() const
{
    FatalErrorInFunction
        << "attempt to access null interpolation weights"
        << abort(FatalError);

    return scalarListList::null();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

FOR_ALL_FIELD_TYPES(IMPLEMENT_FIELD_MAPPER_OPERATOR, generalFieldMapper)


IMPLEMENT_FIELD_MAPPER_OPERATOR(label, generalFieldMapper)


// ************************************************************************* //
