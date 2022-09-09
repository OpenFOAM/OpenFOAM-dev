/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "patchToPatchFvPatchFieldMapper.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::patchToPatchFvPatchFieldMapper::map
(
    Field<Type>& f,
    const Field<Type>& mapF
) const
{
    f = forward_ ? pToP_.srcToTgt(mapF) : pToP_.tgtToSrc(mapF);
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatchFvPatchFieldMapper::map
(
    const Field<Type>& mapF
) const
{
    tmp<Field<Type>> tf(new Field<Type>(size_));
    map(tf.ref(), mapF);
    return tf;
}


template<class Type>
void Foam::patchToPatchFvPatchFieldMapper::operator()
(
    Field<Type>& f,
    const tmp<Field<Type>>& tmapF
) const
{
    map(f, tmapF());
    tmapF.clear();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchToPatchFvPatchFieldMapper::operator()
(
    const tmp<Field<Type>>& tmapF
) const
{
    tmp<Foam::Field<Type>> tf(map(tmapF()));
    tmapF.clear();
    return tf;
}


// ************************************************************************* //
