/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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
void Foam::patchToPatchFvPatchFieldMapper::operator()
(
    Field<Type>& f,
    const tmp<Field<Type>>& tmapF
) const
{
    operator()(f, tmapF());
    tmapF.clear();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchToPatchFvPatchFieldMapper::operator()
(
    const tmp<Field<Type>>& tmapF
) const
{
    tmp<Foam::Field<Type>> tf(operator()(tmapF()));
    tmapF.clear();
    return tf;
}


template<class Type>
void Foam::patchToPatchLeftOverFvPatchFieldMapper::map
(
    Field<Type>& f,
    const Field<Type>& mapF
) const
{
    f = pToP_.srcToTgt(mapF, f);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchToPatchLeftOverFvPatchFieldMapper::map
(
    const Field<Type>& mapF
) const
{
    FatalErrorInFunction
        << "Not a valid operation for this mapper, which should only be "
        << "used for modifying an existing, valid, field"
        << exit(FatalError);
    return tmp<Field<Type>>(nullptr);
}


template<class Type>
void Foam::patchToPatchNormalisedFvPatchFieldMapper::map
(
    Field<Type>& f,
    const Field<Type>& mapF
) const
{
    f = pToP_.srcToTgt(mapF);

    pS_.stabilise(f);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::patchToPatchNormalisedFvPatchFieldMapper::map
(
    const Field<Type>& mapF
) const
{
    tmp<Field<Type>> tf(new Field<Type>());
    map(tf.ref(), mapF);
    return tf;
}


// ************************************************************************* //
