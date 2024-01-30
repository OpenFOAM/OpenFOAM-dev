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

#include "fieldMapper.H"

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::fieldMapper::operator()
(
    Field<Type>& f,
    const tmp<Field<Type>>& tmapF
) const
{
    operator()(f, tmapF());
    tmapF.clear();
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fieldMapper::operator()
(
    const tmp<Field<Type>>& tmapF
) const
{
    tmp<Foam::Field<Type>> tf(operator()(tmapF()));
    tmapF.clear();
    return tf;
}


template<class Type, class FieldOp>
void Foam::fieldMapper::operator()
(
    Field<Type>& f,
    const Field<Type>& mapF,
    const FieldOp& unmappedOp
) const
{
    operator()
    (
        f,
        mapF,
        static_cast<const FieldFunctor<Type>&>
        (
            FieldOpFunctor<Type, FieldOp>(unmappedOp)
        )
    );
}


template<class Type, class FieldOp>
Foam::fieldMapper::TmpFieldTypeIfFieldOp<Type, FieldOp>
Foam::fieldMapper::operator()
(
    const Field<Type>& mapF,
    const FieldOp& unmappedOp
) const
{
    operator()
    (
        mapF,
        static_cast<const FieldFunctor<Type>&>
        (
            FieldOpFunctor<Type, FieldOp>(unmappedOp)
        )
    );
}


// ************************************************************************* //
