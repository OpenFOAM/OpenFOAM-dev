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

#include "LagrangianSubFieldsFwd.H"
#include "LagrangianSubMesh.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
Foam::SubList<Type> Foam::LagrangianSubMesh::sub
(
    const List<Type>& list
) const
{
    return SubList<Type>(list, size(), start());
}


template<class Type>
Foam::SubField<Type> Foam::LagrangianSubMesh::sub
(
    const Field<Type>& field
) const
{
    return SubField<Type>(sub(static_cast<const List<Type>&>(field)));
}


template<class Type, template<class> class PrimitiveField>
Foam::tmp<Foam::LagrangianSubSubField<Type>>
Foam::LagrangianSubMesh::sub
(
    const DimensionedField<Type, LagrangianMesh, PrimitiveField>& field
) const
{
    const word subFieldName = field.name() + ':' + name(group());

    return tmp<LagrangianSubSubField<Type>>
    (
        new LagrangianSubSubField<Type>
        (
            IOobject
            (
                subFieldName,
                field.instance(),
                field.local(),
                field.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                field.db().cacheTemporaryObject(subFieldName)
            ),
            *this,
            field.dimensions(),
            sub(static_cast<const Field<Type>&>(field))
        )
    );
}


// ************************************************************************* //
