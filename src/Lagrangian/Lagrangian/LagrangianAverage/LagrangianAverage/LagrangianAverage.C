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

#include "LagrangianAverage.H"
#include "LagrangianFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianAverage<Type>::LagrangianAverage
(
    const word& name,
    const LagrangianMesh& mesh,
    const dimensionSet& dimensions
)
:
    name_(name),
    mesh_(mesh),
    dimensions_(dimensions)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Type>
template<template<class> class WeightPF, template<class> class PsiPF>
Foam::autoPtr<Foam::LagrangianAverage<Type>>
Foam::LagrangianAverage<Type>::New
(
    const word& type,
    const word& name,
    const DimensionedField<scalar, LagrangianMesh, WeightPF>& weight,
    const DimensionedField<Type, LagrangianMesh, PsiPF>& psi
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Lagrangian average type " << type
            << " for field " << psi.name() << nl << nl
            << "Valid average types : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    autoPtr<LagrangianAverage<Type>> result =
        cstrIter()
        (
            name,
            psi.mesh(),
            psi.dimensions(),
            NullObjectRef<Field<scalar>>()
        );

    result->add
    (
        weight.mesh().subAll().sub(weight),
        psi.mesh().subAll().sub(psi),
        false
    );

    return result;
}


template<class Type>
template<class CellMesh, template<class> class PsiPF>
typename Foam::LagrangianAverageNewReturnType<CellMesh, Type>::type
Foam::LagrangianAverage<Type>::New
(
    const word& type,
    const word& name,
    const DimensionedField<scalar, CellMesh>& cellWeightSum,
    const DimensionedField<Type, LagrangianMesh, PsiPF>& weightPsi
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Lagrangian average type " << type
            << " for field " << weightPsi.name() << nl << nl
            << "Valid average types : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    autoPtr<LagrangianAverage<Type>> result =
        cstrIter()
        (
            name,
            weightPsi.mesh(),
            weightPsi.dimensions()/cellWeightSum.dimensions(),
            cellWeightSum
        );

    result->add
    (
        NullObjectRef<LagrangianSubSubField<scalar>>(),
        weightPsi.mesh().subAll().sub(weightPsi),
        false
    );

    return result;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::LagrangianAverage<Type>::~LagrangianAverage()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void Foam::LagrangianAverage<Type>::remove
(
    const LagrangianSubSubField<Type>& weightPsi
)
{
    remove(NullObjectRef<LagrangianSubSubField<scalar>>(), weightPsi);
}


template<class Type>
void Foam::LagrangianAverage<Type>::add
(
    const LagrangianSubSubField<Type>& weightPsi,
    const bool cache
)
{
    add(NullObjectRef<LagrangianSubSubField<scalar>>(), weightPsi, cache);
}


template<class Type>
void Foam::LagrangianAverage<Type>::correct
(
    const LagrangianSubSubField<Type>& weightPsi,
    const bool cache
)
{
    correct(NullObjectRef<LagrangianSubSubField<scalar>>(), weightPsi, cache);
}


template<class Type>
Foam::tmp<Foam::LagrangianSubField<Type>>
Foam::LagrangianAverage<Type>::interpolate
(
    const LagrangianSubMesh& subMesh
) const
{
    tmp<LagrangianSubField<Type>> tresult =
        LagrangianSubField<Type>::New(name_, subMesh, dimensions_);

    interpolate(tresult.ref());

    return tresult;
}


// ************************************************************************* //
