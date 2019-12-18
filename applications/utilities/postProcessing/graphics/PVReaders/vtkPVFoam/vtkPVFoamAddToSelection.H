/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#ifndef vtkPVFoamAddToSelection_H
#define vtkPVFoamAddToSelection_H

// OpenFOAM includes
#include "IOobjectList.H"
#include "SortableList.H"
#include "surfaceFields.H"

// VTK includes
#include "vtkDataArraySelection.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::label Foam::vtkPVFoam::addToSelection
(
    vtkDataArraySelection *select,
    const IOobjectList& objectLst,
    const string& suffix
)
{
    SortableList<word> names(objectLst.names(Type::typeName));

    forAll(names, nameI)
    {
        if (suffix.size())
        {
            select->AddArray
            (
                (names[nameI] + suffix).c_str()
            );
        }
        else
        {
            select->AddArray
            (
                (names[nameI]).c_str()
            );
        }
    }

    return names.size();
}


template<template<class> class patchType, class meshType>
void Foam::vtkPVFoam::addFieldsToSelection
(
    vtkDataArraySelection *select,
    const IOobjectList& objects,
    const string& suffix
)
{
    addToSelection<GeometricField<scalar, patchType, meshType>>
    (
        select,
        objects,
        suffix
    );

    addToSelection<GeometricField<vector, patchType, meshType>>
    (
        select,
        objects,
        suffix
    );

    addToSelection<GeometricField<sphericalTensor, patchType, meshType>>
    (
        select,
        objects,
        suffix
    );

    addToSelection<GeometricField<symmTensor, patchType, meshType>>
    (
        select,
        objects,
        suffix
    );

    addToSelection<GeometricField<tensor, patchType, meshType>>
    (
        select,
        objects,
        suffix
    );
}


template<template<class> class patchType, class meshType>
void Foam::vtkPVFoam::addInternalFieldsToSelection
(
    vtkDataArraySelection *select,
    const IOobjectList& objects,
    const string& suffix
)
{
    addToSelection
    <
        typename
        GeometricField<scalar, patchType, meshType>::Internal
    >
    (
        select,
        objects,
        suffix
    );

    addToSelection
    <
        typename
        GeometricField<vector, patchType, meshType>::Internal
    >
    (
        select,
        objects,
        suffix
    );

    addToSelection
    <
        typename
        GeometricField<sphericalTensor, patchType, meshType>::Internal
    >
    (
        select,
        objects,
        suffix
    );

    addToSelection
    <
        typename
        GeometricField<symmTensor, patchType, meshType>::Internal
    >
    (
        select,
        objects,
        suffix
    );

    addToSelection
    <
        typename
        GeometricField<tensor, patchType, meshType>::Internal
    >
    (
        select,
        objects,
        suffix
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
