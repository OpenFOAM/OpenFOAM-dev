/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "vtkPVFoam.H"

// OpenFOAM includes
#include "IOobjectList.H"
#include "SortableList.H"

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


template<class meshType>
void Foam::vtkPVFoam::addFieldsToSelection
(
    vtkDataArraySelection *select,
    const IOobjectList& objects,
    const string& suffix
)
{
    #define ADD_TO_SELECTION(Type, nullArg)                                    \
        addToSelection<GeometricField<Type, meshType>>                         \
        (                                                                      \
            select,                                                            \
            objects,                                                           \
            suffix                                                             \
        );
    FOR_ALL_FIELD_TYPES(ADD_TO_SELECTION);
    #undef ADD_TO_SELECTION
}


template<class meshType>
void Foam::vtkPVFoam::addInternalFieldsToSelection
(
    vtkDataArraySelection *select,
    const IOobjectList& objects,
    const string& suffix
)
{
    #define ADD_TO_SELECTION(Type, nullArg)                                    \
        addToSelection<typename GeometricField<Type, meshType>::Internal>      \
        (                                                                      \
            select,                                                            \
            objects,                                                           \
            suffix                                                             \
        );
    FOR_ALL_FIELD_TYPES(ADD_TO_SELECTION);
    #undef ADD_TO_SELECTION
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
