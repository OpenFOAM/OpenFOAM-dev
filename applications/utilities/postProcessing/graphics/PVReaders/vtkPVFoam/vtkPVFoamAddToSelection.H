/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
