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

#include "reconstruct.H"
#include "fvcReconstruct.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(reconstruct, 0);
    addToRunTimeSelectionTable(functionObject, reconstruct, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::reconstruct::calcReconstruction()
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (foundObject<SurfaceFieldType>(fieldName_))
    {
        return store
        (
            resultName_,
            fvc::reconstruct(lookupObject<SurfaceFieldType>(fieldName_)),
            mesh_.changing() && mesh_.solution().cache(resultName_)
        );
    }
    else
    {
        return false;
    }
}


bool Foam::functionObjects::reconstruct::calc()
{
    bool processed = false;

    processed = processed || calcReconstruction<scalar>();
    processed = processed || calcReconstruction<vector>();

    if (!processed)
    {
        cannotFindObject(fieldName_);
    }

    return processed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::reconstruct::reconstruct
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::reconstruct::~reconstruct()
{}


// ************************************************************************* //
