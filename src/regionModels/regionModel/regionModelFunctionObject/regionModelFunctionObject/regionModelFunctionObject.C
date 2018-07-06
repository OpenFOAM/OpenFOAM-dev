/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "regionModelFunctionObject.H"
#include "regionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
    defineTypeNameAndDebug(regionModelFunctionObject, 0);
    defineRunTimeSelectionTable(regionModelFunctionObject, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionModels::regionModelFunctionObject::regionModelFunctionObject
(
    regionModel& region
)
:
    dict_(dictionary::null),
    regionModel_(region),
    modelType_("modelType")
{}


Foam::regionModels::regionModelFunctionObject::regionModelFunctionObject
(
    const dictionary& dict,
    regionModel& region,
    const word& type
)
:
    dict_(dict),
    regionModel_(region),
    modelType_(type)
{}


Foam::regionModels::regionModelFunctionObject::regionModelFunctionObject
(
    const regionModelFunctionObject& rmfo
)
:
    dict_(rmfo.dict_),
    regionModel_(rmfo.regionModel_),
    modelType_(rmfo.modelType_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::regionModelFunctionObject::~regionModelFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionModels::regionModelFunctionObject::preEvolveRegion()
{}


void Foam::regionModels::regionModelFunctionObject::postEvolveRegion()
{
    if (regionModel_.regionMesh().time().writeTime())
    {
        write();
    }
}


void Foam::regionModels::regionModelFunctionObject::write() const
{}


// ************************************************************************* //
