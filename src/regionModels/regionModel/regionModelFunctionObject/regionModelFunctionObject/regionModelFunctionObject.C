/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    regionModel& owner
)
:
    dict_(dictionary::null),
    owner_(owner),
    modelType_("modelType")
{}


Foam::regionModels::regionModelFunctionObject::regionModelFunctionObject
(
    const dictionary& dict,
    regionModel& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    modelType_(type)
{}


Foam::regionModels::regionModelFunctionObject::regionModelFunctionObject
(
    const regionModelFunctionObject& rmfo
)
:
    dict_(rmfo.dict_),
    owner_(rmfo.owner_),
    modelType_(rmfo.modelType_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionModels::regionModelFunctionObject::~regionModelFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionModels::regionModelFunctionObject::preEvolveRegion()
{
    // do nothing
}


void Foam::regionModels::regionModelFunctionObject::postEvolveRegion()
{
    if (owner_.regionMesh().time().outputTime())
    {
        write();
    }
}


void Foam::regionModels::regionModelFunctionObject::write() const
{
    // do nothing
}


// ************************************************************************* //
