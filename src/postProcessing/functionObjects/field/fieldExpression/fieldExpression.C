/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "fieldExpression.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldExpression, 0);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldExpression::fieldExpression
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const word& fieldName,
    const word& resultName
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_(fieldName),
    resultName_(resultName)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldExpression::~fieldExpression()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldExpression::read(const dictionary& dict)
{
    if (fieldName_ == word::null || dict.found("field"))
    {
        dict.lookup("field") >> fieldName_;
    }

    if (dict.found("result"))
    {
        dict.lookup("result") >> resultName_;
    }

    return true;
}


bool Foam::functionObjects::fieldExpression::write(const bool postProcess)
{
    return fvMeshFunctionObject::writeField(resultName_);
}


// ************************************************************************* //
