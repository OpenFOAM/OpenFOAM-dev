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

#include "helpFunctionObject.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace helpTypes
    {
        defineTypeNameAndDebug(helpFunctionObject, 0);
        addNamedToRunTimeSelectionTable
        (
            helpType,
            helpFunctionObject,
            dictionary,
            functionObject
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::helpTypes::helpFunctionObject::helpFunctionObject()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::helpTypes::helpFunctionObject::~helpFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::helpTypes::helpFunctionObject::init()
{
    helpType::init();

    argList::validArgs.append("functionObject");
}


void Foam::helpTypes::helpFunctionObject::execute
(
    const argList& args,
    const fvMesh& mesh
)
{
    word function(word::null);

    if (args.optionReadIfPresent("browse", function))
    {
        displayDoc(function, ".*[fF]unctionObject.*", true);
    }
    else
    {
        displayDocOptions(".*[fF]unctionObject.*", true);
    }
}


// ************************************************************************* //
