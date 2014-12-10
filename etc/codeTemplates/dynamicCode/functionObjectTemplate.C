/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "functionObjectTemplate.H"
#include "Time.H"
#include "fvCFD.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(${typeName}FunctionObject, 0);


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode
${localCode}
//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const objectRegistry& ${typeName}FunctionObject::obr() const
{
    return obr_;
}


const fvMesh& ${typeName}FunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

${typeName}FunctionObject::${typeName}FunctionObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool
)
:
    name_(name),
    obr_(obr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

${typeName}FunctionObject::~${typeName}FunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ${typeName}FunctionObject::read(const dictionary& dict)
{
    if (${verbose:-false})
    {
        Info<<"read ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${codeRead}
//}}} end code
}


void ${typeName}FunctionObject::execute()
{
    if (${verbose:-false})
    {
        Info<<"execute ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${codeExecute}
//}}} end code
}


void ${typeName}FunctionObject::end()
{
    if (${verbose:-false})
    {
        Info<<"end ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${codeEnd}
//}}} end code
}


void ${typeName}FunctionObject::timeSet()
{
    if (${verbose:-false})
    {
        Info<<"timeSet ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin codeTime
    ${codeTimeSet}
//}}} end code
}


void ${typeName}FunctionObject::write()
{
    if (${verbose:-false})
    {
        Info<<"write ${typeName} sha1: ${SHA1sum}\n";
    }

//{{{ begin code
    ${code}
//}}} end code
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
