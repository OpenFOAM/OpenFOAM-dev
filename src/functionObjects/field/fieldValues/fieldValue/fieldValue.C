/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fieldValue.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldValue, 0);
    defineRunTimeSelectionTable(fieldValue, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValue::fieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const word& valueType
)
:
    writeFiles(name, runTime, dict, name),
    dict_(dict),
    regionName_(word::null),
    resultDict_(fileName("name"), dictionary::null)
{
    read(dict);
    resetName(valueType);
}


Foam::functionObjects::fieldValue::fieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const word& valueType
)
:
    writeFiles(name, obr, dict, name),
    dict_(dict),
    regionName_(word::null),
    resultDict_(fileName("name"), dictionary::null)
{
    read(dict);
    resetName(valueType);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValue::~fieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValue::read(const dictionary& dict)
{
    dict_ = dict;
    writeFiles::read(dict);

    dict.lookup("fields") >> fields_;
    dict.lookup("writeFields") >> writeFields_;

    return true;
}


bool Foam::functionObjects::fieldValue::execute()
{
    return true;
}


bool Foam::functionObjects::fieldValue::write()
{
    writeFiles::write();

    Log << type() << " " << name() << " write:" << nl;

    return true;
}


// ************************************************************************* //
