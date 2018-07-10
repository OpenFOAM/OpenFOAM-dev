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

#include "cloudInfo.H"
#include "kinematicCloud.H"
#include "dictionary.H"
#include "PstreamReduceOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cloudInfo, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        cloudInfo,
        dictionary
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::cloudInfo::writeFileHeader(const label i)
{
    writeHeader(file(), "Cloud information");
    writeCommented(file(), "Time");
    writeTabbed(file(), "nParcels");
    writeTabbed(file(), "mass");
    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cloudInfo::cloudInfo
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    regionFunctionObject(name, runTime, dict),
    logFiles(obr_, name)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cloudInfo::~cloudInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::cloudInfo::read(const dictionary& dict)
{
    regionFunctionObject::read(dict);

    logFiles::resetNames(dict.lookup("clouds"));

    Info<< type() << " " << name() << ": ";
    if (names().size())
    {
        Info<< "applying to clouds:" << nl;
        forAll(names(), i)
        {
            Info<< "    " << names()[i] << nl;
        }
        Info<< endl;
    }
    else
    {
        Info<< "no clouds to be processed" << nl << endl;
    }

    return true;
}


bool Foam::functionObjects::cloudInfo::execute()
{
    return true;
}


bool Foam::functionObjects::cloudInfo::write()
{
    logFiles::write();

    forAll(names(), i)
    {
        const word& cloudName = names()[i];

        const kinematicCloud& cloud =
            obr_.lookupObject<kinematicCloud>(cloudName);

        label nParcels = returnReduce(cloud.nParcels(), sumOp<label>());
        scalar massInSystem =
            returnReduce(cloud.massInSystem(), sumOp<scalar>());

        if (Pstream::master())
        {
            writeTime(file(i));
            file(i)
                << tab
                << nParcels << tab
                << massInSystem << endl;
        }
    }

    return true;
}


// ************************************************************************* //
