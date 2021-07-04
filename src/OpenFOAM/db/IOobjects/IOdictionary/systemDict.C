/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "systemDict.H"
#include "IOdictionary.H"
#include "Time.H"

// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

Foam::IOobject Foam::systemDictIO
(
    const word& dictName,
    const argList& args,
    const objectRegistry& ob,
    const word& regionName
)
{
    fileName dictPath = dictName;

    if (args.optionFound("dict"))
    {
        dictPath = args["dict"];

        if
        (
            isDir
            (
                dictPath.isAbsolute()
              ? dictPath
              : ob.time().globalPath()/dictPath
            )
        )
        {
            dictPath = dictPath/dictName;
        }
    }

    Info<< "Reading " << dictPath << nl << endl;

    if (args.optionFound("dict") && !dictPath.isName())
    {
        return
            IOobject
            (
                dictPath.isAbsolute()
              ? dictPath
              : ob.time().globalPath()/dictPath,
                ob,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            );
    }
    else
    {
        return
            IOobject
            (
                dictPath,
                ob.time().system(),
                regionName == polyMesh::defaultRegion ? word::null : regionName,
                ob,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            );
    }
}


Foam::IOdictionary Foam::systemDict
(
    const word& dictName,
    const argList& args,
    const objectRegistry& ob,
    const word& regionName
)
{
    return IOdictionary(systemDictIO(dictName, args, ob, regionName));
}


// ************************************************************************* //
