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
    const word& regionName,
    const word& systemName
)
{
    fileName dictPath = fileName::null;

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

    if (dictPath.size())
    {
        Info<< "Reading " << dictPath << nl << endl;

        return
            IOobject
            (
                dictPath,
                ob.time(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            );
    }
    else
    {
        Info<< "Reading " << dictName << nl << endl;

        return
            IOobject
            (
                dictName,
                systemName == word::null ? ob.time().system() : systemName,
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
    const word& regionName,
    const word& systemName
)
{
    return
        IOdictionary
        (
            systemDictIO(dictName, args, ob, regionName, systemName)
        );
}


// ************************************************************************* //
