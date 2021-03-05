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

Foam::dictionary Foam::systemDict
(
    const word& dictName,
    const argList& args,
    const objectRegistry& mesh
)
{
    fileName dictPath = "";
    if (args.optionFound("dict"))
    {
        dictPath = args["dict"];
        if (isDir(dictPath))
        {
            dictPath = dictPath / dictName;
        }
    }

    if (dictPath.size())
    {
        Info<< "Reading " << dictPath << nl << endl;

        return IOdictionary
        (
            IOobject
            (
                dictPath,
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
    }
    else
    {
        Info<< "Reading " << dictName << nl << endl;

        return IOdictionary
        (
            IOobject
            (
                dictName,
                mesh.time().system(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
    }
}


// ************************************************************************* //
