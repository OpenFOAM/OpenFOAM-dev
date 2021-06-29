/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "engineTime.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::engineTime> Foam::engineTime::New
(
    const word& name,
    const fileName& rootPath,
    const fileName& globalCaseName,
    const fileName& caseName,
    const fileName& systemName,
    const fileName& constantName,
    const fileName& dictName
)
{
    IFstream engineDictFile(rootPath/globalCaseName/constantName/dictName);

    dictionary engineDict(engineDictFile);

    const word engineType
    (
        engineDict.lookupOrDefault<word>("engineType", "crankConRod")
    );

    Info<< "Selecting engine type " << engineType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(engineType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown engine type "
            << engineType << nl << nl
            << "Valid engine types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<engineTime>
    (
        cstrIter()
        (
            name,
            rootPath,
            caseName,
            systemName,
            constantName,
            dictName
        )
    );
}


// ************************************************************************* //
