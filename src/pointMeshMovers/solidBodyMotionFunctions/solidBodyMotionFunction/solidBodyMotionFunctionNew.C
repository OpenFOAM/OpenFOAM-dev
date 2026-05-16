/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidBodyMotionFunction> Foam::solidBodyMotionFunction::New
(
    const dictionary& dict,
    const Time& runTime,
    const word& name
)
{
    const dictionary& coeffDict
    (
        dict.isDict(name)
      ? dict.subDict(name)
      : dict
    );

    const word motionType
    (
        dict.isDict(name)
      ? coeffDict.lookup("type")
      : dict.lookup(name)
    );

    Info<< indentOrNl
        << "Selecting solid-body motion function " << motionType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(motionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(coeffDict)
            << "Unknown solidBodyMotionFunction type "
            << motionType << nl << nl
            << "Valid solidBodyMotionFunctions are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<solidBodyMotionFunction>
    (
        cstrIter()(name, coeffDict, runTime)
    );
}


// ************************************************************************* //
