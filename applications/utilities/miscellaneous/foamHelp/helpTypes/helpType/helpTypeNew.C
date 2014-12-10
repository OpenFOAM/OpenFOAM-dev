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

#include "helpType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::helpType> Foam::helpType::New
(
    const word& helpTypeName
)
{
    Info<< "Selecting helpType " << helpTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(helpTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        // special treatment for -help
        // exit without stack trace
        if (helpTypeName == "-help")
        {
            FatalErrorIn("helpType::New()")
                << "Valid helpType selections are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }
        else
        {
            FatalErrorIn("helpType::New()")
                << "Unknown helpType type " << helpTypeName << nl
                << "Valid helpType selections are:" << nl
                << dictionaryConstructorTablePtr_->sortedToc() << nl
                << abort(FatalError);
        }
    }

    return autoPtr<helpType>(cstrIter()());
}


// ************************************************************************* //
