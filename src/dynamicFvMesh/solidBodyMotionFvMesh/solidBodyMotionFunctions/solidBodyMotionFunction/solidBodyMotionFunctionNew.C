/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
{
    const word motionType(SBMFCoeffs.lookup("solidBodyMotionFunction"));

    Info<< "Selecting solid-body motion function " << motionType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(motionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "solidBodyMotionFunction::New"
            "("
            "    const dictionary& SBMFCoeffs,"
            "    const Time& runTime"
            ")"
        )   << "Unknown solidBodyMotionFunction type "
            << motionType << nl << nl
            << "Valid solidBodyMotionFunctions are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<solidBodyMotionFunction>(cstrIter()(SBMFCoeffs, runTime));
}


// ************************************************************************* //
