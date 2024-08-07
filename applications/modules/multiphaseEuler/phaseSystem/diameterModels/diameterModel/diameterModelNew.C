/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::diameterModel> Foam::diameterModel::New
(
    const dictionary& dict,
    const phaseModel& phase
)
{
    word diameterModelType
    (
        dict.lookup("diameterModel")
    );

    Info << "Selecting diameterModel for phase "
        << phase.name()
        << ": "
        << diameterModelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(diameterModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
           << "Unknown diameterModel type "
           << diameterModelType << endl << endl
           << "Valid diameterModel types are : " << endl
           << dictionaryConstructorTablePtr_->sortedToc()
           << exit(FatalIOError);
    }

    return cstrIter()
    (
        dict.optionalSubDict(diameterModelType + "Coeffs"),
        phase
    );
}


// ************************************************************************* //
