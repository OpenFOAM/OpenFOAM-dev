/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
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

#include "contactAngleModel.H"
#include "constantContactAngle.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::contactAngleModel> Foam::contactAngleModel::New
(
    const dictionary& dict
)
{
    if (dict.isDict("contactAngle"))
    {
        const dictionary& contactAngleDict =
            contactAngleModel::contactAngleDict(dict);

        word contactAngleModelType(contactAngleDict.lookup("type"));

        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(contactAngleModelType);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown contactAngleModel type "
                << contactAngleModelType << endl << endl
                << "Valid contactAngleModel types are : " << endl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return cstrIter()(contactAngleDict);
    }
    else
    {
        return autoPtr<contactAngleModel>
        (
            new contactAngleModels::constant(dict)
        );
    }
}


// ************************************************************************* //
