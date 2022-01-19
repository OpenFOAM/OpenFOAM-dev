/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2022 OpenFOAM Foundation
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

#include "chemistryTabulationMethod.H"
#include "noChemistryTabulation.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::chemistryTabulationMethod>
Foam::chemistryTabulationMethod::New
(
    const IOdictionary& dict,
    const odeChemistryModel& chemistry
)
{
    if (dict.found("tabulation"))
    {
        const dictionary& tabulationDict(dict.subDict("tabulation"));

        const word methodName(tabulationDict.lookup("method"));

        Info<< "Selecting chemistry tabulation method " << methodName << endl;

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(methodName);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown " << typeName_() << " type " << methodName << endl
                << "Valid " << typeName_() << " types are:"
                << dictionaryConstructorTablePtr_->sortedToc() << endl
                << exit(FatalError);
        }

        return autoPtr<chemistryTabulationMethod>
        (
            cstrIter()(dict, chemistry)
        );
    }
    else
    {
        return autoPtr<chemistryTabulationMethod>
        (
            new chemistryTabulationMethods::none(dict, chemistry)
        );
    }
}


// ************************************************************************* //
