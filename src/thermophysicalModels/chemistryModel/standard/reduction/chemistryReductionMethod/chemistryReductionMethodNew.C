/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "chemistryReductionMethod.H"
#include "noChemistryReduction.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::chemistryReductionMethod>
Foam::chemistryReductionMethod::New
(
    const dictionary& dict,
    chemistryModels::standard& chemistry
)
{
    if (dict.found("reduction"))
    {
        const dictionary& reductionDict(dict.subDict("reduction"));

        const word methodName(reductionDict.lookup("method"));

        Info<< indentOrNl
            << "Selecting chemistry reduction method " << methodName << endl;

        typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(methodName);

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown " << typeName_() << " type "
                << methodName << nl << nl
                << "Valid " << typeName_() << " types are :" << nl
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalIOError);
        }

        return autoPtr<chemistryReductionMethod>
        (
            cstrIter()(dict, chemistry)
        );
    }
    else
    {
        return autoPtr<chemistryReductionMethod>
        (
            new chemistryReductionMethods::none(dict, chemistry)
        );
    }
}


// ************************************************************************* //
