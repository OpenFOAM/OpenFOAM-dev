/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2026 OpenFOAM Foundation
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

#include "interfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::interfaceCompositionModel>
Foam::interfaceCompositionModel::New
(
    const dictionary& modelDict,
    const phaseInterface& interface,
    const bool outer
)
{
    const word interfaceCompositionModelType(modelDict.lookup("type"));

    Info<< indentOrNl << "Selecting " << typeName << ' '
        << interfaceCompositionModelType;
    if (outer) Info<< " for " << interface.name();
    Info<< endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(interfaceCompositionModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(modelDict)
            << "Unknown interfaceCompositionModel type "
            << interfaceCompositionModelType << endl << endl
            << "Valid interfaceCompositionModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    printDictionary print(modelDict);

    return cstrIter()(modelDict, interface);
}


Foam::autoPtr<Foam::interfaceCompositionModel>
Foam::interfaceCompositionModel::New
(
    const UPtrList<const dictionary>& subDicts,
    const phaseInterface& interface
)
{
    return
        New
        (
            modelSubDict<interfaceCompositionModel>(subDicts),
            interface,
            true
        );
}


// ************************************************************************* //
