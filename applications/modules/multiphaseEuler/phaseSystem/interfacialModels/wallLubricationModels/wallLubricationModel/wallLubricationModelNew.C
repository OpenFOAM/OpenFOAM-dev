/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2026 OpenFOAM Foundation
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

#include "wallLubricationModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::wallLubricationModel> Foam::wallLubricationModel::New
(
    const dictionary& modelDict,
    const phaseInterface& interface,
    const bool outer
)
{
    const word wallLubricationModelType(modelDict.lookup("type"));

    Info<< indentOrNl << "Selecting " << typeName << ' '
        << wallLubricationModelType;
    if (outer) Info<< " for " << interface.name();
    Info<< endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(wallLubricationModelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(modelDict)
            << "Unknown wallLubricationModel type "
            << wallLubricationModelType << endl << endl
            << "Valid wallLubricationModel types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    printDictionary print(modelDict);

    return cstrIter()(modelDict, interface);
}


Foam::autoPtr<Foam::wallLubricationModel> Foam::wallLubricationModel::New
(
    const UPtrList<const dictionary>& subDicts,
    const phaseInterface& interface
)
{
    return New(modelSubDict<wallLubricationModel>(subDicts), interface, true);
}


// ************************************************************************* //
