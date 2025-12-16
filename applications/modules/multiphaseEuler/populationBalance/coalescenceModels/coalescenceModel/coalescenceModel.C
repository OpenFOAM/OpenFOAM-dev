/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "coalescenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace populationBalance
{
    defineTypeNameAndDebug(coalescenceModel, 0);
    defineRunTimeSelectionTable(coalescenceModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::populationBalance::coalescenceModel>
Foam::populationBalance::coalescenceModel::New
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
{
    const bool haveModelDict = dict.isDict(typeName);

    word modelType;
    const dictionary* modelDictPtr = nullptr;
    if (haveModelDict)
    {
        modelDictPtr = &dict.subDict(typeName);
        modelType = modelDictPtr->lookup<word>("type");
    }
    else
    {
        modelType = dict.lookup<word>(typeName);
        modelDictPtr = &dict.optionalSubDict(modelType + "Coeffs");
    }
    const dictionary& modelDict = *modelDictPtr;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type "
            << modelType << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(popBal, modelDict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::populationBalance::coalescenceModel::coalescenceModel
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    popBal_(popBal)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::populationBalance::coalescenceModel::precompute()
{}


bool Foam::populationBalance::coalescenceModel::coalesces() const
{
    return true;
}


// ************************************************************************* //
