/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

#include "noPhaseChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::twoPhaseChangeModel> Foam::twoPhaseChangeModel::New
(
    const immiscibleIncompressibleTwoPhaseMixture& mixture
)
{
    typeIOobject<IOdictionary> twoPhaseChangeModelIO
    (
        IOobject
        (
            phaseChangePropertiesName,
            mixture.U().time().constant(),
            mixture.U().db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word modelType(twoPhaseChangeModels::noPhaseChange::typeName);

    if (twoPhaseChangeModelIO.headerOk())
    {
        IOdictionary(twoPhaseChangeModelIO).lookup
        (
            twoPhaseChangeModel::typeName
        ) >> modelType;
    }
    else
    {
        Info<< "No phase change: "
            << twoPhaseChangeModelIO.name()
            << " not found" << endl;
    }

    Info<< "Selecting phaseChange model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << twoPhaseChangeModel::typeName<< " type "
            << modelType << nl << nl
            << "Valid  twoPhaseChangeModels are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<twoPhaseChangeModel>(cstrIter()(mixture));
}


// ************************************************************************* //
