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

#include "radiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::radiationModel> Foam::radiationModel::New
(
    const volScalarField& T
)
{
    const IOdictionary radiationProperties
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const word modelType(radiationProperties.lookup("radiationModel"));

    Info<< "Selecting radiationModel " << modelType << endl;

    TConstructorTable::iterator cstrIter =
        TConstructorTablePtr_->find(modelType);

    if (cstrIter == TConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(radiationProperties)
            << "Unknown radiationModel type "
            << modelType << nl << nl
            << "Valid radiationModel types are:" << nl
            << TConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<radiationModel>(cstrIter()(T));
}


Foam::autoPtr<Foam::radiationModel> Foam::radiationModel::New
(
    const dictionary& dict,
    const volScalarField& T
)
{
    const word modelType(dict.lookup("radiationModel"));

    Info<< "Selecting radiationModel " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown radiationModel type "
            << modelType << nl << nl
            << "Valid radiationModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<radiationModel>(cstrIter()(dict, T));
}


// ************************************************************************* //
