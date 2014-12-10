/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "relaxationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationModel, 0);
defineRunTimeSelectionTable(relaxationModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

relaxationModel::relaxationModel
(
    const word& type,
    const dictionary& relaxationDict,
    const Time& runTime
)
:
    dictionary(relaxationDict),
    runTime_(runTime),
    coeffDict_(subDict(type + "Coeffs"))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<relaxationModel> relaxationModel::New
(
    const dictionary& relaxationDict,
    const Time& runTime
)
{
    word relaxationModelTypeName
    (
        relaxationDict.lookup("relaxationModel")
    );

    Info<< nl << "Selecting relaxationModel "
        << relaxationModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(relaxationModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "relaxationModel::New(const dictionary&, "
            "const conformalVoronoiMesh&)"
        )   << "Unknown relaxationModel type "
            << relaxationModelTypeName
            << endl << endl
            << "Valid relaxationModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<relaxationModel>(cstrIter()(relaxationDict, runTime));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

relaxationModel::~relaxationModel()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
