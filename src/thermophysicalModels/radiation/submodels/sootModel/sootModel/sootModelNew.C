/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "error.H"
#include "sootModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::radiation::sootModel>
Foam::radiation::sootModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    word modelType("none");

    if (dict.found("sootModel"))
    {
        dict.lookup("sootModel") >> modelType;

        Info<< "Selecting sootModel " << modelType << endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "sootModel::New(const dictionary&, const fvMesh&)"
        )   << "Unknown sootModel type "
            << modelType << nl << nl
            << "Valid sootModel types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << exit(FatalError);
    }

    const label tempOpen = modelType.find('<');

    const word className = modelType(0, tempOpen);

    return autoPtr<sootModel>(cstrIter()(dict, mesh, className));
}


// ************************************************************************* //
