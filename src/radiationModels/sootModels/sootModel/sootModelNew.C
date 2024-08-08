/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2024 OpenFOAM Foundation
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
#include "noSoot.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::radiationModels::sootModel>
Foam::radiationModels::sootModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    // Get the soot model type name
    word modelType(sootModels::noSoot::typeName);
    if (dict.found(sootModel::typeName))
    {
        dict.lookup(sootModel::typeName) >> modelType;
        Info<< "Selecting soot model " << modelType << endl;
    }

    // Backwards compatibility check
    const wordList cmpts(basicThermo::splitThermoName(modelType, 3));
    if (cmpts.size() == 3)
    {
        modelType = cmpts[0];

        WarningInFunction
            << "Template parameters are no longer required when selecting a "
            << sootModel::typeName << ". This information is now "
            << "obtained directly from the thermodynamics. Actually selecting "
            << "combustion model " << modelType << "." << endl;
    }

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorInFunction(dict)
            << "Unknown " << sootModel::typeName << " type "
            << modelType << nl << nl
            << "Valid " << sootModel::typeName << " types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<sootModel>
    (
        cstrIter()(dict.optionalSubDict(modelType + "Coeffs"), mesh, modelType)
    );
}


// ************************************************************************* //
