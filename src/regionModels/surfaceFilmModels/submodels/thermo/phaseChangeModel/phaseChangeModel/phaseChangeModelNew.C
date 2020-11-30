/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<phaseChangeModel> phaseChangeModel::New
(
    surfaceFilmRegionModel& model,
    const dictionary& dict
)
{
    if
    (
        !dict.found("phaseChangeModel")
     && !dict.found(phaseChangeModel::typeName)
    )
    {
        return autoPtr<phaseChangeModel>(new noPhaseChange(model, dict));
    }

    const dictionary& phaseChangeDict
    (
        dict.found("phaseChangeModel")
      ? dict
      : dict.subDict(phaseChangeModel::typeName)
    );

    const word modelType
    (
        dict.found("phaseChangeModel")
      ? phaseChangeDict.lookup("phaseChangeModel")
      : phaseChangeDict.lookup("model")
    );

    Info<< "    Selecting phaseChangeModel " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown phaseChangeModel type " << modelType
            << nl << nl << "Valid phaseChangeModel types are:" << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<phaseChangeModel>(cstrIter()(model, phaseChangeDict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // end namespace surfaceFilmModels
} // end namespace regionModels
} // end namespace Foam

// ************************************************************************* //
