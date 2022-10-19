/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2022 OpenFOAM Foundation
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

#include "filmMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmSubModels
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<momentumTransportModel> momentumTransportModel::New
(
    surfaceFilm& model,
    const dictionary& dict
)
{
    dict.lookupEntryBackwardsCompatible
    (
        {momentumTransportModel::typeName, "turbulence"},
        false,
        true
    );

    const dictionary& momentumTransportDict
    (
        dict.found(momentumTransportModel::typeName)
      ? dict.subDict(momentumTransportModel::typeName)
      : dict
    );

    const word modelType
    (
        dict.found(momentumTransportModel::typeName)
      ? momentumTransportDict.lookup("model")
      : momentumTransportDict.lookup("turbulence")
    );

    Info<< "    Selecting momentumTransportModel " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown momentumTransportModel type " << modelType
            << nl << nl << "Valid momentumTransportModel types are:" << nl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<momentumTransportModel>
    (
        cstrIter()(model, momentumTransportDict)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmSubModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
