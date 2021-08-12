/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "surfaceFilmModel.H"
#include "noFilm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<surfaceFilmModel> surfaceFilmModel::New
(
    const fvMesh& mesh,
    const dimensionedVector& g,
    const word& regionType
)
{
    word modelType;

    {
        typeIOobject<IOdictionary> surfaceFilmPropertiesDictHeader
        (
            regionType + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (surfaceFilmPropertiesDictHeader.headerOk())
        {
            IOdictionary surfaceFilmPropertiesDict
            (
                surfaceFilmPropertiesDictHeader
            );

            surfaceFilmPropertiesDict.lookup("surfaceFilmModel") >> modelType;
        }
        else
        {
            modelType = surfaceFilmModels::noFilm::typeName;
        }
    }

    Info<< "Selecting surfaceFilmModel " << modelType << endl;

    meshConstructorTable::iterator cstrIter =
        meshConstructorTablePtr_->find(modelType);

    if (cstrIter == meshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown surfaceFilmModel type " << modelType
            << nl << nl << "Valid surfaceFilmModel types are:" << nl
            << meshConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<surfaceFilmModel>
    (
        cstrIter()
        (
            modelType,
            mesh,
            g,
            regionType
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
