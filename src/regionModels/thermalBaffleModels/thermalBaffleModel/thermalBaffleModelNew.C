/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "thermalBaffleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace thermalBaffleModels
{

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<thermalBaffleModel> thermalBaffleModel::New(const fvMesh& mesh)
{
    word modelType;
    {
        IOdictionary thermalBafflePropertiesDict
        (
            IOobject
            (
                "thermalBaffleProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        );

        word modelType =
            thermalBafflePropertiesDict.lookupOrDefault<word>
            (
                "thermalBaffleModel",
                "thermalBaffle"
            );
    }

    meshConstructorTable::iterator cstrIter =
        meshConstructorTablePtr_->find(modelType);

    if (cstrIter == meshConstructorTablePtr_->end())
    {

        FatalErrorIn("thermalBaffleModel::New(const fvMesh&)")
            << "Unknown thermalBaffleModel type " << modelType
            << nl << nl
            <<  "Valid thermalBaffleModel types are:" << nl
            << meshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<thermalBaffleModel>(cstrIter()(modelType, mesh));
}


autoPtr<thermalBaffleModel> thermalBaffleModel::New
(
    const fvMesh& mesh,
    const dictionary& dict
)
{
    word modelType =
        dict.lookupOrDefault<word>("thermalBaffleModel", "thermalBaffle");

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {

        FatalErrorIn
        (
            "thermalBaffleModel::New(const fvMesh&, const dictionary&)"
        )   << "Unknown thermalBaffleModel type " << modelType
            << nl << nl
            <<  "Valid thermalBaffleModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<thermalBaffleModel>(cstrIter()(modelType, mesh, dict));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalBaffleModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
