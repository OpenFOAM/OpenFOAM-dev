/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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

#include "heatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(heatTransferCoefficientModel, 0);
    defineRunTimeSelectionTable(heatTransferCoefficientModel, mesh);
    defineRunTimeSelectionTable(heatTransferCoefficientModel, model);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferCoefficientModel::heatTransferCoefficientModel
(
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh)
{}


Foam::fv::heatTransferCoefficientModel::heatTransferCoefficientModel
(
    const word& modelType,
    const dictionary& dict,
    const interRegionModel& model
)
:
    heatTransferCoefficientModel(modelType, dict, model.mesh())
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fv::heatTransferCoefficientModel>
Foam::fv::heatTransferCoefficientModel::New
(
    const dictionary& dict,
    const fvMesh& mesh
)
{
    word heatTransferCoefficientModelType(dict.lookup(typeName));

    Info<< "Selecting " << typeName << " "
        << heatTransferCoefficientModelType << endl;

    meshConstructorTable::iterator cstrIter =
        meshConstructorTablePtr_->find(heatTransferCoefficientModelType);

    if (cstrIter == meshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type "
            << heatTransferCoefficientModelType << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << meshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return
        cstrIter()
        (
            dict.optionalSubDict(heatTransferCoefficientModelType + "Coeffs"),
            mesh
        );
}


Foam::autoPtr<Foam::fv::heatTransferCoefficientModel>
Foam::fv::heatTransferCoefficientModel::New
(
    const dictionary& dict,
    const interRegionModel& model
)
{
    word heatTransferCoefficientModelType(dict.lookup(typeName));

    Info<< "Selecting " << typeName << " "
        << heatTransferCoefficientModelType << endl;

    modelConstructorTable::iterator cstrIter =
        modelConstructorTablePtr_->find(heatTransferCoefficientModelType);

    if (cstrIter == modelConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type "
            << heatTransferCoefficientModelType << endl << endl
            << "Valid " << typeName << " types are : " << endl
            << modelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return
        cstrIter()
        (
            dict.optionalSubDict(heatTransferCoefficientModelType + "Coeffs"),
            model
        );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferCoefficientModel::~heatTransferCoefficientModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::fv::heatTransferCoefficientModel::read(const dictionary& dict)
{
    return true;
}


// ************************************************************************* //
