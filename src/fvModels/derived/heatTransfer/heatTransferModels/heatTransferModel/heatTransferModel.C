/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
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

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(heatTransferModel, 0);
    defineRunTimeSelectionTable(heatTransferModel, mesh);
    defineRunTimeSelectionTable(heatTransferModel, model);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferModel::readCoeffs()
{
    typeIOobject<volScalarField> AoVIO
    (
        "AoV",
        mesh().time().constant(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (coeffs().found("AoV"))
    {
        AoV_ = dimensionedScalar("AoV", dimless/dimLength, coeffs());
        AoVPtr_.clear();
    }
    else if (AoVIO.headerOk())
    {
        AoV_ = dimensionedScalar("AoV", dimless/dimLength, NaN);
        AoVPtr_.set(new volScalarField(AoVIO, mesh()));
    }
    else
    {
        FatalIOErrorInFunction(coeffs())
            << "Area per unit volume (AoV) not found. A uniform AoV "
            << "value should be specified, or a non-uniform field should "
            << "exist at " << AoVIO.objectPath()
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferModel::heatTransferModel
(
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    coeffs_(dict.optionalSubDict(modelType + "Coeffs")),
    AoV_("AoV", dimless/dimLength, NaN),
    AoVPtr_(nullptr)
{
    readCoeffs();
}


Foam::fv::heatTransferModel::heatTransferModel
(
    const word& modelType,
    const dictionary& dict,
    const interRegionModel& model
)
:
    heatTransferModel(modelType, dict, model.mesh())
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fv::heatTransferModel>
Foam::fv::heatTransferModel::New(const dictionary& dict, const fvMesh& mesh)
{
    word heatTransferModelType(dict.lookup("type"));

    Info<< "Selecting heatTransferModel "
        << heatTransferModelType << endl;

    meshConstructorTable::iterator cstrIter =
        meshConstructorTablePtr_->find(heatTransferModelType);

    if (cstrIter == meshConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown heatTransferModelType type "
            << heatTransferModelType << endl << endl
            << "Valid heatTransferModel types are : " << endl
            << meshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, mesh);
}


Foam::autoPtr<Foam::fv::heatTransferModel>
Foam::fv::heatTransferModel::New
(
    const dictionary& dict,
    const interRegionModel& model
)
{
    word heatTransferModelType(dict.lookup("type"));

    Info<< "Selecting heatTransferModel "
        << heatTransferModelType << endl;

    modelConstructorTable::iterator cstrIter =
        modelConstructorTablePtr_->find(heatTransferModelType);

    if (cstrIter == modelConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown heatTransferModelType type "
            << heatTransferModelType << endl << endl
            << "Valid heatTransferModel types are : " << endl
            << modelConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(dict, model);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferModel::~heatTransferModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::fv::heatTransferModel::AoV() const
{
    if (!AoVPtr_.valid())
    {
        return volScalarField::New(typedName("AoV"), mesh_, AoV_);
    }
    else
    {
        return AoVPtr_();
    }
}


bool Foam::fv::heatTransferModel::read(const dictionary& dict)
{
    coeffs_ = dict.optionalSubDict(type() + "Coeffs");

    readCoeffs();

    return true;
}


// ************************************************************************* //
