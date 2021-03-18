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

#include "function2.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(function2, 0);
    addToRunTimeSelectionTable(heatTransferModel, function2, model);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferModels::function2::readCoeffs()
{
    UName_ = coeffs().lookupOrDefault<word>("U", "U");
    UNbrName_ = coeffs().lookupOrDefault<word>("UNbr", "U");

    htcFunc_.reset(Function2<scalar>::New("htcFunc", coeffs()).ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferModels::function2::function2
(
    const dictionary& dict,
    const interRegionModel& model
)
:
    heatTransferModel(typeName, dict, model),
    model_(model),
    UName_(word::null),
    UNbrName_(word::null),
    htcFunc_(),
    htc_
    (
        IOobject
        (
            type() + ":htc",
            model.mesh().time().timeName(),
            model.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        model.mesh(),
        dimensionedScalar(dimPower/dimTemperature/dimArea, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferModels::function2::~function2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::heatTransferModels::function2::correct()
{
    const fvMesh& mesh = model_.mesh();
    const fvMesh& nbrMesh = model_.nbrMesh();

    const volVectorField& U = mesh.lookupObject<volVectorField>(UName_);

    const volVectorField& UNbr =
        nbrMesh.lookupObject<volVectorField>(UNbrName_);
    const scalarField UMagNbr(mag(UNbr));
    const scalarField UMagNbrMapped(model_.interpolate(UMagNbr));

    htc_.primitiveFieldRef() = htcFunc_->value(mag(U()), UMagNbrMapped);
    htc_.correctBoundaryConditions();
}


bool Foam::fv::heatTransferModels::function2::read(const dictionary& dict)
{
    if (heatTransferModel::read(dict))
    {
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
