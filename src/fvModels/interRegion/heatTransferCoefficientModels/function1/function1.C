/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "function1.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace heatTransferCoefficientModels
{
    defineTypeNameAndDebug(function1, 0);
    addToRunTimeSelectionTable(heatTransferCoefficientModel, function1, mesh);
    addToRunTimeSelectionTable(heatTransferCoefficientModel, function1, model);
}
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferCoefficientModels::function1::readCoeffs
(
    const dictionary& dict
)
{
    UName_ = dict.lookupOrDefault<word>("U", "U");

    htcFunc_.reset
    (
        Function1<scalar>::New
        (
            "htcFunc",
            dimVelocity,
            dimPower/dimArea/dimTemperature,
            dict
        ).ptr()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferCoefficientModels::function1::function1
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatTransferCoefficientModel(typeName, dict, mesh),
    UName_(word::null),
    htcFunc_(nullptr)
{
    readCoeffs(dict);
}


Foam::fv::heatTransferCoefficientModels::function1::function1
(
    const dictionary& dict,
    const interRegionModel& model
)
:
    function1(dict, model.mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferCoefficientModels::function1::~function1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::fv::heatTransferCoefficientModels::function1::htc() const
{
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    tmp<volScalarField> tHtc =
        volScalarField::New
        (
            typedName("htc"),
            mesh_,
            dimPower/dimTemperature/dimArea,
            zeroGradientFvPatchScalarField::typeName
        );

    tHtc->primitiveFieldRef() = htcFunc_->value(mag(U));
    tHtc->correctBoundaryConditions();

    return tHtc;
}


void Foam::fv::heatTransferCoefficientModels::function1::correct()
{}


bool Foam::fv::heatTransferCoefficientModels::function1::read
(
    const dictionary& dict
)
{
    if (heatTransferCoefficientModel::read(dict))
    {
        readCoeffs(dict);
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
