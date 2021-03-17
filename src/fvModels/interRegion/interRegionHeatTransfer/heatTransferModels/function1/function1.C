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

#include "function1.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(function1, 0);
    addToRunTimeSelectionTable(heatTransferModel, function1, mesh);
    addToRunTimeSelectionTable(heatTransferModel, function1, model);
}
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferModels::function1::readCoeffs()
{
    UName_ = coeffs().lookupOrDefault<word>("U", "U");

    htcFunc_.reset(Function1<scalar>::New("htcFunc", coeffs()).ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferModels::function1::function1
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatTransferModel(typeName, dict, mesh),
    UName_(word::null),
    htcFunc_(nullptr)
{
    readCoeffs();
}


Foam::fv::heatTransferModels::function1::function1
(
    const dictionary& dict,
    const interRegionModel& model
)
:
    function1(dict, model.mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferModels::function1::~function1()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::fv::heatTransferModels::function1::htc() const
{
    const volVectorField& U = mesh().lookupObject<volVectorField>(UName_);

    tmp<volScalarField> tHtc =
        volScalarField::New
        (
            type() + ":htc",
            mesh(),
            dimPower/dimTemperature/dimArea,
            zeroGradientFvPatchScalarField::typeName
        );

    tHtc->primitiveFieldRef() = htcFunc_->value(mag(U));
    tHtc->correctBoundaryConditions();

    return tHtc;
}


void Foam::fv::heatTransferModels::function1::correct()
{}


bool Foam::fv::heatTransferModels::function1::read(const dictionary& dict)
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
