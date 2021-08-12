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

#include "constant.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(heatTransferModel, constant, mesh);
    addToRunTimeSelectionTable(heatTransferModel, constant, model);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferModels::constant::readCoeffs()
{
    typeIOobject<volScalarField> htcIO
    (
        "htc",
        mesh().time().constant(),
        mesh(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (coeffs().found("htc"))
    {
        htc_ =
            dimensionedScalar
            (
                "htc",
                dimPower/dimTemperature/dimArea,
                coeffs()
            );
        htcPtr_.clear();
    }
    else if (htcIO.headerOk())
    {
        htc_ = dimensionedScalar("htc", dimPower/dimTemperature/dimArea, NaN);
        htcPtr_.set(new volScalarField(htcIO, mesh()));
    }
    else
    {
        FatalIOErrorInFunction(coeffs())
            << "Heat transfer coefficient (htc) not found. A uniform htc "
            << "value should be specified, or a non-uniform field should "
            << "exist at " << htcIO.objectPath()
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferModels::constant::constant
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatTransferModel(typeName, dict, mesh),
    htc_("htc", dimPower/dimTemperature/dimArea, NaN),
    htcPtr_(nullptr)
{
    readCoeffs();
}


Foam::fv::heatTransferModels::constant::constant
(
    const dictionary& dict,
    const interRegionModel& model
)
:
    constant(dict, model.mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferModels::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::fv::heatTransferModels::constant::htc() const
{
    if (!htcPtr_.valid())
    {
        return volScalarField::New(type() + ":htc", mesh(), htc_);
    }
    else
    {
        return htcPtr_();
    }
}


void Foam::fv::heatTransferModels::constant::correct()
{}


bool Foam::fv::heatTransferModels::constant::read(const dictionary& dict)
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
