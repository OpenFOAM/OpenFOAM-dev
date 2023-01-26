/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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
namespace heatTransferCoefficientModels
{
    defineTypeNameAndDebug(constant, 0);
    addToRunTimeSelectionTable(heatTransferCoefficientModel, constant, mesh);
    addToRunTimeSelectionTable(heatTransferCoefficientModel, constant, model);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferCoefficientModels::constant::readCoeffs
(
    const dictionary& dict
)
{
    typeIOobject<volScalarField> htcIO
    (
        "htc",
        mesh_.time().constant(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (dict.found("htc"))
    {
        htc_ =
            dimensionedScalar
            (
                "htc",
                dimPower/dimTemperature/dimArea,
                dict
            );
        htcPtr_.clear();
    }
    else if (htcIO.headerOk())
    {
        htc_ = dimensionedScalar("htc", dimPower/dimTemperature/dimArea, NaN);
        htcPtr_.set(new volScalarField(htcIO, mesh_));
    }
    else
    {
        FatalIOErrorInFunction(dict)
            << "Heat transfer coefficient (htc) not found. A uniform htc "
            << "value should be specified, or a non-uniform field should "
            << "exist at " << htcIO.objectPath()
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferCoefficientModels::constant::constant
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatTransferCoefficientModel(typeName, dict, mesh),
    htc_("htc", dimPower/dimTemperature/dimArea, NaN),
    htcPtr_(nullptr)
{
    readCoeffs(dict);
}


Foam::fv::heatTransferCoefficientModels::constant::constant
(
    const dictionary& dict,
    const interRegionModel& model
)
:
    constant(dict, model.mesh())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferCoefficientModels::constant::~constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::fv::heatTransferCoefficientModels::constant::htc() const
{
    if (!htcPtr_.valid())
    {
        return volScalarField::New(typedName("htc"), mesh_, htc_);
    }
    else
    {
        return htcPtr_();
    }
}


void Foam::fv::heatTransferCoefficientModels::constant::correct()
{}


bool Foam::fv::heatTransferCoefficientModels::constant::read
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
