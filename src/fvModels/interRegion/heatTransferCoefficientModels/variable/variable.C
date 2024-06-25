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

#include "variable.H"
#include "fluidThermophysicalTransportModel.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace heatTransferCoefficientModels
{
    defineTypeNameAndDebug(variable, 0);
    addToRunTimeSelectionTable(heatTransferCoefficientModel, variable, mesh);
    addToRunTimeSelectionTable(heatTransferCoefficientModel, variable, model);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferCoefficientModels::variable::readCoeffs
(
    const dictionary& dict
)
{
    UName_ = dict.lookupOrDefault<word>("U", "U");

    a_ = dict.lookup<scalar>("a");
    b_ = dict.lookup<scalar>("b");
    c_ = dict.lookup<scalar>("c");
    L_ = dimensionedScalar("L", dimLength, dict);
    Pr_ = dimensionedScalar("Pr", dimless, dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferCoefficientModels::variable::variable
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatTransferCoefficientModel(typeName, dict, mesh),
    UName_(word::null),
    a_(NaN),
    b_(NaN),
    c_(NaN),
    L_("L", dimLength, NaN),
    Pr_("Pr", dimless, NaN),
    htc_
    (
        IOobject
        (
            typedName("htc"),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower/dimTemperature/dimArea, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    readCoeffs(dict);
}


Foam::fv::heatTransferCoefficientModels::variable::variable
(
    const dictionary& dict,
    const interRegionModel& model
)
:
    variable(dict, model.mesh())
{
    readCoeffs(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferCoefficientModels::variable::~variable()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::heatTransferCoefficientModels::variable::correct()
{
    const fluidThermophysicalTransportModel& ttm =
        mesh_.lookupType<fluidThermophysicalTransportModel>();

    const compressibleMomentumTransportModel& mtm =
        ttm.momentumTransport();

    const volVectorField& U =
        mesh_.lookupObject<volVectorField>(UName_);

    const volScalarField Re(mag(U)*L_/mtm.nuEff());
    const volScalarField Nu(a_*pow(Re, b_)*pow(Pr_, c_));

    htc_ = Nu*ttm.kappaEff()/L_;
    htc_.correctBoundaryConditions();
}


bool Foam::fv::heatTransferCoefficientModels::variable::read
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
