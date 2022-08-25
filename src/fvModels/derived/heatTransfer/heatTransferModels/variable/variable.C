/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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
#include "thermophysicalTransportModel.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
namespace heatTransferModels
{
    defineTypeNameAndDebug(variable, 0);
    addToRunTimeSelectionTable(heatTransferModel, variable, mesh);
    addToRunTimeSelectionTable(heatTransferModel, variable, model);
}
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::heatTransferModels::variable::readCoeffs()
{
    UName_ = coeffs().lookupOrDefault<word>("U", "U");

    a_ = coeffs().lookup<scalar>("a");
    b_ = coeffs().lookup<scalar>("b");
    c_ = coeffs().lookup<scalar>("c");
    L_ = dimensionedScalar("L", dimLength, coeffs());
    Pr_ = dimensionedScalar("Pr", dimless, coeffs());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::heatTransferModels::variable::variable
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    heatTransferModel(typeName, dict, mesh),
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
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower/dimTemperature/dimArea, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    readCoeffs();
}


Foam::fv::heatTransferModels::variable::variable
(
    const dictionary& dict,
    const interRegionModel& model
)
:
    variable(dict, model.mesh())
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::heatTransferModels::variable::~variable()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::heatTransferModels::variable::correct()
{
    const thermophysicalTransportModel& ttm =
        mesh().lookupObject<thermophysicalTransportModel>
        (
            thermophysicalTransportModel::typeName
        );
    const compressibleMomentumTransportModel& mtm =
        ttm.momentumTransport();

    const volVectorField& U =
        mesh().lookupObject<volVectorField>(UName_);

    const volScalarField Re(mag(U)*L_/mtm.nuEff());
    const volScalarField Nu(a_*pow(Re, b_)*pow(Pr_, c_));

    htc_ = Nu*ttm.kappaEff()/L_;
    htc_.correctBoundaryConditions();
}


bool Foam::fv::heatTransferModels::variable::read(const dictionary& dict)
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
