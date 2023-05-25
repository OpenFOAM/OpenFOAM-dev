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

#include "transport.H"
#include "surfaceInterpolate.H"
#include "fvmDdt.H"
#include "fvcLaplacian.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiModels
{
    defineTypeNameAndDebug(transport, 0);
    addToRunTimeSelectionTable(XiModel, transport, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiModels::transport::transport
(
    const dictionary& XiProperties,
    const psiuMulticomponentThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su,
    const volScalarField& rho,
    const volScalarField& b,
    const surfaceScalarField& phi
)
:
    XiModel(XiProperties, thermo, turbulence, Su, rho, b, phi),
    XiShapeCoef(XiModelCoeffs_.lookup<scalar>("XiShapeCoef")),
    XiEqModel_(XiEqModel::New(XiProperties, thermo, turbulence, Su)),
    XiGModel_(XiGModel::New(XiProperties, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiModels::transport::~transport()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiModels::transport::Db() const
{
    return XiGModel_->Db();
}


void Foam::XiModels::transport::correct
(
    const fv::convectionScheme<scalar>& mvConvection
)
{
    volScalarField XiEqEta(XiEqModel_->XiEq());
    volScalarField GEta(XiGModel_->G());

    volScalarField R(GEta*XiEqEta/(XiEqEta - 0.999));

    volScalarField XiEqStar(R/(R - GEta));

    volScalarField XiEq
    (
        1.0 + (1.0 + (2*XiShapeCoef)*(0.5 - b_))*(XiEqStar - 1.0)
    );

    volScalarField G(R*(XiEq - 1.0)/XiEq);

    const objectRegistry& db = b_.db();
    const volScalarField& betav = db.lookupObject<volScalarField>("betav");
    const volScalarField& mgb = db.lookupObject<volScalarField>("mgb");
    const surfaceScalarField& phiSt =
        db.lookupObject<surfaceScalarField>("phiSt");
    const volScalarField& Db = db.lookupObject<volScalarField>("Db");
    const surfaceScalarField& nf = db.lookupObject<surfaceScalarField>("nf");

    surfaceScalarField phiXi
    (
        "phiXi",
        phiSt
      + (
          - fvc::interpolate(fvc::laplacian(Db, b_)/mgb)*nf
          + fvc::interpolate(rho_)*fvc::interpolate(Su_*(1.0/Xi_ - Xi_))*nf
        )
    );

    solve
    (
        betav*fvm::ddt(rho_, Xi_)
      + mvConvection.fvmDiv(phi_, Xi_)
      + fvm::div(phiXi, Xi_)
      - fvm::Sp(fvc::div(phiXi), Xi_)
     ==
        betav*rho_*R
      - fvm::Sp(betav*rho_*(R - G), Xi_)
    );

    // Correct boundedness of Xi
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    Xi_.max(1.0);
    Xi_ = min(Xi_, 2.0*XiEq);
}


bool Foam::XiModels::transport::read(const dictionary& XiProperties)
{
    XiModel::read(XiProperties);

    XiModelCoeffs_.lookup("XiShapeCoef") >> XiShapeCoef;

    return true;
}


// ************************************************************************* //
