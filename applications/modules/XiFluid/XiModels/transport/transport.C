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

#include "transport.H"
#include "XiEqModel.H"
#include "XiProfile.H"
#include "XiGModel.H"
#include "fvmDiv.H"
#include "fvcLaplacian.H"
#include "fvModels.H"
#include "fvConstraints.H"
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


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::XiModels::transport::readCoeffs(const dictionary& dict)
{
    XiModel::readCoeffs(dict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiModels::transport::transport
(
    const dictionary& dict,
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& turbulence,
    const volScalarField& Su
)
:
    XiModel(thermo, turbulence, Su),
    XiEqModel_(XiEqModel::New(dict, thermo, turbulence, Su)),
    XiProfile_(XiProfile::New(dict, b_)),
    XiGModel_(XiGModel::New(dict, thermo, turbulence, Su))
{
    readCoeffs(dict);
    Xi_.writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::XiModels::transport::~transport()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiModels::transport::Db() const
{
    return XiEqModel_->Db();
}


void Foam::XiModels::transport::correct()
{
    const fvMesh& mesh(thermo_.mesh());

    const Foam::fvModels& fvModels(Foam::fvModels::New(mesh));
    const Foam::fvConstraints& fvConstraints(Foam::fvConstraints::New(mesh));

    const volScalarField XiEqEta(XiEqModel_->XiEq());
    const volScalarField GEta(XiGModel_->G());

    const volScalarField R(GEta*XiEqEta/(XiEqEta - 0.999));

    const volScalarField XiEqStar(R/(R - GEta));

    const volScalarField XiEq(1 + XiProfile_->profile()*(XiEqStar - 1));

    const volScalarField G(R*(XiEq - 1)/XiEq);

    const volScalarField& mgb = mesh.lookupObject<volScalarField>("mgb");
    const surfaceScalarField& phiSt =
        mesh.lookupObject<surfaceScalarField>("phiSt");
    const volScalarField& Db = mesh.lookupObject<volScalarField>("Db");
    const volVectorField& n = mesh.lookupObject<volVectorField>("n");
    const surfaceScalarField& nf = mesh.lookupObject<surfaceScalarField>("nf");

    const surfaceScalarField phiXi
    (
        "phiXi",
        phiSt
      + (
          - fvc::interpolate(fvc::laplacian(Db, b_)/mgb)*nf
          + fvc::interpolate(rho_)*fvc::interpolate(Su_*(1/Xi_ - Xi_))*nf
        )
    );

    const surfaceScalarField& phi = turbulence_.alphaRhoPhi();

    const volVectorField& U(turbulence_.U());

    const volVectorField Ut(U + Su_*Xi_*n);
    const volScalarField sigmat((n & n)*fvc::div(Ut) - (n & fvc::grad(Ut) & n));

    const volScalarField sigmas
    (
        ((n & n)*fvc::div(U) - (n & fvc::grad(U) & n))/Xi_
      + (
            (n & n)*fvc::div(Su_*n)
          - (n & fvc::grad(Su_*n) & n)
        )*(Xi_ + scalar(1))/(2*Xi_)
    );

    fvScalarMatrix XiEqn
    (
        fvm::ddt(rho_, Xi_)
      + fvm::div(phi + phiXi, Xi_, "div(phiXi,Xi)")
      - fvm::Sp(fvc::div(phiXi), Xi_)
     ==
        rho_*R
      - fvm::Sp(rho_*(R - G), Xi_)
      - fvm::Sp
        (
            rho_*max
            (
                sigmat - sigmas,
                dimensionedScalar(sigmat.dimensions(), 0)
            ),
            Xi_
        )
      + fvModels.source(rho_, Xi_)
    );

    XiEqn.relax();

    fvConstraints.constrain(XiEqn);

    XiEqn.solve();

    fvConstraints.constrain(Xi_);

    // Limit range of Xi for realisability and stability
    Xi_.max(1);
    Xi_ = min(Xi_, 2*XiEq);
}


// ************************************************************************* //
