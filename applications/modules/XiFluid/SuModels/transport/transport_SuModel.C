/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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

#include "transport_SuModel.H"
#include "laminarFlameSpeed.H"
#include "fvmDiv.H"
#include "fvcLaplacian.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace SuModels
{
    defineTypeNameAndDebug(transport, 0);
    addToRunTimeSelectionTable(SuModel, transport, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::SuModels::transport::readCoeffs(const dictionary& dict)
{
    SuModel::readCoeffs(dict);

    sigmaExt_.read(dict);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SuModels::transport::transport
(
    const dictionary& dict,
    const psiuMulticomponentThermo& thermo,
    const fluidThermoThermophysicalTransportModel& turbulence
)
:
    unstrained(dict, thermo, turbulence),
    sigmaExt_("sigmaExt", dimless/dimTime, dict),
    SuMin_(0.01*Su_.average()),
    SuMax_(4*Su_.average())
{
    Su_.writeOpt() = IOobject::AUTO_WRITE;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SuModels::transport::~transport()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SuModels::transport::correct()
{
    const fvMesh& mesh(thermo_.mesh());

    const Foam::fvModels& fvModels(Foam::fvModels::New(mesh));
    const Foam::fvConstraints& fvConstraints(Foam::fvConstraints::New(mesh));

    const volScalarField& rho = turbulence_.rho();
    const volScalarField& b = thermo_.Y("b");
    const volScalarField& mgb = mesh.lookupObject<volScalarField>("mgb");
    const volScalarField& Xi = mesh.lookupObject<volScalarField>("Xi");
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
          - fvc::interpolate(fvc::laplacian(Db, b)/mgb)*nf
          + fvc::interpolate(rho)*fvc::interpolate(Su_*(1/Xi - Xi))*nf
        )
    );

    const surfaceScalarField& phi = turbulence_.alphaRhoPhi();
    const volVectorField& U(turbulence_.U());

    const volScalarField sigmas
    (
        ((n & n)*fvc::div(U) - (n & fvc::grad(U) & n))/Xi
      + (
            (n & n)*fvc::div(Su_*n)
          - (n & fvc::grad(Su_*n) & n)
        )*(Xi + scalar(1))/(2*Xi)
    );

    const volScalarField Su0(Su0_()());

    const volScalarField SuInf
    (
        Su0*max(scalar(1) - sigmas/sigmaExt_, scalar(0.01))
    );

    const volScalarField Rc
    (
        (sigmas*SuInf*(Su0 - SuInf) + sqr(SuMin_)*sigmaExt_)
        /(sqr(Su0 - SuInf) + sqr(SuMin_))
    );

    fvScalarMatrix SuEqn
    (
        fvm::ddt(rho, Su_)
      + fvm::div(phi + phiXi, Su_, "div(phiXi,Su)")
      - fvm::Sp(fvc::div(phiXi), Su_)
      ==
      - fvm::SuSp(-rho*Rc*Su0/Su_, Su_)
      - fvm::SuSp(rho*(sigmas + Rc), Su_)
      + fvModels.source(rho, Su_)
    );

    SuEqn.relax();

    fvConstraints.constrain(SuEqn);

    SuEqn.solve();

    fvConstraints.constrain(Su_);

    // Limit the maximum and minimum Su
    Su_.min(SuMax_);
    Su_.max(SuMin_);
}


// ************************************************************************* //
