/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "PDRkEpsilon.H"
#include "PDRDragModel.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(PDRkEpsilon, 0);
addToRunTimeSelectionTable(RASModel, PDRkEpsilon, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

PDRkEpsilon::PDRkEpsilon
(
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi,
    const fluidThermo& thermophysicalModel,
    const word& turbulenceModelName,
    const word& modelName
)
:
    kEpsilon(rho, U, phi, thermophysicalModel, turbulenceModelName, modelName),

    C4_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C4",
            coeffDict_,
            0.1
        )
    )
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

PDRkEpsilon::~PDRkEpsilon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool PDRkEpsilon::read()
{
    if (RASModel::read())
    {
        C4_.readIfPresent(coeffDict_);
        return true;
    }
    else
    {
        return false;
    }
}


void PDRkEpsilon::correct()
{
    if (!turbulence_)
    {
        // Re-calculate viscosity
        mut_ = rho_*Cmu_*sqr(k_)/epsilon_;
        mut_.correctBoundaryConditions();

        // Re-calculate thermal diffusivity
        alphat_ = mut_/Prt_;
        alphat_.correctBoundaryConditions();

        return;
    }

    RASModel::correct();

    volScalarField divU(fvc::div(phi_/fvc::interpolate(rho_)));

    if (mesh_.moving())
    {
        divU += fvc::div(mesh_.phi());
    }

    tmp<volTensorField> tgradU = fvc::grad(U_);
    volScalarField G(GName(), mut_*(tgradU() && dev(twoSymm(tgradU()))));
    tgradU.clear();

    // Update espsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Add the blockage generation term so that it is included consistently
    // in both the k and epsilon equations
    const volScalarField& betav =
        U_.db().lookupObject<volScalarField>("betav");

    const volScalarField& Lobs =
        U_.db().lookupObject<volScalarField>("Lobs");

    const PDRDragModel& drag =
        U_.db().lookupObject<PDRDragModel>("PDRDragModel");

    volScalarField GR(drag.Gk());

    volScalarField LI
        (C4_*(Lobs + dimensionedScalar("minLength", dimLength, VSMALL)));

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        betav*fvm::ddt(rho_, epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1_*betav*G*epsilon_/k_
      + 1.5*pow(Cmu_, 3.0/4.0)*GR*sqrt(k_)/LI
      - fvm::SuSp(((2.0/3.0)*C1_)*betav*rho_*divU, epsilon_)
      - fvm::Sp(C2_*betav*rho_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        betav*fvm::ddt(rho_, k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        betav*G + GR
      - fvm::SuSp((2.0/3.0)*betav*rho_*divU, k_)
      - fvm::Sp(betav*rho_*epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);

    // Re-calculate viscosity
    mut_ = rho_*Cmu_*sqr(k_)/epsilon_;
    mut_.correctBoundaryConditions();

    // Re-calculate thermal diffusivity
    alphat_ = mut_/Prt_;
    alphat_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
