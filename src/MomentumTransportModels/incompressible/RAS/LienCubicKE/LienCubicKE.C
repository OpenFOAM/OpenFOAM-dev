/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
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

#include "LienCubicKE.H"
#include "wallDist.H"
#include "bound.H"
#include "makeMomentumTransportModel.H"

makeMomentumTransportModelTypes
(
    geometricOneField,
    geometricOneField,
    incompressibleMomentumTransportModel
)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LienCubicKE, 0);
addToRunTimeSelectionTable
(
    RASincompressibleMomentumTransportModel,
    LienCubicKE,
    dictionary
);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void LienCubicKE::boundEpsilon()
{
    epsilon_ = max(epsilon_, Cmu_*sqr(k_)/(nutMaxCoeff_*nu()));
}


tmp<volScalarField> LienCubicKE::fMu() const
{
    const volScalarField yStar(sqrt(k_)*y()/nu());

    return
        (scalar(1) - exp(-Anu_*yStar))
       *(scalar(1) + (2*kappa_/(pow(Cmu_, 0.75))/(yStar + small)));
}


tmp<volScalarField> LienCubicKE::f2() const
{
    tmp<volScalarField> Rt = sqr(k_)/(nu()*epsilon_);

    return scalar(1) - 0.3*exp(-sqr(Rt));
}


tmp<volScalarField> LienCubicKE::E(const volScalarField& f2) const
{
    const volScalarField yStar(sqrt(k_)*y()/nu());
    const volScalarField le
    (
        kappa_*y()/(scalar(1) + (2*kappa_/(pow(Cmu_, 0.75))/(yStar + small)))
    );

    return
        (Ceps2_*pow(Cmu_, 0.75))
       *(f2*sqrt(k_)*epsilon_/le)*exp(-AE_*sqr(yStar));
}


void LienCubicKE::correctNut()
{
    correctNonlinearStress(fvc::grad(U_));
}


void LienCubicKE::correctNonlinearStress(const volTensorField& gradU)
{
    volSymmTensorField S(symm(gradU));
    volTensorField W(skew(gradU));

    volScalarField sBar((k_/epsilon_)*sqrt(2.0)*mag(S));
    volScalarField wBar((k_/epsilon_)*sqrt(2.0)*mag(W));

    volScalarField Cmu((2.0/3.0)/(Cmu1_ + sBar + Cmu2_*wBar));
    volScalarField fMu(this->fMu());

    boundEpsilon();
    nut_ = fMu*Cmu*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    nonlinearStress_ =
        fMu*k_
       *(
            // Quadratic terms
            sqr(k_/epsilon_)/(Cbeta_ + pow3(sBar))
           *(
                Cbeta1_*dev(innerSqr(S))
              + Cbeta2_*twoSymm(S&W)
              + Cbeta3_*dev(symm(W&W))
            )

            // Cubic terms
          - pow3(Cmu*k_/epsilon_)
           *(
                (Cgamma1_*magSqr(S) - Cgamma2_*magSqr(W))*S
              + Cgamma4_*twoSymm((innerSqr(S)&W))
            )
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LienCubicKE::LienCubicKE
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    nonlinearEddyViscosity<incompressible::RASModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    Ceps1_("Ceps1", coeffDict(), 1.44),
    Ceps2_("Ceps2", coeffDict(), 1.92),
    sigmak_("sigmak", coeffDict(), 1.0),
    sigmaEps_("sigmaEps", coeffDict(), 1.3),
    Cmu_("Cmu", coeffDict(), 0.09),
    Cmu1_("Cmu1", coeffDict(), 1.25),
    Cmu2_("Cmu2", coeffDict(), 0.9),
    Cbeta_("Cbeta", coeffDict(), 1000.0),
    Cbeta1_("Cbeta1", coeffDict(), 3.0),
    Cbeta2_("Cbeta2", coeffDict(), 15.0),
    Cbeta3_("Cbeta3", coeffDict(), -19.0),
    Cgamma1_("Cgamma1", coeffDict(), 16.0),
    Cgamma2_("Cgamma2", coeffDict(), 16.0),
    Cgamma4_("Cgamma4", coeffDict(), -80.0),
    kappa_("kappa", coeffDict(), 0.41),
    Anu_("Anu", coeffDict(), 0.0198),
    AE_("AE", coeffDict(), 0.00375),

    k_
    (
        IOobject
        (
            this->groupName("k"),
            runTime_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            this->groupName("epsilon"),
            runTime_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    bound(k_, kMin_);
    boundEpsilon();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool LienCubicKE::read()
{
    if (nonlinearEddyViscosity<incompressible::RASModel>::read())
    {
        Ceps1_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        Cmu_.readIfPresent(coeffDict());
        Cmu1_.readIfPresent(coeffDict());
        Cmu2_.readIfPresent(coeffDict());
        Cbeta_.readIfPresent(coeffDict());
        Cbeta1_.readIfPresent(coeffDict());
        Cbeta2_.readIfPresent(coeffDict());
        Cbeta3_.readIfPresent(coeffDict());
        Cgamma1_.readIfPresent(coeffDict());
        Cgamma2_.readIfPresent(coeffDict());
        Cgamma4_.readIfPresent(coeffDict());
        kappa_.readIfPresent(coeffDict());
        Anu_.readIfPresent(coeffDict());
        AE_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void LienCubicKE::correct()
{
    if (!turbulence_)
    {
        return;
    }

    nonlinearEddyViscosity<incompressible::RASModel>::correct();

    tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    volScalarField G
    (
        GName(),
        (nut_*twoSymm(gradU) - nonlinearStress_) && gradU
    );


    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    const volScalarField f2(this->f2());

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        Ceps1_*G*epsilon_/k_
      - fvm::Sp(Ceps2_*f2*epsilon_/k_, epsilon_)
      + E(f2)
    );

    epsEqn.ref().relax();
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    boundEpsilon();


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
      ==
        G
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn.ref().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity and non-linear stress
    correctNonlinearStress(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
