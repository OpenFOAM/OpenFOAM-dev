/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "LienCubicKELowRe.H"
#include "wallDist.H"
#include "wallFvPatch.H"
#include "bound.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LienCubicKELowRe, 0);
addToRunTimeSelectionTable(RASModel, LienCubicKELowRe, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> LienCubicKELowRe::fMu()
{
    return
        (scalar(1) - exp(-Am_*yStar_))
       /(scalar(1) - exp(-Aepsilon_*yStar_) + SMALL);
}


void LienCubicKELowRe::correctNut()
{
    nut_ =
        Cmu_*fMu()*sqr(k_)/epsilon_
        // C5 term, implicit
      + max
        (
            C5viscosity_,
            dimensionedScalar("0", C5viscosity_.dimensions(), 0.0)
        );

    nut_.correctBoundaryConditions();
}


void LienCubicKELowRe::correctNonlinearStress(const volTensorField& gradU)
{
    nonlinearStress_ = symm
    (
        // quadratic terms
        pow3(k_)/sqr(epsilon_)
       *(
            Ctau1_/fEta_
           *(
                (gradU & gradU)
              + (gradU & gradU)().T()
            )
          + Ctau2_/fEta_*(gradU & gradU.T())
          + Ctau3_/fEta_*(gradU.T() & gradU)
        )
        // cubic term C4
      - 20.0*pow4(k_)/pow3(epsilon_)
       *pow3(Cmu_)
       *(
            ((gradU & gradU) & gradU.T())
          + ((gradU & gradU.T()) & gradU.T())
          - ((gradU.T() & gradU) & gradU)
          - ((gradU.T() & gradU.T()) & gradU)
        )
        // cubic term C5, explicit part
      + min
        (
            C5viscosity_,
            dimensionedScalar("0", C5viscosity_.dimensions(), 0.0)
        )*gradU
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LienCubicKELowRe::LienCubicKELowRe
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
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
        transport,
        propertiesName
    ),

    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),
    A1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A1",
            coeffDict_,
            1.25
        )
    ),
    A2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A2",
            coeffDict_,
            1000.0
        )
    ),
    Ctau1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ctau1",
            coeffDict_,
            -4.0
        )
    ),
    Ctau2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ctau2",
            coeffDict_,
            13.0
        )
    ),
    Ctau3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ctau3",
            coeffDict_,
            -2.0
        )
    ),
    alphaKsi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaKsi",
            coeffDict_,
            0.9
        )
    ),
    CmuWall_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    Am_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Am",
            coeffDict_,
            0.016
        )
    ),
    Aepsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Aepsilon",
            coeffDict_,
            0.263
        )
    ),
    Amu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Amu",
            coeffDict_,
            0.00222
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            runTime_.timeName(),
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
            IOobject::groupName("epsilon", U.group()),
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    y_(wallDist::New(mesh_).y()),

    eta_
    (
        k_/bound(epsilon_, epsilonMin_)
       *sqrt(2.0*magSqr(0.5*(fvc::grad(U) + T(fvc::grad(U)))))
    ),
    ksi_
    (
        k_/epsilon_
       *sqrt(2.0*magSqr(0.5*(fvc::grad(U) - T(fvc::grad(U)))))
    ),
    Cmu_(2.0/(3.0*(A1_ + eta_ + alphaKsi_*ksi_))),
    fEta_(A2_ + pow3(eta_)),

    C5viscosity_
    (
       -2.0*pow3(Cmu_)*pow4(k_)/pow3(epsilon_)
       *(
            magSqr(fvc::grad(U) + T(fvc::grad(U)))
          - magSqr(fvc::grad(U) - T(fvc::grad(U)))
        )
    ),

    yStar_(sqrt(k_)*y_/nu() + SMALL)
{
    bound(k_, kMin_);
    // already bounded: bound(epsilon_, epsilonMin_);

    if (type == typeName)
    {
        correctNut();
        correctNonlinearStress(fvc::grad(U));
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool LienCubicKELowRe::read()
{
    if (nonlinearEddyViscosity<incompressible::RASModel>::read())
    {
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        A1_.readIfPresent(coeffDict());
        A2_.readIfPresent(coeffDict());
        Ctau1_.readIfPresent(coeffDict());
        Ctau2_.readIfPresent(coeffDict());
        Ctau3_.readIfPresent(coeffDict());
        alphaKsi_.readIfPresent(coeffDict());
        CmuWall_.readIfPresent(coeffDict());
        kappa_.readIfPresent(coeffDict());
        Am_.readIfPresent(coeffDict());
        Aepsilon_.readIfPresent(coeffDict());
        Amu_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void LienCubicKELowRe::correct()
{
    nonlinearEddyViscosity<incompressible::RASModel>::correct();

    if (!turbulence_)
    {
        return;
    }

    tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    yStar_ = sqrt(k_)*y_/nu() + SMALL;
    tmp<volScalarField> Rt = sqr(k_)/(nu()*epsilon_);

    const volScalarField f2
    (
        scalar(1) - 0.3*exp(-sqr(Rt))
    );

    volScalarField G
    (
        GName(),
        (nut_*twoSymm(gradU) - nonlinearStress_) && gradU
    );

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1_*G*epsilon_/k_
        // E-term
      + C2_*f2*pow(Cmu_, 0.75)*pow(k_, scalar(0.5))
       /(kappa_*y_*(scalar(1) - exp(-Aepsilon_*yStar_)))
       *exp(-Amu_*sqr(yStar_))*epsilon_
      - fvm::Sp(C2_*f2*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();
    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


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

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity

    eta_ = k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU + gradU.T())));
    ksi_ = k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU - gradU.T())));
    Cmu_ = 2.0/(3.0*(A1_ + eta_ + alphaKsi_*ksi_));
    fEta_ = A2_ + pow3(eta_);

    C5viscosity_ =
       -2.0*pow3(Cmu_)*pow4(k_)/pow3(epsilon_)
       *(magSqr(gradU + gradU.T()) - magSqr(gradU - gradU.T()));

    correctNut();
    correctNonlinearStress(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
