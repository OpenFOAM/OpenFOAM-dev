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

#include "NonlinearKEShih.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"
#include "nutkWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(NonlinearKEShih, 0);
addToRunTimeSelectionTable(RASModel, NonlinearKEShih, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

NonlinearKEShih::NonlinearKEShih
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

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

    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa_",
            coeffDict_,
            0.41
        )
    ),
    E_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "E",
            coeffDict_,
            9.8
        )
    ),

    k_
    (
        IOobject
        (
            "k",
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
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

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
    fEta_(A2_ + pow(eta_, 3.0)),

    nut_("nut", Cmu_*sqr(k_)/epsilon_),

    nonlinearStress_
    (
        "nonlinearStress",
        symm
        (
            pow3(k_)/sqr(epsilon_)
           *(
                Ctau1_/fEta_
               *(
                    (fvc::grad(U) & fvc::grad(U))
                  + (fvc::grad(U) & fvc::grad(U))().T()
                )
              + Ctau2_/fEta_*(fvc::grad(U) & T(fvc::grad(U)))
              + Ctau3_/fEta_*(T(fvc::grad(U)) & fvc::grad(U))
            )
        )
    )
{
    bound(k_, kMin_);
    // already bounded: bound(epsilon_, epsilonMin_);

    #include "wallNonlinearViscosityI.H"

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> NonlinearKEShih::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)) + nonlinearStress_,
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> NonlinearKEShih::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_))) + nonlinearStress_
        )
    );
}


tmp<fvVectorMatrix> NonlinearKEShih::divDevReff(volVectorField& U) const
{
    return
    (
        fvc::div(nonlinearStress_)
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> NonlinearKEShih::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
        fvc::div(rho*nonlinearStress_)
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool NonlinearKEShih::read()
{
    if (RASModel::read())
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

        kappa_.readIfPresent(coeffDict());
        E_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void NonlinearKEShih::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    // generation term
    tmp<volScalarField> S2 = symm(gradU) && gradU;

    volScalarField G
    (
        GName(),
        Cmu_*sqr(k_)/epsilon_*S2 - (nonlinearStress_ && gradU)
    );

    #include "nonLinearWallFunctionsI.H"

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1_*G*epsilon_/k_
      - fvm::Sp(C2_*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

    #include "wallDissipationI.H"

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

    eta_ = k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU + T(gradU))));
    ksi_ = k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU - T(gradU))));
    Cmu_ = 2.0/(3.0*(A1_ + eta_ + alphaKsi_*ksi_));
    fEta_ = A2_ + pow(eta_, 3.0);

    nut_ = Cmu_*sqr(k_)/epsilon_;

    #include "wallNonlinearViscosityI.H"

    nonlinearStress_ = symm
    (
        pow(k_, 3.0)/sqr(epsilon_)
       *(
            Ctau1_/fEta_
           *(
                (gradU & gradU)
              + (gradU & gradU)().T()
            )
          + Ctau2_/fEta_*(gradU & T(gradU))
          + Ctau3_/fEta_*(T(gradU) & gradU)
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
