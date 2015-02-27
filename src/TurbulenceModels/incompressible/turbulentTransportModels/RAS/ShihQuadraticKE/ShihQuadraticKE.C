/*---------------------------------------------------------------------------*\
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

#include "ShihQuadraticKE.H"
#include "bound.H"
#include "wallFvPatch.H"
#include "nutkWallFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ShihQuadraticKE, 0);
addToRunTimeSelectionTable(RASModel, ShihQuadraticKE, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void ShihQuadraticKE::correctNut()
{
    correctNonlinearStress(fvc::grad(U_));
}


void ShihQuadraticKE::correctNonlinearStress(const volTensorField& gradU)
{
    volSymmTensorField S(symm(gradU));
    volTensorField W(skew(gradU));

    volScalarField sBar((k_/epsilon_)*sqrt(2.0)*mag(symm(gradU)));
    volScalarField wBar((k_/epsilon_)*sqrt(2.0)*mag(skew(gradU)));

    volScalarField Cmu(2.0/(3.0*(Cmu1_ + sBar + Cmu2_*wBar)));

    nut_ = Cmu*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    nonlinearStress_ =
        pow3(k_)/((A1_ + pow3(sBar))*sqr(epsilon_))
       *(
           beta1_*dev(symm(S&S))
         + beta2_*symm((W&S) - (S&W))
         + beta3_*dev(symm(W&W))
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ShihQuadraticKE::ShihQuadraticKE
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
    Cmu1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A1",
            coeffDict_,
            1.25
        )
    ),
    Cmu2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaKsi",
            coeffDict_,
            0.9
        )
    ),
    A1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A2",
            coeffDict_,
            1000.0
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            3.0
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            15.0
        )
    ),
    beta3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta3",
            coeffDict_,
            -19.0
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
    )
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);

    if (type == typeName)
    {
        correctNut();
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool ShihQuadraticKE::read()
{
    if (nonlinearEddyViscosity<incompressible::RASModel>::read())
    {
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        Cmu1_.readIfPresent(coeffDict());
        Cmu2_.readIfPresent(coeffDict());
        A1_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        beta3_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void ShihQuadraticKE::correct()
{
    nonlinearEddyViscosity<incompressible::RASModel>::correct();

    if (!turbulence_)
    {
        return;
    }

    tmp<volTensorField> tgradU = fvc::grad(U_);
    const volTensorField& gradU = tgradU();

    volScalarField G
    (
        GName(),
        (nut_*twoSymm(gradU) - nonlinearStress_) && gradU
    );

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

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

    epsEqn().boundaryManipulate(epsilon_.boundaryField());

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


    // Re-calculate viscosity and non-linear stress
    correctNonlinearStress(gradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
