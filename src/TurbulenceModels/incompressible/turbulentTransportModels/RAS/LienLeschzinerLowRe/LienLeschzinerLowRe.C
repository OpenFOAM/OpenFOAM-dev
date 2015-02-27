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

#include "LienLeschzinerLowRe.H"
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

defineTypeNameAndDebug(LienLeschzinerLowRe, 0);
addToRunTimeSelectionTable(RASModel, LienLeschzinerLowRe, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<volScalarField> LienLeschzinerLowRe::fMu()
{
    return
        (scalar(1) - exp(-Am_*yStar_))
       /(scalar(1) - exp(-Aepsilon_*yStar_) + SMALL);
}


void LienLeschzinerLowRe::correctNut()
{
    nut_ = Cmu_*fMu()*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LienLeschzinerLowRe::LienLeschzinerLowRe
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
    eddyViscosity<incompressible::RASModel>
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
    Cmu_
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

    yStar_(sqrt(k_)*y_/nu() + SMALL)
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

bool LienLeschzinerLowRe::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        Cmu_.readIfPresent(coeffDict());
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


void LienLeschzinerLowRe::correct()
{
    eddyViscosity<incompressible::RASModel>::correct();

    if (!turbulence_)
    {
        return;
    }

    tmp<volTensorField> tgradU = fvc::grad(U_);
    volScalarField G(GName(), nut_*(tgradU() && twoSymm(tgradU())));
    tgradU.clear();

    scalar Cmu75 = pow(Cmu_.value(), 0.75);
    yStar_ = sqrt(k_)*y_/nu() + SMALL;
    tmp<volScalarField> Rt = sqr(k_)/(nu()*epsilon_);
    const volScalarField f2(scalar(1) - 0.3*exp(-sqr(Rt)));

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1_*G*epsilon_/k_

        // E-term
        + C2_*f2*Cmu75*sqrt(k_)
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

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
