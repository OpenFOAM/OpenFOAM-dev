/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "v2f.H"
#include "fixedValueFvPatchField.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(v2f, 0);
addToRunTimeSelectionTable(RASModel, v2f, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

wordList v2f::RBoundaryTypes() const
{
    const volScalarField::GeometricBoundaryField& bf(k_.boundaryField());

    wordList bTypes
    (
        bf.size(),
        zeroGradientFvPatchField<symmTensor>::typeName
    );

    forAll(bf, patchI)
    {
        if (bf[patchI].fixesValue())
        {
            bTypes[patchI] = fixedValueFvPatchField<symmTensor>::typeName;
        }
    }

    return bTypes;
}


tmp<volScalarField> v2f::davidsonCorrectNut
(
    const tmp<volScalarField>& value
) const
{
    return min(CmuKEps_*sqr(k_)/epsilon_, value);
}


tmp<volScalarField> v2f::Ts() const
{
    return max(k_/epsilon_, 6.0*sqrt(nu()/epsilon_));
}


tmp<volScalarField> v2f::Ls() const
{
    return CL_*max(pow(k_, 1.5)/epsilon_, Ceta_*pow025(pow3(nu())/epsilon_));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

v2f::v2f
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.22
        )
    ),
    CmuKEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CmuKEps",
            coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.4
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            0.3
        )
    ),
    CL_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CL",
            coeffDict_,
            0.23
        )
    ),
    Ceta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceta",
            coeffDict_,
            70.0
        )
    ),
    Ceps2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps2",
            coeffDict_,
            1.9
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
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
    v2_
    (
        IOobject
        (
            "v2",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    f_
    (
        IOobject
        (
            "f",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    v2Min_(dimensionedScalar("v2Min", v2_.dimensions(), SMALL)),
    fMin_(dimensionedScalar("fMin", f_.dimensions(), 0.0))
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);
    bound(v2_, v2Min_);
    bound(f_, fMin_);

    nut_ = davidsonCorrectNut(Cmu_*v2_*Ts());
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> v2f::R() const
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
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            RBoundaryTypes()
        )
    );
}


tmp<volSymmTensorField> v2f::devReff() const
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
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> v2f::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> v2f::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}


bool v2f::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        CmuKEps_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        CL_.readIfPresent(coeffDict());
        Ceta_.readIfPresent(coeffDict());
        Ceps2_.readIfPresent(coeffDict());
        sigmaK_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void v2f::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    // use N=6 so that f=0 at walls
    const dimensionedScalar N("N", dimless, 6.0);

    const volTensorField gradU(fvc::grad(U_));
    const volScalarField S2(2*magSqr(dev(symm(gradU))));

    const volScalarField G(GName(), nut_*S2);
    const volScalarField T(Ts());
    const volScalarField L2(type() + ":L2", sqr(Ls()));
    const volScalarField alpha
    (
        "v2f:alpha",
        1.0/T*((C1_ - N)*v2_ - 2.0/3.0*k_*(C1_ - 1.0))
    );


    tmp<volScalarField> Ceps1 =
        1.4*(1.0 + 0.05*min(sqrt(k_/v2_), scalar(100.0)));

    // Update epsilon (and possibly G) at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        Ceps1*G/T
      - fvm::Sp(Ceps2_/T, epsilon_)
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


    // Relaxation function equation
    tmp<fvScalarMatrix> fEqn
    (
      - fvm::laplacian(f_)
     ==
      - fvm::Sp(1.0/L2, f_)
      - 1.0/L2/k_*(alpha - C2_*G)
    );

    fEqn().relax();
    solve(fEqn);
    bound(f_, fMin_);


    // Turbulence stress normal to streamlines equation
    tmp<fvScalarMatrix> v2Eqn
    (
        fvm::ddt(v2_)
      + fvm::div(phi_, v2_)
      - fvm::laplacian(DkEff(), v2_)
      ==
        min(k_*f_, -alpha + C2_*G)
      - fvm::Sp(N*epsilon_/k_, v2_)
    );

    v2Eqn().relax();
    solve(v2Eqn);
    bound(v2_, v2Min_);


    // Re-calculate viscosity
    nut_ = davidsonCorrectNut(Cmu_*v2_*T);
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
