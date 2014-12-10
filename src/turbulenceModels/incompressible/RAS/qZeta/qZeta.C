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

#include "qZeta.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(qZeta, 0);
addToRunTimeSelectionTable(RASModel, qZeta, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> qZeta::fMu() const
{
    const volScalarField Rt(q_*k_/(2.0*nu()*zeta_));

    if (anisotropic_)
    {
        return exp((-scalar(2.5) + Rt/20.0)/pow(scalar(1) + Rt/130.0, 3.0));
    }
    else
    {
        return
            exp(-6.0/sqr(scalar(1) + Rt/50.0))
           *(scalar(1) + 3.0*exp(-Rt/10.0));
    }
}


tmp<volScalarField> qZeta::f2() const
{
    tmp<volScalarField> Rt = q_*k_/(2.0*nu()*zeta_);
    return scalar(1) - 0.3*exp(-sqr(Rt));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

qZeta::qZeta
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
            0.09
        )
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
    sigmaZeta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaZeta",
            coeffDict_,
            1.3
        )
    ),
    anisotropic_
    (
        Switch::lookupOrAddToDict
        (
            "anisotropic",
            coeffDict_,
            false
        )
    ),

    qMin_("qMin", dimVelocity, SMALL),
    zetaMin_("zetaMin", dimVelocity/dimTime, SMALL),

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

    q_
    (
        IOobject
        (
            "q",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt(k_),
        k_.boundaryField().types()
    ),

    zeta_
    (
        IOobject
        (
            "zeta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        epsilon_/(2.0*bound(q_, qMin_)),
        epsilon_.boundaryField().types()
    ),

    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    )
{
    bound(k_, kMin_);
    bound(epsilon_, epsilonMin_);
    // already bounded: bound(q_, qMin_);
    bound(zeta_, zetaMin_);

    nut_ = Cmu_*fMu()*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> qZeta::R() const
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
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> qZeta::devReff() const
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


tmp<fvVectorMatrix> qZeta::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> qZeta::divDevRhoReff
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


bool qZeta::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmaZeta_.readIfPresent(coeffDict());
        anisotropic_.readIfPresent("anisotropic", coeffDict());

        qMin_.readIfPresent(*this);
        zetaMin_.readIfPresent(*this);

        return true;
    }
    else
    {
        return false;
    }
}


void qZeta::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    tmp<volScalarField> S2 = 2*magSqr(symm(fvc::grad(U_)));

    volScalarField G(GName(), nut_/(2.0*q_)*S2);
    const volScalarField E(nu()*nut_/q_*fvc::magSqrGradGrad(U_));


    // Zeta equation

    tmp<fvScalarMatrix> zetaEqn
    (
        fvm::ddt(zeta_)
      + fvm::div(phi_, zeta_)
      - fvm::laplacian(DzetaEff(), zeta_)
     ==
        (2.0*C1_ - 1)*G*zeta_/q_
      - fvm::SuSp((2.0*C2_*f2() - dimensionedScalar(1.0))*zeta_/q_, zeta_)
      + E
    );

    zetaEqn().relax();
    solve(zetaEqn);
    bound(zeta_, zetaMin_);


    // q equation

    tmp<fvScalarMatrix> qEqn
    (
        fvm::ddt(q_)
      + fvm::div(phi_, q_)
      - fvm::laplacian(DqEff(), q_)
     ==
        G - fvm::Sp(zeta_/q_, q_)
    );

    qEqn().relax();
    solve(qEqn);
    bound(q_, qMin_);


    // Re-calculate k and epsilon
    k_ = sqr(q_);
    k_.correctBoundaryConditions();

    epsilon_ = 2*q_*zeta_;
    epsilon_.correctBoundaryConditions();


    // Re-calculate viscosity
    nut_ = Cmu_*fMu()*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
