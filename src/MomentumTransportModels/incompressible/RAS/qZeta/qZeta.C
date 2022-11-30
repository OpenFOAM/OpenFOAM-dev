/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

defineTypeNameAndDebug(qZeta, 0);
addToRunTimeSelectionTable(RASModel, qZeta, dictionary);

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> qZeta::fMu() const
{
    const volScalarField Rt(q_*k_/(2.0*nu()*zeta_));

    if (anisotropic_)
    {
        return exp((-scalar(2.5) + Rt/20.0)/pow3(scalar(1) + Rt/130.0));
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


void qZeta::correctNut()
{
    nut_ = Cmu_*fMu()*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

qZeta::qZeta
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
    eddyViscosity<incompressible::RASModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
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

    qMin_("qMin", sqrt(kMin_)),
    zetaMin_("zetaMin", epsilonMin_/(2*qMin_)),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
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
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            runTime_.name(),
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
            IOobject::groupName("q", alphaRhoPhi.group()),
            runTime_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt(bound(k_, kMin_)),
        k_.boundaryField().types()
    ),

    zeta_
    (
        IOobject
        (
            IOobject::groupName("zeta", alphaRhoPhi.group()),
            runTime_.name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        bound(epsilon_, epsilonMin_)/(2.0*q_),
        epsilon_.boundaryField().types()
    )
{
    bound(zeta_, zetaMin_);

    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool qZeta::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
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
    if (!turbulence_)
    {
        return;
    }

    eddyViscosity<incompressible::RASModel>::correct();

    volScalarField G(GName(), nut_/(2.0*q_)*2*magSqr(symm(fvc::grad(U_))));
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

    zetaEqn.ref().relax();
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

    qEqn.ref().relax();
    solve(qEqn);
    bound(q_, qMin_);


    // Re-calculate k and epsilon
    k_ = sqr(q_);
    k_.correctBoundaryConditions();

    epsilon_ = 2*q_*zeta_;
    epsilon_.correctBoundaryConditions();

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
