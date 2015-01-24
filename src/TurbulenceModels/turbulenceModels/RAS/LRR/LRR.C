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

#include "LRR.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void LRR<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = this->Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
LRR<BasicTurbulenceModel>::LRR
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    ReynoldsStress<RASModel<BasicTurbulenceModel> >
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

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    Clrr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clrr1",
            this->coeffDict_,
            1.8
        )
    ),
    Clrr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clrr2",
            this->coeffDict_,
            0.6
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.25
        )
    ),
    Ceps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ceps",
            this->coeffDict_,
            0.15
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*tr(this->R_)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    if (type == typeName)
    {
        this->boundNormalStress(this->R_);
        bound(epsilon_, this->epsilonMin_);
        k_ = 0.5*tr(this->R_);
        correctNut();
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool LRR<BasicTurbulenceModel>::read()
{
    if (ReynoldsStress<RASModel<BasicTurbulenceModel> >::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        Clrr1_.readIfPresent(this->coeffDict());
        Clrr2_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());
        Ceps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void LRR<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volSymmTensorField& R = this->R_;

    ReynoldsStress<RASModel<BasicTurbulenceModel> >::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    volSymmTensorField P(-twoSymm(R & gradU));
    volScalarField G(this->GName(), 0.5*mag(tr(P)));

    // Update epsilon and G at the wall
    epsilon_.boundaryField().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(Ceps_*alpha*rho*(k_/epsilon_)*R, epsilon_)
     ==
        C1_*alpha*rho*G*epsilon_/k_
      - fvm::Sp(C2_*alpha*rho*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

    epsEqn().boundaryManipulate(epsilon_.boundaryField());

    solve(epsEqn);
    bound(epsilon_, this->epsilonMin_);


    // Reynolds stress equation

    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];
                P[faceCelli] *= min
                (
                    G[faceCelli]/(0.5*mag(tr(P[faceCelli])) + SMALL),
                    1.0
                );
            }
        }
    }


    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(alpha, rho, R)
      + fvm::div(alphaRhoPhi, R)
      - fvm::laplacian(Cs_*alpha*rho*(k_/epsilon_)*R, R)
      + fvm::Sp(Clrr1_*alpha*rho*epsilon_/k_, R)
      ==
        alpha*rho*P
      - (2.0/3.0*(1 - Clrr1_)*I)*alpha*rho*epsilon_
      - Clrr2_*alpha*rho*dev(P)
    );

    REqn().relax();
    solve(REqn);

    this->boundNormalStress(R);

    k_ = 0.5*tr(R);

    correctNut();

    // Correct wall shear-stresses when applying wall-functions
    this->correctWallShearStress(R);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
