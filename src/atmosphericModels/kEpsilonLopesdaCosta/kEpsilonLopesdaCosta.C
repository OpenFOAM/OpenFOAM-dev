/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2023 OpenFOAM Foundation
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

#include "kEpsilonLopesdaCosta.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "porosityForce.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void kEpsilonLopesdaCosta<BasicMomentumTransportModel>::setPorosityCoefficient
(
    volScalarField::Internal& C,
    const porosityModels::powerLawLopesdaCosta& pm
)
{
    if (pm.dict().found(C.name()))
    {
        const scalar Cpm = pm.dict().lookup<scalar>(C.name());

        const labelList& cells = this->mesh_.cellZones()[pm.zoneName()];

        forAll(cells, i)
        {
            const label celli = cells[i];
            C[celli] = Cpm;
        }
    }
}


template<class BasicMomentumTransportModel>
void kEpsilonLopesdaCosta<BasicMomentumTransportModel>::setCdAv
(
    volScalarField::Internal& C,
    const porosityModels::powerLawLopesdaCosta& pm
)
{
    if (pm.dict().found(C.name()))
    {
        const scalarField& Av = pm.Av();

        const scalar Cpm = pm.dict().lookup<scalar>(C.name());

        const labelList& cells = this->mesh_.cellZones()[pm.zoneName()];

        forAll(cells, i)
        {
            const label celli = cells[i];
            C[celli] = Cpm*Av[celli];
        }
    }
}


template<class BasicMomentumTransportModel>
void kEpsilonLopesdaCosta<BasicMomentumTransportModel>::
setPorosityCoefficients()
{
    const PtrListDictionary<fvModel>& fvModels
    (
        Foam::fvModels::New(this->mesh_)
    );

    forAll(fvModels, i)
    {
        if (isA<fv::porosityForce>(fvModels[i]))
        {
            const fv::porosityForce& eps =
                refCast<const fv::porosityForce>(fvModels[i]);

            if (isA<porosityModels::powerLawLopesdaCosta>(eps.model()))
            {
                const porosityModels::powerLawLopesdaCosta& pm =
                    refCast<const porosityModels::powerLawLopesdaCosta>
                    (
                        eps.model()
                    );

                setPorosityCoefficient(Cmu_, pm);
                setPorosityCoefficient(C1_, pm);
                setPorosityCoefficient(C2_, pm);
                setPorosityCoefficient(sigmak_, pm);
                setPorosityCoefficient(sigmaEps_, pm);

                setCdAv(CdAv_, pm);
                setPorosityCoefficient(betap_, pm);
                setPorosityCoefficient(betad_, pm);
                setPorosityCoefficient(C4_, pm);
                setPorosityCoefficient(C5_, pm);
            }
        }
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField>
kEpsilonLopesdaCosta<BasicMomentumTransportModel>::boundEpsilon()
{
    tmp<volScalarField> tCmuk2(Cmu_*sqr(k_));
    epsilon_ = max(epsilon_, tCmuk2()/(this->nutMaxCoeff_*this->nu()));
    return tCmuk2;
}


template<class BasicMomentumTransportModel>
void kEpsilonLopesdaCosta<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = boundEpsilon()/epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kEpsilonLopesdaCosta<BasicMomentumTransportModel>::kSource
(
    const volScalarField::Internal& magU,
    const volScalarField::Internal& magU3
) const
{
    return fvm::Su(CdAv_*(betap_*magU3 - betad_*magU*k_()), k_);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
kEpsilonLopesdaCosta<BasicMomentumTransportModel>::epsilonSource
(
    const volScalarField::Internal& magU,
    const volScalarField::Internal& magU3
) const
{
    return fvm::Su
    (
        CdAv_
       *(C4_*betap_*epsilon_()/k_()*magU3 - C5_*betad_*magU*epsilon_()),
        epsilon_
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kEpsilonLopesdaCosta<BasicMomentumTransportModel>::kEpsilonLopesdaCosta
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
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
        IOobject
        (
            "Cmu",
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        IOobject
        (
            "C1",
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        IOobject
        (
            "C2",
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    sigmak_
    (
        IOobject
        (
            "sigmak",
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        IOobject
        (
            "sigmaEps",
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),

    CdAv_
    (
        IOobject
        (
            "CdAv",
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(dimless/dimLength, 0)
    ),
    betap_
    (
        IOobject
        (
            "betap",
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(dimless, 0)
    ),
    betad_
    (
        IOobject
        (
            "betad",
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(dimless, 0)
    ),
    C4_
    (
        IOobject
        (
            "C4",
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(dimless, 0)
    ),
    C5_
    (
        IOobject
        (
            "C5",
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(dimless, 0)
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    boundEpsilon();

    if (type == typeName)
    {
        this->printCoeffs(type);
    }

    setPorosityCoefficients();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kEpsilonLopesdaCosta<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void kEpsilonLopesdaCosta<BasicMomentumTransportModel>::correct()
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
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    volScalarField::Internal magU(mag(U));
    volScalarField::Internal magU3(pow3(magU));

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*G*epsilon_()/k_()
      - fvm::SuSp(((2.0/3.0)*C1_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + epsilonSource(magU, magU3)
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    boundEpsilon();

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      + kSource(magU, magU3)
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
