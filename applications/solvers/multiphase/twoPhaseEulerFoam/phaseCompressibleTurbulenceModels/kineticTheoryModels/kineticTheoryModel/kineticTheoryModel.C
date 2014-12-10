/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "kineticTheoryModel.H"
#include "mathematicalConstants.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::kineticTheoryModel::kineticTheoryModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<PhaseCompressibleTurbulenceModel<phaseModel> > >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        phase,
        propertiesName
    ),

    phase_(phase),

    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            this->coeffDict_
        )
    ),
    conductivityModel_
    (
        kineticTheoryModels::conductivityModel::New
        (
            this->coeffDict_
        )
    ),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            this->coeffDict_
        )
    ),
    granularPressureModel_
    (
        kineticTheoryModels::granularPressureModel::New
        (
            this->coeffDict_
        )
    ),
    frictionalStressModel_
    (
        kineticTheoryModels::frictionalStressModel::New
        (
            this->coeffDict_
        )
    ),

    equilibrium_(this->coeffDict_.lookup("equilibrium")),
    e_("e", dimless, this->coeffDict_.lookup("e")),
    alphaMax_("alphaMax", dimless, this->coeffDict_.lookup("alphaMax")),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        this->coeffDict_.lookup("alphaMinFriction")
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        this->coeffDict_.lookup("residualAlpha")
    ),

    Theta_
    (
        IOobject
        (
            IOobject::groupName("Theta", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    lambda_
    (
        IOobject
        (
            IOobject::groupName("lambda", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0), 0.0)
    ),

    gs0_
    (
        IOobject
        (
            IOobject::groupName("gs0", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 0.0)
    ),

    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", phase.name()),
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::kineticTheoryModel::read()
{
    if
    (
        eddyViscosity
        <
            RASModel<PhaseCompressibleTurbulenceModel<phaseModel> >
        >::read()
    )
    {
        this->coeffDict().lookup("equilibrium") >> equilibrium_;
        e_.readIfPresent(this->coeffDict());
        alphaMax_.readIfPresent(this->coeffDict());
        alphaMinFriction_.readIfPresent(this->coeffDict());

        viscosityModel_->read();
        conductivityModel_->read();
        radialModel_->read();
        granularPressureModel_->read();
        frictionalStressModel_->read();

        return true;
    }
    else
    {
        return false;
    }
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::kineticTheoryModel::k() const
{
    notImplemented("kineticTheoryModel::k()");
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::kineticTheoryModel::epsilon() const
{
    notImplemented("kineticTheoryModel::epsilon()");
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::kineticTheoryModel::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (this->nut_)*dev(twoSymm(fvc::grad(this->U_)))
          - (lambda_*fvc::div(this->phi_))*symmTensor::I
        )
    );
}


/*
Foam::tmp<Foam::volScalarField>
Foam::RASModels::kineticTheoryModel::pp() const
{

    // Particle pressure coefficient
    // Coefficient in front of Theta (Eq. 3.22, p. 45)
    volScalarField PsCoeff
    (
        granularPressureModel_->granularPressureCoeff
        (
            alpha,
            gs0,
            rho,
            e_
        )
    );

    // Frictional pressure
    volScalarField pf
    (
        frictionalStressModel_->frictionalPressure
        (
            alpha,
            alphaMinFriction_,
            alphaMax_
        )
    );

    // Return total particle pressure
    return PsCoeff*Theta_ + pf;
}
*/


Foam::tmp<Foam::volScalarField>
Foam::RASModels::kineticTheoryModel::pPrime() const
{
    // Local references
    const volScalarField& alpha = this->alpha_;
    const volScalarField& rho = phase_.rho();

    return
    (
        Theta_
       *granularPressureModel_->granularPressureCoeffPrime
        (
            alpha,
            radialModel_->g0(alpha, alphaMinFriction_, alphaMax_),
            radialModel_->g0prime(alpha, alphaMinFriction_, alphaMax_),
            rho,
            e_
        )
     +  frictionalStressModel_->frictionalPressurePrime
        (
            alpha,
            alphaMinFriction_,
            alphaMax_
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::kineticTheoryModel::pPrimef() const
{
    // Local references
    const volScalarField& alpha = this->alpha_;
    const volScalarField& rho = phase_.rho();

    return fvc::interpolate
    (
        Theta_
       *granularPressureModel_->granularPressureCoeffPrime
        (
            alpha,
            radialModel_->g0(alpha, alphaMinFriction_, alphaMax_),
            radialModel_->g0prime(alpha, alphaMinFriction_, alphaMax_),
            rho,
            e_
        )
     +  frictionalStressModel_->frictionalPressurePrime
        (
            alpha,
            alphaMinFriction_,
            alphaMax_
        )
    );
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::kineticTheoryModel::devRhoReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (this->rho_*this->nut_)
           *dev(twoSymm(fvc::grad(this->U_)))
          - ((this->rho_*lambda_)*fvc::div(this->phi_))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::kineticTheoryModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(this->rho_*this->nut_, U)
      - fvc::div
        (
            (this->rho_*this->nut_)*dev2(T(fvc::grad(U)))
          + ((this->rho_*lambda_)*fvc::div(this->phi_))
           *dimensioned<symmTensor>("I", dimless, symmTensor::I)
        )
    );
}


void Foam::RASModels::kineticTheoryModel::correct()
{
    // Local references
    volScalarField alpha(max(this->alpha_, scalar(0)));
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volVectorField& Uc_ = phase_.fluid().otherPhase(phase_).U();

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1.0e-6);
    dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));

    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();

    tmp<volTensorField> tgradU(fvc::grad(this->U_));
    const volTensorField& gradU(tgradU());
    volSymmTensorField D(symm(gradU));

    // Calculating the radial distribution function
    gs0_ = radialModel_->g0(alpha, alphaMinFriction_, alphaMax_);

    if (!equilibrium_)
    {
        // Particle viscosity (Table 3.2, p.47)
        nut_ = viscosityModel_->nu(alpha, Theta_, gs0_, rho, da, e_);

        volScalarField ThetaSqrt(sqrt(Theta_));

        // Bulk viscosity  p. 45 (Lun et al. 1984).
        lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

        // Stress tensor, Definitions, Table 3.1, p. 43
        volSymmTensorField tau
        (
            rho*(2.0*nut_*D + (lambda_ - (2.0/3.0)*nut_)*tr(D)*I)
        );

        // Dissipation (Eq. 3.24, p.50)
        volScalarField gammaCoeff
        (
            12.0*(1.0 - sqr(e_))
           *max(sqr(alpha), residualAlpha_)
           *rho*gs0_*(1.0/da)*ThetaSqrt/sqrtPi
        );

        // Drag
        volScalarField beta(phase_.fluid().drag(phase_).K());

        // Eq. 3.25, p. 50 Js = J1 - J2
        volScalarField J1(3.0*beta);
        volScalarField J2
        (
            0.25*sqr(beta)*da*magSqr(U - Uc_)
           /(
               max(alpha, residualAlpha_)*rho
              *sqrtPi*(ThetaSqrt + ThetaSmallSqrt)
            )
        );

        // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
        volScalarField PsCoeff
        (
            granularPressureModel_->granularPressureCoeff
            (
                alpha,
                gs0_,
                rho,
                e_
            )
        );

        // 'thermal' conductivity (Table 3.3, p. 49)
        kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);

        // Construct the granular temperature equation (Eq. 3.20, p. 44)
        // NB. note that there are two typos in Eq. 3.20:
        //     Ps should be without grad
        //     the laplacian has the wrong sign
        fvScalarMatrix ThetaEqn
        (
            1.5*
            (
                fvm::ddt(alpha, rho, Theta_)
              + fvm::div(alphaRhoPhi, Theta_)
              - fvc::Sp(fvc::ddt(alpha, rho) + fvc::div(alphaRhoPhi), Theta_)
            )
          - fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
         ==
            fvm::SuSp(-((PsCoeff*I) && gradU), Theta_)
          + (tau && gradU)
          + fvm::Sp(-gammaCoeff, Theta_)
          + fvm::Sp(-J1, Theta_)
          + fvm::Sp(J2/(Theta_ + ThetaSmall), Theta_)
        );

        ThetaEqn.relax();
        ThetaEqn.solve();
    }
    else
    {
        // Equilibrium => dissipation == production
        // Eq. 4.14, p.82
        volScalarField K1(2.0*(1.0 + e_)*rho*gs0_);
        volScalarField K3
        (
            0.5*da*rho*
            (
                (sqrtPi/(3.0*(3.0 - e_)))
               *(1.0 + 0.4*(1.0 + e_)*(3.0*e_ - 1.0)*alpha*gs0_)
               +1.6*alpha*gs0_*(1.0 + e_)/sqrtPi
            )
        );

        volScalarField K2
        (
            4.0*da*rho*(1.0 + e_)*alpha*gs0_/(3.0*sqrtPi) - 2.0*K3/3.0
        );

        volScalarField K4(12.0*(1.0 - sqr(e_))*rho*gs0_/(da*sqrtPi));

        volScalarField trD
        (
            alpha/(alpha + residualAlpha_)
           *fvc::div(this->phi_)
        );
        volScalarField tr2D(sqr(trD));
        volScalarField trD2(tr(D & D));

        volScalarField t1(K1*alpha + rho);
        volScalarField l1(-t1*trD);
        volScalarField l2(sqr(t1)*tr2D);
        volScalarField l3
        (
            4.0
           *K4
           *alpha
           *(2.0*K3*trD2 + K2*tr2D)
        );

        Theta_ = sqr
        (
            (l1 + sqrt(l2 + l3))
           /(2.0*max(alpha, residualAlpha_)*K4)
        );

        kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);
    }

    Theta_.max(0);
    Theta_.min(100);

    {
        // particle viscosity (Table 3.2, p.47)
        nut_ = viscosityModel_->nu(alpha, Theta_, gs0_, rho, da, e_);

        volScalarField ThetaSqrt(sqrt(Theta_));

        // Bulk viscosity  p. 45 (Lun et al. 1984).
        lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

        // Frictional pressure
        volScalarField pf
        (
            frictionalStressModel_->frictionalPressure
            (
                alpha,
                alphaMinFriction_,
                alphaMax_
            )
        );

        // Add frictional shear viscosity, Eq. 3.30, p. 52
        nut_ += frictionalStressModel_->nu
        (
            alpha,
            alphaMax_,
            pf/rho,
            D
        );

        // Limit viscosity
        nut_.min(100);
    }

    if (debug)
    {
        Info<< typeName << ':' << nl
            << "    max(Theta) = " << max(Theta_).value() << nl
            << "    max(nut) = " << max(nut_).value() << endl;
    }
}

// ************************************************************************* //
