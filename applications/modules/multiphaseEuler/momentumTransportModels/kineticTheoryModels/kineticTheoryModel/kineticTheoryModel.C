/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
#include "phaseSystem.H"
#include "fvcDdt.H"
#include "fvcSup.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::phaseModel&
Foam::RASModels::kineticTheoryModel::continuousPhase() const
{
    const phaseSystem& fluid = phase_.fluid();

    if (continuousPhaseName_ == word::null)
    {
        if (fluid.movingPhases().size() != 2)
        {
            FatalIOErrorInFunction(coeffDict())
                << "Continuous phase name must be specified "
                << "when there are more than two moving phases."
                << exit(FatalIOError);
        }

        forAll(fluid.movingPhases(), movingPhasei)
        {
            const phaseModel& otherPhase = fluid.movingPhases()[movingPhasei];

            if (&otherPhase != &phase_)
            {
                return otherPhase;
            }
        }
    }

    return fluid.phases()[continuousPhaseName_];
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::kineticTheoryModel::kineticTheoryModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    eddyViscosity<RASModel<phaseCompressible::momentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    phase_(refCast<const phaseModel>(viscosity)),

    continuousPhaseName_
    (
        coeffDict().lookupOrDefault("continuousPhase", word::null)
    ),

    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            coeffDict()
        )
    ),
    conductivityModel_
    (
        kineticTheoryModels::conductivityModel::New
        (
            coeffDict()
        )
    ),
    radialModel_
    (
        kineticTheoryModels::radialModel::New
        (
            coeffDict()
        )
    ),
    granularPressureModel_
    (
        kineticTheoryModels::granularPressureModel::New
        (
            coeffDict()
        )
    ),
    frictionalStressModel_
    (
        kineticTheoryModels::frictionalStressModel::New
        (
            coeffDict()
        )
    ),

    equilibrium_(coeffDict().lookup("equilibrium")),
    e_("e", dimless, coeffDict()),
    alphaMinFriction_
    (
        "alphaMinFriction",
        dimless,
        coeffDict()
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        coeffDict()
    ),

    maxNut_
    (
        "maxNut",
        dimensionSet(0, 2, -1, 0, 0),
        coeffDict().lookupOrDefault<scalar>("maxNut", 1000)
    ),

    Theta_
    (
        IOobject
        (
            IOobject::groupName("Theta", phase_.name()),
            U.time().name(),
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
            IOobject::groupName(typedName("lambda"), phase_.name()),
            U.time().name(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
    ),

    gs0_
    (
        IOobject
        (
            IOobject::groupName(typedName("gs0"), phase_.name()),
            U.time().name(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimensionSet(0, 0, 0, 0, 0), 0)
    ),

    kappa_
    (
        IOobject
        (
            IOobject::groupName(typedName("kappa"), phase_.name()),
            U.time().name(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), 0)
    ),

    nuFric_
    (
        IOobject
        (
            IOobject::groupName(typedName("nuFric"), phase_.name()),
            U.time().name(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar(dimensionSet(0, 2, -1, 0, 0), 0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::RASModels::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::kineticTheoryModel::read()
{
    if
    (
        eddyViscosity<RASModel<phaseCompressible::momentumTransportModel>>::
        read()
    )
    {
        const dictionary& coeffDict = this->coeffDict();
        coeffDict.lookup("equilibrium") >> equilibrium_;
        e_.readIfPresent(coeffDict);
        alphaMinFriction_.readIfPresent(coeffDict);

        viscosityModel_->read(coeffDict);
        conductivityModel_->read(coeffDict);
        radialModel_->read(coeffDict);
        granularPressureModel_->read(coeffDict);
        frictionalStressModel_->read(coeffDict);

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
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::kineticTheoryModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::kineticTheoryModel::omega() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::kineticTheoryModel::R() const
{
    return tmp<volSymmTensorField>
    (
        volSymmTensorField::New
        (
            IOobject::groupName("R", U_.group()),
          - (nut_ + nuFric_)*dev(twoSymm(fvc::grad(U_)))
          - (lambda_*fvc::div(phi_))*symmTensor::I
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::kineticTheoryModel::pPrime() const
{
    const volScalarField& rho = phase_.rho();

    tmp<volScalarField> tpPrime
    (
        volScalarField::New
        (
            IOobject::groupName("pPrime", U_.group()),
            Theta_
           *granularPressureModel_->granularPressureCoeffPrime
            (
                alpha_,
                radialModel_->g0
                (
                    alpha_,
                    alphaMinFriction_,
                    phase_.alphaMax()
                ),
                radialModel_->g0prime
                (
                    alpha_,
                    alphaMinFriction_,
                    phase_.alphaMax()
                ),
                rho,
                e_
            )
         +  frictionalStressModel_->frictionalPressurePrime
            (
                phase_,
                alphaMinFriction_,
                phase_.alphaMax()
            )
        )
    );

    volScalarField::Boundary& bpPrime = tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::kineticTheoryModel::pPrimef() const
{
    return surfaceScalarField::New
    (
        IOobject::groupName("pPrimef", U_.group()),
        fvc::interpolate(pPrime())
    );
}


Foam::tmp<Foam::surfaceVectorField>
Foam::RASModels::kineticTheoryModel::devTau() const
{
    const surfaceScalarField rhoNuEff(fvc::interpolate(rho_*(nut_ + nuFric_)));

    return tmp<surfaceVectorField>
    (
        surfaceVectorField::New
        (
            IOobject::groupName("devTau", U_.group()),
          - rhoNuEff
           *(
               fvc::dotInterpolate(mesh().nf(), dev2(T(fvc::grad(U_))))
             + fvc::snGrad(U_)
            )
          - fvc::interpolate((rho_*lambda_)*fvc::div(phi_))*mesh().Sf()
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::kineticTheoryModel::divDevTau
(
    volVectorField& U
) const
{
    const surfaceScalarField rhoNuEff(fvc::interpolate(rho_*(nut_ + nuFric_)));

    return
    (
      - fvc::div
        (
            rhoNuEff
           *fvc::dotInterpolate(mesh().Sf(), dev2(T(fvc::grad(U))))
          + fvc::interpolate((rho_*lambda_)*fvc::div(phi_))*mesh().Sf()
        )
      - fvm::laplacian(rhoNuEff, U)
    );
}


void Foam::RASModels::kineticTheoryModel::correct()
{
    // Local references
    const volScalarField alpha(max(alpha_, scalar(0)));
    const phaseSystem& fluid = phase_.fluid();
    const phaseModel& continuousPhase = this->continuousPhase();
    const volScalarField& rho = phase_.rho();
    const surfaceScalarField& alphaRhoPhi = alphaRhoPhi_;
    const volVectorField& U = U_;

    tmp<volVectorField> tUc(continuousPhase.U());
    const volVectorField& Uc = tUc();

    const scalar sqrtPi = sqrt(constant::mathematical::pi);
    const dimensionedScalar ThetaSmall("ThetaSmall", Theta_.dimensions(), 1e-6);
    const dimensionedScalar ThetaSmallSqrt(sqrt(ThetaSmall));

    tmp<volScalarField> tda(phase_.d());
    const volScalarField& da = tda();

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU(tgradU());
    const volSymmTensorField D(symm(gradU));

    // Calculating the radial distribution function
    gs0_ = radialModel_->g0(alpha, alphaMinFriction_, phase_.alphaMax());

    if (!equilibrium_)
    {
        // Particle viscosity (Table 3.2, p.47)
        nut_ = viscosityModel_->nu(alpha, Theta_, gs0_, rho, da, e_);

        const volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

        // Bulk viscosity  p. 45 (Lun et al. 1984).
        lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1 + e_)*ThetaSqrt/sqrtPi;

        // Stress tensor, Definitions, Table 3.1, p. 43
        const volSymmTensorField tau
        (
            rho*(2*nut_*D + (lambda_ - (2.0/3.0)*nut_)*tr(D)*I)
        );

        // Dissipation (Eq. 3.24, p.50)
        const volScalarField gammaCoeff
        (
            "gammaCoeff",
            12*(1 - sqr(e_))
           *max(sqr(alpha), residualAlpha_)
           *rho*gs0_*(1.0/da)*ThetaSqrt/sqrtPi
        );

        // Drag
        const dispersedPhaseInterface interface(phase_, continuousPhase);
        const volScalarField beta
        (
            fluid.foundInterfacialModel<dragModel>(interface)
          ? fluid.lookupInterfacialModel<dragModel>(interface).K()
          : volScalarField::New
            (
                "beta",
                phase_.mesh(),
                dimensionedScalar(dragModel::dimK, 0)
            )
        );

        // Eq. 3.25, p. 50 Js = J1 - J2
        const volScalarField J1("J1", 3*beta);
        const volScalarField J2
        (
            "J2",
            0.25*sqr(beta)*da*magSqr(U - Uc)
           /(
               max(alpha, residualAlpha_)*rho
              *sqrtPi*(ThetaSqrt + ThetaSmallSqrt)
            )
        );

        // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
        const volScalarField PsCoeff
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

        const Foam::fvModels& fvModels(Foam::fvModels::New(mesh_));
        const Foam::fvConstraints& fvConstraints
        (
            Foam::fvConstraints::New(mesh_)
        );

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
          - fvm::SuSp((PsCoeff*I) && gradU, Theta_)
          + (tau && gradU)
          + fvm::Sp(-gammaCoeff, Theta_)
          + fvm::Sp(-J1, Theta_)
          + fvm::Sp(J2/(Theta_ + ThetaSmall), Theta_)
          + fvModels.source(alpha, rho, Theta_)
        );

        ThetaEqn.relax();
        fvConstraints.constrain(ThetaEqn);
        ThetaEqn.solve();
        fvConstraints.constrain(Theta_);
    }
    else
    {
        // Equilibrium => dissipation == production
        // Eq. 4.14, p.82
        const volScalarField K1("K1", 2*(1 + e_)*rho*gs0_);
        const volScalarField K3
        (
            "K3",
            0.5*da*rho*
            (
                (sqrtPi/(3*(3.0 - e_)))
               *(1 + 0.4*(1 + e_)*(3*e_ - 1)*alpha*gs0_)
               +1.6*alpha*gs0_*(1 + e_)/sqrtPi
            )
        );

        const volScalarField K2
        (
            "K2",
            4*da*rho*(1 + e_)*alpha*gs0_/(3*sqrtPi) - 2*K3/3.0
        );

        const volScalarField K4("K4", 12*(1 - sqr(e_))*rho*gs0_/(da*sqrtPi));

        const volScalarField trD
        (
            "trD",
            alpha/(alpha + residualAlpha_)
           *fvc::div(phi_)
        );
        const volScalarField tr2D("tr2D", sqr(trD));
        const volScalarField trD2("trD2", tr(D & D));

        const volScalarField t1("t1", K1*alpha + rho);
        const volScalarField l1("l1", -t1*trD);
        const volScalarField l2("l2", sqr(t1)*tr2D);
        const volScalarField l3
        (
            "l3",
            4.0
           *K4
           *alpha
           *(2*K3*trD2 + K2*tr2D)
        );

        Theta_ = sqr
        (
            (l1 + sqrt(l2 + l3))
           /(2*max(alpha, residualAlpha_)*K4)
        );

        kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rho, da, e_);
    }

    Theta_.max(0);
    Theta_.min(100);

    {
        // particle viscosity (Table 3.2, p.47)
        nut_ = viscosityModel_->nu(alpha, Theta_, gs0_, rho, da, e_);

        const volScalarField ThetaSqrt("sqrtTheta", sqrt(Theta_));

        // Bulk viscosity  p. 45 (Lun et al. 1984).
        lambda_ = (4.0/3.0)*sqr(alpha)*da*gs0_*(1 + e_)*ThetaSqrt/sqrtPi;

        // Frictional pressure
        const volScalarField pf
        (
            frictionalStressModel_->frictionalPressure
            (
                phase_,
                alphaMinFriction_,
                phase_.alphaMax()
            )
        );

        // Limit viscosity
        nut_.min(maxNut_);

        nuFric_ = min
        (
            frictionalStressModel_->nu
            (
                phase_,
                alphaMinFriction_,
                phase_.alphaMax(),
                pf/rho,
                D
            ),
            maxNut_ - nut_
        );
    }

    if (debug)
    {
        Info<< typeName << ':' << nl
            << "    max(Theta) = " << max(Theta_).value() << nl
            << "    max(nut) = " << max(nut_).value() << endl;
    }
}


// ************************************************************************* //
