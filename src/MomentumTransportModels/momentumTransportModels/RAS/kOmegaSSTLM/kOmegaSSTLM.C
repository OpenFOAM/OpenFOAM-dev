/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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

#include "kOmegaSSTLM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> kOmegaSSTLM<BasicMomentumTransportModel>::F1
(
    const volScalarField& CDkOmega
) const
{
    const volScalarField Ry(this->y()*sqrt(this->k_)/this->nu());
    const volScalarField F3(exp(-pow(Ry/120.0, 8)));

    return max(kOmegaSST<BasicMomentumTransportModel>::F1(CDkOmega), F3);
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField> kOmegaSSTLM<BasicMomentumTransportModel>::Pk
(
    const volInternalScalarField& G
) const
{
    return gammaIntEff_*kOmegaSST<BasicMomentumTransportModel>::Pk(G);
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
kOmegaSSTLM<BasicMomentumTransportModel>::epsilonByk
(
    const volInternalScalarField& F1,
    const volInternalScalarField& F2
) const
{
    return
        min(max(gammaIntEff_, scalar(0.1)), scalar(1))
       *kOmegaSST<BasicMomentumTransportModel>::epsilonByk(F1, F2);
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField> kOmegaSSTLM<BasicMomentumTransportModel>::Fthetat
(
    const volInternalScalarField& Us,
    const volInternalScalarField& Omega,
    const volInternalScalarField& nu
) const
{
    const volInternalScalarField& omega = this->omega_();
    const volInternalScalarField& y = this->y()();

    const volInternalScalarField yBydelta
    (
        sqr(Us)
       /max(375*Omega*nu*ReThetat_(), sqr(deltaU_))
    );
    const volInternalScalarField ReOmega(sqr(y)*omega/nu);
    const volInternalScalarField Fwake(exp(-sqr(ReOmega/1e5)));

    return volInternalScalarField::New
    (
        this->groupName("Fthetat"),
        min
        (
            max
            (
                Fwake*exp(-pow4(yBydelta)),
                (1 - sqr((gammaInt_() - 1.0/ce2_)/(1 - 1.0/ce2_)))
            ),
            scalar(1)
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
kOmegaSSTLM<BasicMomentumTransportModel>::ReThetac() const
{
    tmp<volInternalScalarField> tReThetac
    (
        volInternalScalarField::New
        (
            this->groupName("ReThetac"),
            this->mesh_,
            dimless
        )
    );
    volInternalScalarField& ReThetac = tReThetac.ref();

    forAll(ReThetac, celli)
    {
        const scalar ReThetat = ReThetat_[celli];

        ReThetac[celli] =
            ReThetat <= 1870
          ?
            ReThetat
          - 396.035e-2
          + 120.656e-4*ReThetat
          - 868.230e-6*sqr(ReThetat)
          + 696.506e-9*pow3(ReThetat)
          - 174.105e-12*pow4(ReThetat)
          :
            ReThetat - 593.11 - 0.482*(ReThetat - 1870);
    }

    return tReThetac;
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField> kOmegaSSTLM<BasicMomentumTransportModel>::Flength
(
    const volInternalScalarField& nu
) const
{
    tmp<volInternalScalarField> tFlength
    (
        volInternalScalarField::New
        (
            this->groupName("Flength"),
            this->mesh_,
            dimless
        )
    );
    volInternalScalarField& Flength = tFlength.ref();

    const volInternalScalarField& omega = this->omega_();
    const volInternalScalarField& y = this->y()();

    forAll(ReThetat_, celli)
    {
        const scalar ReThetat = ReThetat_[celli];

        if (ReThetat < 400)
        {
            Flength[celli] =
                398.189e-1
              - 119.270e-4*ReThetat
              - 132.567e-6*sqr(ReThetat);
        }
        else if (ReThetat < 596)
        {
            Flength[celli] =
                263.404
              - 123.939e-2*ReThetat
              + 194.548e-5*sqr(ReThetat)
              - 101.695e-8*pow3(ReThetat);
        }
        else if (ReThetat < 1200)
        {
            Flength[celli] = 0.5 - 3e-4*(ReThetat - 596);
        }
        else
        {
            Flength[celli] = 0.3188;
        }

        const scalar Fsublayer =
            exp(-sqr(sqr(y[celli])*omega[celli]/(200*nu[celli])));

        Flength[celli] = Flength[celli]*(1 - Fsublayer) + 40*Fsublayer;
    }

    return tFlength;
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField>
kOmegaSSTLM<BasicMomentumTransportModel>::ReThetat0
(
    const volInternalScalarField& Us,
    const volInternalScalarField& dUsds,
    const volInternalScalarField& nu
) const
{
    tmp<volInternalScalarField> tReThetat0
    (
        volInternalScalarField::New
        (
            this->groupName("ReThetat0"),
            this->mesh_,
            dimless
        )
    );
    volInternalScalarField& ReThetat0 = tReThetat0.ref();

    const volScalarField& k = this->k_;

    label maxIter = 0;

    forAll(ReThetat0, celli)
    {
        const scalar Tu
        (
            max(100*sqrt((2.0/3.0)*k[celli])/Us[celli], scalar(0.027))
        );

        // Initialise lambda to zero.
        // If lambda were cached between time-steps convergence would be faster
        // starting from the previous time-step value.
        scalar lambda = 0;

        scalar lambdaErr;
        scalar thetat;
        label iter = 0;

        do
        {
            // Previous iteration lambda for convergence test
            const scalar lambda0 = lambda;

            if (Tu <= 1.3)
            {
                const scalar Flambda =
                    dUsds[celli] <= 0
                  ?
                    1
                  - (
                     - 12.986*lambda
                     - 123.66*sqr(lambda)
                     - 405.689*pow3(lambda)
                    )*exp(-pow(Tu/1.5, 1.5))
                  :
                    1
                  + 0.275*(1 - exp(-35*lambda))
                   *exp(-Tu/0.5);

                thetat =
                    (1173.51 - 589.428*Tu + 0.2196/sqr(Tu))
                   *Flambda*nu[celli]
                   /Us[celli];
            }
            else
            {
                const scalar Flambda =
                    dUsds[celli] <= 0
                  ?
                    1
                  - (
                      -12.986*lambda
                      -123.66*sqr(lambda)
                      -405.689*pow3(lambda)
                    )*exp(-pow(Tu/1.5, 1.5))
                  :
                    1
                  + 0.275*(1 - exp(-35*lambda))
                   *exp(-2*Tu);

                thetat =
                    331.50*pow((Tu - 0.5658), -0.671)
                   *Flambda*nu[celli]/Us[celli];
            }

            lambda = sqr(thetat)/nu[celli]*dUsds[celli];
            lambda = max(min(lambda, 0.1), -0.1);

            lambdaErr = mag(lambda - lambda0);

            maxIter = max(maxIter, ++iter);

        } while (lambdaErr > lambdaErr_);

        ReThetat0[celli] = max(thetat*Us[celli]/nu[celli], scalar(20));
    }

    if (maxIter > maxLambdaIter_)
    {
        WarningInFunction
            << "Number of lambda iterations exceeds maxLambdaIter("
            << maxLambdaIter_ << ')'<< endl;
    }

    return tReThetat0;
}


template<class BasicMomentumTransportModel>
tmp<volInternalScalarField> kOmegaSSTLM<BasicMomentumTransportModel>::Fonset
(
    const volInternalScalarField& Rev,
    const volInternalScalarField& ReThetac,
    const volInternalScalarField& RT
) const
{
    const volInternalScalarField Fonset1(Rev/(2.193*ReThetac));

    const volInternalScalarField Fonset2
    (
        min(max(Fonset1, pow4(Fonset1)), scalar(2))
    );

    const volInternalScalarField Fonset3(max(1 - pow3(RT/2.5), scalar(0)));

    return volInternalScalarField::New
    (
        this->groupName("Fonset"),
        max(Fonset2 - Fonset3, scalar(0))
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kOmegaSSTLM<BasicMomentumTransportModel>::kOmegaSSTLM
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
    kOmegaSST<BasicMomentumTransportModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    ca1_("ca1", this->typeDict(type), 2),
    ca2_("ca2", this->typeDict(type), 0.06),
    ce1_("ce1", this->typeDict(type), 1),
    ce2_("ce2", this->typeDict(type), 50),
    cThetat_("cThetat", this->typeDict(type), 0.03),
    sigmaThetat_("sigmaThetat", this->typeDict(type), 2),
    lambdaErr_(this->typeDict(type).lookupOrDefault("lambdaErr", 1e-6)),
    maxLambdaIter_(this->typeDict(type).lookupOrDefault("maxLambdaIter", 10)),
    deltaU_("deltaU", dimensions::velocity, small),

    ReThetat_
    (
        IOobject
        (
            this->groupName("ReThetat"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimless
    ),

    gammaInt_
    (
        IOobject
        (
            this->groupName("gammaInt"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimless
    ),

    gammaIntEff_
    (
        IOobject
        (
            this->groupName("gammaIntEff"),
            this->runTime_.name(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(dimless, 0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kOmegaSSTLM<BasicMomentumTransportModel>::read()
{
    if (kOmegaSST<BasicMomentumTransportModel>::read())
    {
        ca1_.readIfPresent(this->typeDict());
        ca2_.readIfPresent(this->typeDict());
        ce1_.readIfPresent(this->typeDict());
        ce2_.readIfPresent(this->typeDict());
        sigmaThetat_.readIfPresent(this->typeDict());
        cThetat_.readIfPresent(this->typeDict());
        this->typeDict().readIfPresent("lambdaErr", lambdaErr_);
        this->typeDict().readIfPresent("maxLambdaIter", maxLambdaIter_);

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void kOmegaSSTLM<BasicMomentumTransportModel>::correctReThetatGammaInt()
{
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const tmp<volScalarField> tnu = this->nu();
    const volInternalScalarField& nu = tnu()();
    const volInternalScalarField& y = this->y()();
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    // Fields derived from the velocity gradient
    tmp<volInternalTensorField> tgradU = fvi::grad(U);
    const volInternalScalarField Omega(sqrt(2*magSqr(skew(tgradU()))));
    const volInternalScalarField S(sqrt(2*magSqr(symm(tgradU()))));
    const volInternalScalarField Us(max(mag(U()), deltaU_));
    const volInternalScalarField dUsds((U() & (U() & tgradU()))/sqr(Us));
    tgradU.clear();

    const volInternalScalarField Fthetat(this->Fthetat(Us, Omega, nu));

    {
        const volInternalScalarField t(500*nu/sqr(Us));
        const volInternalScalarField Pthetat
        (
            alpha()*rho()*(cThetat_/t)*(1 - Fthetat)
        );

        // Transition onset momentum-thickness Reynolds number equation
        tmp<fvScalarMatrix> ReThetatEqn
        (
            fvm::ddt(alpha, rho, ReThetat_)
          + fvm::div(alphaRhoPhi, ReThetat_)
          - fvm::laplacian(alpha*rho*DReThetatEff(), ReThetat_)
         ==
            Pthetat*ReThetat0(Us, dUsds, nu) - fvm::Sp(Pthetat, ReThetat_)
          + fvModels.source(alpha, rho, ReThetat_)
        );

        ReThetatEqn.ref().relax();
        fvConstraints.constrain(ReThetatEqn.ref());
        solve(ReThetatEqn);
        fvConstraints.constrain(ReThetat_);
        bound(ReThetat_, 0);
    }

    const volInternalScalarField ReThetac(this->ReThetac());
    const volInternalScalarField Rev(sqr(y)*S/nu);
    const volInternalScalarField RT(k()/(nu*omega()));

    {
        const volInternalScalarField Pgamma
        (
            alpha()*rho()
           *ca1_*Flength(nu)*S*sqrt(gammaInt_()*Fonset(Rev, ReThetac, RT))
        );

        const volInternalScalarField Fturb(exp(-pow4(0.25*RT)));

        const volInternalScalarField Egamma
        (
            alpha()*rho()*ca2_*Omega*Fturb*gammaInt_()
        );

        // Intermittency equation
        tmp<fvScalarMatrix> gammaIntEqn
        (
            fvm::ddt(alpha, rho, gammaInt_)
          + fvm::div(alphaRhoPhi, gammaInt_)
          - fvm::laplacian(alpha*rho*DgammaIntEff(), gammaInt_)
        ==
            Pgamma - fvm::Sp(ce1_*Pgamma, gammaInt_)
          + Egamma - fvm::Sp(ce2_*Egamma, gammaInt_)
          + fvModels.source(alpha, rho, gammaInt_)
        );

        gammaIntEqn.ref().relax();
        fvConstraints.constrain(gammaIntEqn.ref());
        solve(gammaIntEqn);
        fvConstraints.constrain(gammaInt_);
        bound(gammaInt_, 0);
    }

    const volInternalScalarField Freattach(exp(-pow4(RT/20.0)));
    const volInternalScalarField gammaSep
    (
        min(2*max(Rev/(3.235*ReThetac) - 1, scalar(0))*Freattach, scalar(2))
       *Fthetat
    );

    gammaIntEff_ = max(gammaInt_(), gammaSep);
}


template<class BasicMomentumTransportModel>
void kOmegaSSTLM<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Correct ReThetat and gammaInt
    correctReThetatGammaInt();

    // Correct k and omega
    kOmegaSST<BasicMomentumTransportModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
