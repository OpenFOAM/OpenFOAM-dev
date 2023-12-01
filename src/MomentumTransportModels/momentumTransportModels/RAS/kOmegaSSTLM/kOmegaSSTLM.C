/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
tmp<volScalarField::Internal> kOmegaSSTLM<BasicMomentumTransportModel>::Pk
(
    const volScalarField::Internal& G
) const
{
    return gammaIntEff_*kOmegaSST<BasicMomentumTransportModel>::Pk(G);
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSSTLM<BasicMomentumTransportModel>::epsilonByk
(
    const volScalarField::Internal& F1,
    const volScalarField::Internal& F2
) const
{
    return
        min(max(gammaIntEff_, scalar(0.1)), scalar(1))
       *kOmegaSST<BasicMomentumTransportModel>::epsilonByk(F1, F2);
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal> kOmegaSSTLM<BasicMomentumTransportModel>::Fthetat
(
    const volScalarField::Internal& Us,
    const volScalarField::Internal& Omega,
    const volScalarField::Internal& nu
) const
{
    const volScalarField::Internal& omega = this->omega_();
    const volScalarField::Internal& y = this->y()();

    const volScalarField::Internal delta(375*Omega*nu*ReThetat_()*y/sqr(Us));
    const volScalarField::Internal ReOmega(sqr(y)*omega/nu);
    const volScalarField::Internal Fwake(exp(-sqr(ReOmega/1e5)));

    return volScalarField::Internal::New
    (
        this->groupName("Fthetat"),
        min
        (
            max
            (
                Fwake*exp(-pow4((y/delta))),
                (1 - sqr((gammaInt_() - 1.0/ce2_)/(1 - 1.0/ce2_)))
            ),
            scalar(1)
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
kOmegaSSTLM<BasicMomentumTransportModel>::ReThetac() const
{
    tmp<volScalarField::Internal> tReThetac
    (
        volScalarField::Internal::New
        (
            this->groupName("ReThetac"),
            this->mesh_,
            dimless
        )
    );
    volScalarField::Internal& ReThetac = tReThetac.ref();

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
tmp<volScalarField::Internal> kOmegaSSTLM<BasicMomentumTransportModel>::Flength
(
    const volScalarField::Internal& nu
) const
{
    tmp<volScalarField::Internal> tFlength
    (
        volScalarField::Internal::New
        (
            this->groupName("Flength"),
            this->mesh_,
            dimless
        )
    );
    volScalarField::Internal& Flength = tFlength.ref();

    const volScalarField::Internal& omega = this->omega_();
    const volScalarField::Internal& y = this->y()();

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
tmp<volScalarField::Internal>
kOmegaSSTLM<BasicMomentumTransportModel>::ReThetat0
(
    const volScalarField::Internal& Us,
    const volScalarField::Internal& dUsds,
    const volScalarField::Internal& nu
) const
{
    tmp<volScalarField::Internal> tReThetat0
    (
        volScalarField::Internal::New
        (
            this->groupName("ReThetat0"),
            this->mesh_,
            dimless
        )
    );
    volScalarField::Internal& ReThetat0 = tReThetat0.ref();

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
tmp<volScalarField::Internal> kOmegaSSTLM<BasicMomentumTransportModel>::Fonset
(
    const volScalarField::Internal& Rev,
    const volScalarField::Internal& ReThetac,
    const volScalarField::Internal& RT
) const
{
    const volScalarField::Internal Fonset1(Rev/(2.193*ReThetac));

    const volScalarField::Internal Fonset2
    (
        min(max(Fonset1, pow4(Fonset1)), scalar(2))
    );

    const volScalarField::Internal Fonset3(max(1 - pow3(RT/2.5), scalar(0)));

    return volScalarField::Internal::New
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

    ca1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ca1",
            this->coeffDict_,
            2
        )
    ),
    ca2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ca2",
            this->coeffDict_,
            0.06
        )
    ),
    ce1_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ce1",
            this->coeffDict_,
            1
        )
    ),
    ce2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ce2",
            this->coeffDict_,
            50
        )
    ),
    cThetat_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "cThetat",
            this->coeffDict_,
            0.03
        )
    ),
    sigmaThetat_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "sigmaThetat",
            this->coeffDict_,
            2
        )
    ),
    lambdaErr_
    (
        this->coeffDict_.lookupOrDefault("lambdaErr", 1e-6)
    ),
    maxLambdaIter_
    (
        this->coeffDict_.lookupOrDefault("maxLambdaIter", 10)
    ),
    deltaU_("deltaU", dimVelocity, small),

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
        this->mesh_
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
        this->mesh_
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
        ca1_.readIfPresent(this->coeffDict());
        ca2_.readIfPresent(this->coeffDict());
        ce1_.readIfPresent(this->coeffDict());
        ce2_.readIfPresent(this->coeffDict());
        sigmaThetat_.readIfPresent(this->coeffDict());
        cThetat_.readIfPresent(this->coeffDict());
        this->coeffDict().readIfPresent("lambdaErr", lambdaErr_);
        this->coeffDict().readIfPresent("maxLambdaIter", maxLambdaIter_);

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
    const volScalarField::Internal& nu = tnu()();
    const volScalarField::Internal& y = this->y()();
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    // Fields derived from the velocity gradient
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal Omega(sqrt(2*magSqr(skew(tgradU()()))));
    const volScalarField::Internal S(sqrt(2*magSqr(symm(tgradU()()))));
    const volScalarField::Internal Us(max(mag(U()), deltaU_));
    const volScalarField::Internal dUsds((U() & (U() & tgradU()()))/sqr(Us));
    tgradU.clear();

    const volScalarField::Internal Fthetat(this->Fthetat(Us, Omega, nu));

    {
        const volScalarField::Internal t(500*nu/sqr(Us));
        const volScalarField::Internal Pthetat
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

    const volScalarField::Internal ReThetac(this->ReThetac());
    const volScalarField::Internal Rev(sqr(y)*S/nu);
    const volScalarField::Internal RT(k()/(nu*omega()));

    {
        const volScalarField::Internal Pgamma
        (
            alpha()*rho()
           *ca1_*Flength(nu)*S*sqrt(gammaInt_()*Fonset(Rev, ReThetac, RT))
        );

        const volScalarField::Internal Fturb(exp(-pow4(0.25*RT)));

        const volScalarField::Internal Egamma
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

    const volScalarField::Internal Freattach(exp(-pow4(RT/20.0)));
    const volScalarField::Internal gammaSep
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
