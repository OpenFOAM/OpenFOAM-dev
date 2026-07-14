/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "kkLOmega.H"
#include "bound.H"
#include "wallDist.H"
#include "makeMomentumTransportModel.H"

makeMomentumTransportModelTypes
(
    geometricOneField,
    geometricOneField,
    incompressibleMomentumTransportModel
)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kkLOmega, 0);
addToRunTimeSelectionTable
(
    RASincompressibleMomentumTransportModel,
    kkLOmega,
    dictionary
);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> kkLOmega::fv(const volScalarField& Ret) const
{
    return(1 - exp(-sqrt(Ret)/Av_));
}


tmp<volScalarField> kkLOmega::fINT() const
{
    return
    (
        min
        (
            kt_/(Cint_*(kl_ + kt_ + kMin_)),
            dimensionedScalar(dimless, 1)
        )
    );
}


tmp<volScalarField> kkLOmega::fSS(const volScalarField& Omega) const
{
    return(exp(-sqr(Css_*nu()*Omega/(kt_ + kMin_))));
}


tmp<volScalarField> kkLOmega::Cmu(const volScalarField& S) const
{
    return(1/(A0_ + As_*(S/(omega_ + omegaMin_))));
}


tmp<volScalarField> kkLOmega::BetaTS(const volScalarField& ReOmega) const
{
    return(scalar(1) - exp(-sqr(max(ReOmega - CtsCrit_, scalar(0)))/Ats_));
}


tmp<volScalarField> kkLOmega::fTaul
(
    const volScalarField& lambdaEff,
    const volScalarField& ktL,
    const volScalarField& Omega
) const
{
    return
    (
        scalar(1)
      - exp
        (
           -CtauL_*ktL
           /(
                sqr(lambdaEff*Omega)
              + dimensionedScalar
                (
                    "rootVSmall",
                    sqr(dimensions::velocity),
                    rootVSmall
                )
            )
        )
    );
}


tmp<volScalarField> kkLOmega::alphaT
(
    const volScalarField& lambdaEff,
    const volScalarField& fv,
    const volScalarField& ktS
) const
{
    return(fv*CmuStd_*sqrt(ktS)*lambdaEff);
}


tmp<volScalarField> kkLOmega::fOmega
(
    const volScalarField& lambdaEff,
    const volScalarField& lambdaT
) const
{
    return
    (
        scalar(1)
      - exp
        (
           -0.41
           *pow4
            (
                lambdaEff
              / (
                    lambdaT
                  + dimensionedScalar
                    (
                        "ROTvSmall",
                        lambdaT.dimensions(),
                        rootVSmall
                    )
                )
            )
        )
    );
}


tmp<volScalarField> kkLOmega::phiBP(const volScalarField& Omega) const
{
    return
    (
        min
        (
            max
            (
                kt_/nu()
             / (
                    Omega
                  + dimensionedScalar
                    (
                        "ROTvSmall",
                        Omega.dimensions(),
                        rootVSmall
                    )
                )
              - CbpCrit_,
                scalar(0)
            ),
            scalar(50)
        )
    );
}


tmp<volScalarField> kkLOmega::phiNAT
(
    const volScalarField& ReOmega,
    const volScalarField& fNatCrit
) const
{
    return
    (
        max
        (
            ReOmega
          - CnatCrit_
          / (
                fNatCrit + dimensionedScalar(dimless, rootVSmall)
            ),
            scalar(0)
        )
    );
}


tmp<volScalarField> kkLOmega::D(const volScalarField& k) const
{
    return nu()*magSqr(fvc::grad(sqrt(k)));
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void kkLOmega::correctNut()
{
    // Currently this function is not implemented due to the complexity of
    // evaluating nut.  Better calculate nut at the end of correct()
    NotImplemented;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kkLOmega::kkLOmega
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

    A0_("A0", typeDict(type), 4.04),
    As_("As", typeDict(type), 2.12),
    Av_("Av", typeDict(type), 6.75),
    Abp_("Abp", typeDict(type), 0.6),
    Anat_("Anat", typeDict(type), 200),
    Ats_("Ats", typeDict(type), 200),
    CbpCrit_("CbpCrit", typeDict(type), 1.2),
    Cnc_("Cnc", typeDict(type), 0.1),
    CnatCrit_("CnatCrit", typeDict(type), 1250),
    Cint_("Cint", typeDict(type), 0.75),
    CtsCrit_("CtsCrit", typeDict(type), 1000),
    CrNat_("CrNat", typeDict(type), 0.02),
    C11_("C11", typeDict(type), 3.4e-6),
    C12_("C12", typeDict(type), 1.0e-10),
    CR_("CR", typeDict(type), 0.12),
    CalphaTheta_("CalphaTheta", typeDict(type), 0.035),
    Css_("Css", typeDict(type), 1.5),
    CtauL_("CtauL", typeDict(type), 4360),
    Cw1_("Cw1", typeDict(type), 0.44),
    Cw2_("Cw2", typeDict(type), 0.92),
    Cw3_("Cw3", typeDict(type), 0.3),
    CwR_("CwR", typeDict(type), 1.5),
    Clambda_("Clambda", typeDict(type), 2.495),
    CmuStd_("CmuStd", typeDict(type), 0.09),
    Prtheta_("Prtheta", typeDict(type), 0.85),
    Sigmak_("Sigmak", typeDict(type), 1),
    Sigmaw_("Sigmaw", typeDict(type), 1.17),
    omegaMin_("omegaMin", dimensions::turbulentOmega, typeDict(type), small),
    kt_
    (
        IOobject
        (
            this->groupName("kt"),
            runTime_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensions::turbulentKineticEnergy
    ),
    kl_
    (
        IOobject
        (
            this->groupName("kl"),
            runTime_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensions::turbulentKineticEnergy
    ),
    omega_
    (
        IOobject
        (
            this->groupName("omega"),
            runTime_.name(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensions::turbulentOmega
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.name(),
            mesh_
        ),
        kt_*omega_
    )
{
    bound(kt_, kMin_);
    bound(kl_, kMin_);
    bound(omega_, omegaMin_);
    epsilon_ = kt_*omega_ + D(kl_) + D(kt_);

    if (type == typeName)
    {
        // Evaluating nut_ is complex so start from the field read from file
        nut_.correctBoundaryConditions();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool kkLOmega::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {
        A0_.readIfPresent(typeDict());
        As_.readIfPresent(typeDict());
        Av_.readIfPresent(typeDict());
        Abp_.readIfPresent(typeDict());
        Anat_.readIfPresent(typeDict());
        Abp_.readIfPresent(typeDict());
        Ats_.readIfPresent(typeDict());
        CbpCrit_.readIfPresent(typeDict());
        Cnc_.readIfPresent(typeDict());
        CnatCrit_.readIfPresent(typeDict());
        Cint_.readIfPresent(typeDict());
        CtsCrit_.readIfPresent(typeDict());
        CrNat_.readIfPresent(typeDict());
        C11_.readIfPresent(typeDict());
        C12_.readIfPresent(typeDict());
        CR_.readIfPresent(typeDict());
        CalphaTheta_.readIfPresent(typeDict());
        Css_.readIfPresent(typeDict());
        CtauL_.readIfPresent(typeDict());
        Cw1_.readIfPresent(typeDict());
        Cw2_.readIfPresent(typeDict());
        Cw3_.readIfPresent(typeDict());
        CwR_.readIfPresent(typeDict());
        Clambda_.readIfPresent(typeDict());
        CmuStd_.readIfPresent(typeDict());
        Prtheta_.readIfPresent(typeDict());
        Sigmak_.readIfPresent(typeDict());
        Sigmaw_.readIfPresent(typeDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kkLOmega::validate()
{}


void kkLOmega::correct()
{
    eddyViscosity<incompressible::RASModel>::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField y(this->y());

    const volScalarField lambdaT(sqrt(kt_)/(omega_ + omegaMin_));

    const volScalarField lambdaEff(min(Clambda_*y, lambdaT));

    const volScalarField fw
    (
        pow
        (
            lambdaEff
           /(lambdaT + dimensionedScalar(dimensions::length, rootVSmall)),
            2.0/3.0
        )
    );

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU = tgradU();

    const volScalarField Omega(sqrt(2.0)*mag(skew(gradU)));

    const volScalarField S2(2*magSqr(dev(symm(gradU))));

    const volScalarField ktS(fSS(Omega)*fw*kt_);

    const volScalarField nuts
    (
        fv(sqr(fw)*kt_/nu()/(omega_ + omegaMin_))
       *fINT()
       *Cmu(sqrt(S2))*sqrt(ktS)*lambdaEff
    );
    const volScalarField Pkt(this->GName(), nuts*S2);

    const volScalarField ktL(kt_ - ktS);
    const volScalarField ReOmega(sqr(y)*Omega/nu());

    const volScalarField nutl
    (
        min
        (
            C11_*fTaul(lambdaEff, ktL, Omega)*Omega*sqr(lambdaEff)
           *sqrt(ktL)*lambdaEff/nu()
          + C12_*BetaTS(ReOmega)*sqr(sqr(lambdaEff/Clambda_)*Omega)/nu()
          , 0.5*(kl_ + ktL)/(sqrt(S2) + omegaMin_)
        )
    );

    const volScalarField Pkl(nutl*S2);

    const volScalarField alphaTEff
    (
        alphaT(lambdaEff, fv(sqr(fw)*kt_/nu()/(omega_ + omegaMin_)), ktS)
    );

    // By pass source term divided by kl_

    const dimensionedScalar fwMin("small", dimless, rootVSmall);

    const volScalarField Rbp
    (
        CR_*(1 - exp(-phiBP(Omega)()/Abp_))*omega_
       /(fw + fwMin)
    );

    const volScalarField fNatCrit(1 - exp(-Cnc_*sqrt(kl_)*y/nu()));

    // Natural source term divided by kl_
    const volScalarField Rnat
    (
        CrNat_*(1 - exp(-phiNAT(ReOmega, fNatCrit)/Anat_))*Omega
    );


    omega_.boundaryFieldRef().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DomegaEff(alphaTEff), omega_)
     ==
        Cw1_*Pkt*omega_/(kt_ + kMin_)
      - fvm::SuSp
        (
            (1 - CwR_/(fw + fwMin))*kl_*(Rbp + Rnat)/(kt_ + kMin_)
          , omega_
        )
      - fvm::Sp(Cw2_*sqr(fw)*omega_, omega_)
      + (
            Cw3_*fOmega(lambdaEff, lambdaT)*alphaTEff*sqr(fw)*sqrt(kt_)
        )()()/pow3(y())
    );

    omegaEqn.ref().relax();
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());

    solve(omegaEqn);
    bound(omega_, omegaMin_);


    const volScalarField Dl(D(kl_));

    // Laminar kinetic energy equation
    tmp<fvScalarMatrix> klEqn
    (
        fvm::ddt(kl_)
      + fvm::div(phi_, kl_)
      - fvm::laplacian(nu(), kl_)
     ==
        Pkl
      - fvm::Sp(Rbp + Rnat + Dl/(kl_ + kMin_), kl_)
    );

    klEqn.ref().relax();
    klEqn.ref().boundaryManipulate(kl_.boundaryFieldRef());

    solve(klEqn);
    bound(kl_, kMin_);


    const volScalarField Dt(D(kt_));

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> ktEqn
    (
        fvm::ddt(kt_)
      + fvm::div(phi_, kt_)
      - fvm::laplacian(DkEff(alphaTEff), kt_)
     ==
        Pkt
      + (Rbp + Rnat)*kl_
      - fvm::Sp(omega_ + Dt/(kt_+ kMin_), kt_)
    );

    ktEqn.ref().relax();
    ktEqn.ref().boundaryManipulate(kt_.boundaryFieldRef());

    solve(ktEqn);
    bound(kt_, kMin_);


    // Update total fluctuation kinetic energy dissipation rate
    epsilon_ = kt_*omega_ + Dl + Dt;


    // Re-calculate turbulent viscosity
    nut_ = nuts + nutl;
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
