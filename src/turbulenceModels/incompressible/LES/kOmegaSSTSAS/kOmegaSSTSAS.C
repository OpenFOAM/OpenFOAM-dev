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

#include "kOmegaSSTSAS.H"
#include "addToRunTimeSelectionTable.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kOmegaSSTSAS, 0);
addToRunTimeSelectionTable(LESModel, kOmegaSSTSAS, dictionary);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void kOmegaSSTSAS::updateSubGridScaleFields(const volScalarField& S2)
{
    nuSgs_ == a1_*k_/max(a1_*omega_, F2()*sqrt(S2));
    nuSgs_.correctBoundaryConditions();
}


tmp<volScalarField> kOmegaSSTSAS::F1(const volScalarField& CDkOmega) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}


tmp<volScalarField> kOmegaSSTSAS::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


tmp<volScalarField> kOmegaSSTSAS::Lvk2
(
    const volScalarField& S2
) const
{
    return max
    (
        kappa_*sqrt(S2)
       /(
            mag(fvc::laplacian(U()))
          + dimensionedScalar
            (
                "ROOTVSMALL",
                dimensionSet(0, -1 , -1, 0, 0, 0, 0),
                ROOTVSMALL
            )
        ),
        Cs_*delta()
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaSSTSAS::kOmegaSSTSAS
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    LESModel(modelName, U, phi, transport, turbulenceModelName),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85034
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.85616
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.4403
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),
    Cs_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cs",
            coeffDict_,
            0.262
        )
    ),
    alphaPhi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaPhi",
            coeffDict_,
            0.666667
        )
    ),
    zetaTilda2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "zetaTilda2",
            coeffDict_,
            1.755
        )
    ),
    FSAS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "FSAS",
            coeffDict_,
            1.25
        )
    ),

    omegaMin_("omegaMin", dimless/dimTime, SMALL),
    y_(mesh_),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
         )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            *this,
            0.41
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

    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    nuSgs_
    (
        IOobject
        (
            "nuSgs",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    )
{
    omegaMin_.readIfPresent(*this);

    bound(k_, kMin_);
    bound(omega_, omegaMin_);

    updateSubGridScaleFields(2.0*magSqr(symm(fvc::grad(U))));

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kOmegaSSTSAS::correct(const tmp<volTensorField>& gradU)
{
    LESModel::correct(gradU);

    if (mesh_.changing())
    {
        y_.correct();
    }

    volScalarField S2(2.0*magSqr(symm(gradU())));
    gradU.clear();

    volVectorField gradK(fvc::grad(k_));
    volVectorField gradOmega(fvc::grad(omega_));
    volScalarField L(sqrt(k_)/(pow025(Cmu_)*omega_));
    volScalarField CDkOmega((2.0*alphaOmega2_)*(gradK & gradOmega)/omega_);
    volScalarField F1(this->F1(CDkOmega));
    volScalarField G(GName(), nuSgs_*S2);

    // Turbulent kinetic energy equation
    {
        fvScalarMatrix kEqn
        (
            fvm::ddt(k_)
          + fvm::div(phi(), k_)
          - fvm::laplacian(DkEff(F1), k_)
        ==
            min(G, c1_*betaStar_*k_*omega_)
          - fvm::Sp(betaStar_*omega_, k_)
        );

        kEqn.relax();
        kEqn.solve();
    }
    bound(k_, kMin_);

    tmp<volScalarField> grad_omega_k = max
    (
        magSqr(gradOmega)/sqr(omega_),
        magSqr(gradK)/sqr(k_)
    );

    // Turbulent frequency equation
    {
        fvScalarMatrix omegaEqn
        (
            fvm::ddt(omega_)
          + fvm::div(phi(), omega_)
          - fvm::laplacian(DomegaEff(F1), omega_)
        ==
            gamma(F1)*S2
          - fvm::Sp(beta(F1)*omega_, omega_)
          - fvm::SuSp       // cross diffusion term
            (
                (F1 - scalar(1))*CDkOmega/omega_,
                omega_
            )
          + FSAS_
           *max
            (
                dimensionedScalar("zero",dimensionSet(0, 0, -2, 0, 0), 0.0),
                zetaTilda2_*kappa_*S2*sqr(L/Lvk2(S2))
              - 2.0/alphaPhi_*k_*grad_omega_k
            )
        );

        omegaEqn.relax();
        omegaEqn.solve();
    }
    bound(omega_, omegaMin_);

    updateSubGridScaleFields(S2);
}


tmp<volScalarField> kOmegaSSTSAS::epsilon() const
{
    return 2.0*nuEff()*magSqr(symm(fvc::grad(U())));
}


tmp<volSymmTensorField> kOmegaSSTSAS::B() const
{
    return ((2.0/3.0)*I)*k() - nuSgs()*twoSymm(fvc::grad(U()));
}


tmp<volSymmTensorField> kOmegaSSTSAS::devReff() const
{
    return -nuEff()*dev(twoSymm(fvc::grad(U())));
}


tmp<fvVectorMatrix> kOmegaSSTSAS::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}


tmp<fvVectorMatrix> kOmegaSSTSAS::divDevRhoReff
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


bool kOmegaSSTSAS::read()
{
    if (LESModel::read())
    {
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        Cs_.readIfPresent(coeffDict());
        alphaPhi_.readIfPresent(coeffDict());
        zetaTilda2_.readIfPresent(coeffDict());
        FSAS_.readIfPresent(coeffDict());

        omegaMin_.readIfPresent(*this);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
