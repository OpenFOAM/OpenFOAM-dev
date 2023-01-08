/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2023 OpenFOAM Foundation
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

#include "compressibleMultiphaseMixture.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "correctContactAngle.H"
#include "Time.H"
#include "subCycle.H"
#include "MULES.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(compressibleMultiphaseMixture, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::compressibleMultiphaseMixture::calcAlphas()
{
    scalar level = 0.0;
    alphas_ == 0.0;

    forAll(phases_, phasei)
    {
        alphas_ += level*phases_[phasei];
        level += 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::compressibleMultiphaseMixture::compressibleMultiphaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "phaseProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    mesh_(U.mesh()),

    p_
    (
        IOobject
        (
            "p",
            U.mesh().time().name(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    T_
    (
        IOobject
        (
            "T",
            U.mesh().time().name(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),

    phases_(lookup("phases"), compressiblePhase::iNew(T_)),

    rho_
    (
        IOobject
        (
            "rho",
            U.mesh().time().name(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar("rho", dimDensity, 0)
    ),

    U_(U),

    phi_(phi),

    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimTime, 0)
    ),

    alphas_
    (
        IOobject
        (
            "alphas",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
    ),

    sigmas_(lookup("sigmas")),
    dimSigma_(1, 0, -2, 0, 0),
    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh_.V()), 1.0/3.0)
    )
{
    calcAlphas();
    alphas_.write();
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::compressibleMultiphaseMixture::correctThermo()
{
    forAll(phases_, phasei)
    {
        phases_[phasei].correct(p_, T_);
    }
}


void Foam::compressibleMultiphaseMixture::correct()
{
    rho_ = phases_[0]*phases_[0].thermo().rho();

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        rho_ += phases_[phasei]*phases_[phasei].thermo().rho();
    }

    forAll(phases_, phasei)
    {
        phases_[phasei].Alpha() =
            phases_[phasei]*phases_[phasei].thermo().rho()/rho_;
    }
}


void Foam::compressibleMultiphaseMixture::correctRho(const volScalarField& dp)
{
    forAll(phases_, phasei)
    {
        phases_[phasei].thermo().rho() += phases_[phasei].thermo().psi()*dp;
    }
}


Foam::tmp<Foam::volScalarField> Foam::compressibleMultiphaseMixture::nu() const
{
    volScalarField mu(phases_[0].Alpha()*phases_[0].thermo().mu());

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        mu += phases_[phasei].Alpha()*phases_[phasei].thermo().mu();
    }

    return mu/rho_;
}


Foam::tmp<Foam::scalarField> Foam::compressibleMultiphaseMixture::nu
(
    const label patchi
) const
{
    scalarField mu
    (
        phases_[0].Alpha().boundaryField()[patchi]
       *phases_[0].thermo().mu(patchi)
    );

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        mu +=
            phases_[phasei].Alpha().boundaryField()[patchi]
           *phases_[phasei].thermo().mu(patchi);
    }

    return mu/rho_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::compressibleMultiphaseMixture::alphaEff
(
    const volScalarField& nut
) const
{
    tmp<volScalarField> talphaEff
    (
        phases_[0]
       *(
           phases_[0].thermo().kappa()
         + phases_[0].thermo().rho()*phases_[0].thermo().Cp()*nut
        )/phases_[0].thermo().Cv()
    );

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        talphaEff.ref() +=
            phases_[phasei]
           *(
               phases_[phasei].thermo().kappa()
             + phases_[phasei].thermo().rho()*phases_[phasei].thermo().Cp()*nut
            )/phases_[phasei].thermo().Cv();
    }

    return talphaEff;
}


Foam::tmp<Foam::volScalarField> Foam::compressibleMultiphaseMixture::rCv() const
{
    tmp<volScalarField> trCv(phases_[0]/phases_[0].thermo().Cv());

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        trCv.ref() += phases_[phasei]/phases_[phasei].thermo().Cv();
    }

    return trCv;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::compressibleMultiphaseMixture::surfaceTensionForce() const
{
    tmp<surfaceScalarField> tstf
    (
        surfaceScalarField::New
        (
            "surfaceTensionForce",
            mesh_,
            dimensionedScalar(dimensionSet(1, -2, -2, 0, 0), 0)
        )
    );

    surfaceScalarField& stf = tstf.ref();

    forAll(phases_, phasei)
    {
        const compressiblePhase& alpha1 = phases_[phasei];

        for (label phasej = phasei+1; phasej<phases_.size(); phasej++)
        {
            const compressiblePhase& alpha2 = phases_[phasej];

            sigmaTable::const_iterator sigma =
                sigmas_.find(interfacePair(alpha1, alpha2));

            if (sigma == sigmas_.end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of sigma values"
                    << exit(FatalError);
            }

            stf += dimensionedScalar(dimSigma_, sigma())
               *fvc::interpolate(K(alpha1, alpha2))*
                (
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }

    return tstf;
}


void Foam::compressibleMultiphaseMixture::solve()
{
    const Time& runTime = mesh_.time();

    const dictionary& alphaControls = mesh_.solution().solverDict("alpha");
    label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));
    scalar cAlpha(alphaControls.lookup<scalar>("cAlpha"));

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum(0.0*rhoPhi_);
        dimensionedScalar totalDeltaT = runTime.deltaT();

        List<volScalarField*> alphaPtrs(phases_.size());
        forAll(phases_, phasei)
        {
            alphaPtrs[phasei] = &phases_[phasei];
        }

        for
        (
            subCycle<volScalarField, subCycleFields> alphaSubCycle
            (
                alphaPtrs,
                nAlphaSubCycles
            );
            !(++alphaSubCycle).end();
        )
        {
            solveAlphas(cAlpha);
            rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi_;
        }

        rhoPhi_ = rhoPhiSum;
    }
    else
    {
        solveAlphas(cAlpha);
    }

    correct();
}


Foam::tmp<Foam::surfaceVectorField> Foam::compressibleMultiphaseMixture::nHatfv
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    /*
    // Cell gradient of alpha
    volVectorField gradAlpha =
        alpha2*fvc::grad(alpha1) - alpha1*fvc::grad(alpha2);

    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);
    */

    surfaceVectorField gradAlphaf
    (
        fvc::interpolate(alpha2)*fvc::interpolate(fvc::grad(alpha1))
      - fvc::interpolate(alpha1)*fvc::interpolate(fvc::grad(alpha2))
    );

    // Face unit interface normal
    return gradAlphaf/(mag(gradAlphaf) + deltaN_);
}


Foam::tmp<Foam::surfaceScalarField> Foam::compressibleMultiphaseMixture::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


Foam::tmp<Foam::volScalarField> Foam::compressibleMultiphaseMixture::K
(
    const phase& alpha1,
    const phase& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    correctContactAngle
    (
        alpha1,
        alpha2,
        U_.boundaryField(),
        deltaN_,
        tnHatfv.ref().boundaryFieldRef()
    );

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField>
Foam::compressibleMultiphaseMixture::nearInterface() const
{
    tmp<volScalarField> tnearInt
    (
        volScalarField::New
        (
            "nearInterface",
            mesh_,
            dimensionedScalar(dimless, 0)
        )
    );

    forAll(phases_, phasei)
    {
        tnearInt.ref() = max
        (
            tnearInt(),
            pos0(phases_[phasei] - 0.01)*pos0(0.99 - phases_[phasei])
        );
    }

    return tnearInt;
}


void Foam::compressibleMultiphaseMixture::solveAlphas
(
    const scalar cAlpha
)
{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi_/mesh_.magSf()));
    phic = min(cAlpha*phic, max(phic));

    UPtrList<const volScalarField> alphas(phases_.size());
    PtrList<surfaceScalarField> alphaPhis(phases_.size());

    forAll(phases_, phasei)
    {
        const compressiblePhase& alpha = phases_[phasei];

        alphas.set(phasei, &alpha);

        alphaPhis.set
        (
            phasei,
            new surfaceScalarField
            (
                phi_.name() + alpha.name(),
                fvc::flux
                (
                    phi_,
                    alpha,
                    alphaScheme
                )
            )
        );

        surfaceScalarField& alphaPhi = alphaPhis[phasei];

        forAll(phases_, phasej)
        {
            compressiblePhase& alpha2 = phases_[phasej];

            if (&alpha2 == &alpha) continue;

            surfaceScalarField phir(phic*nHatf(alpha, alpha2));

            alphaPhi += fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha,
                alpharScheme
            );
        }

        // Limit alphaPhi for each phase
        MULES::limit
        (
            1.0/mesh_.time().deltaT().value(),
            geometricOneField(),
            alpha,
            phi_,
            alphaPhi,
            zeroField(),
            zeroField(),
            oneField(),
            zeroField(),
            false
        );
    }

    MULES::limitSum(alphas, alphaPhis, phi_);

    rhoPhi_ = dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), 0);

    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh_.time().name(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, 0)
    );

    volScalarField divU(fvc::div(fvc::absolute(phi_, U_)));

    forAll(phases_, phasei)
    {
        compressiblePhase& alpha = phases_[phasei];

        surfaceScalarField& alphaPhi = alphaPhis[phasei];

        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                mesh_.time().name(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(alpha.dgdt().dimensions(), 0)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                mesh_.time().name(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            divU.v()*min(alpha.v(), scalar(1))
        );

        {
            const scalarField& dgdt = alpha.dgdt();

            forAll(dgdt, celli)
            {
                if (dgdt[celli] < 0.0 && alpha[celli] > 0.0)
                {
                    Sp[celli] += dgdt[celli]*alpha[celli];
                    Su[celli] -= dgdt[celli]*alpha[celli];
                }
                else if (dgdt[celli] > 0.0 && alpha[celli] < 1.0)
                {
                    Sp[celli] -= dgdt[celli]*(1.0 - alpha[celli]);
                }
            }
        }

        forAll(phases_, phasej)
        {
            const compressiblePhase& alpha2 = phases_[phasej];

            if (&alpha2 == &alpha) continue;

            const scalarField& dgdt2 = alpha2.dgdt();

            forAll(dgdt2, celli)
            {
                if (dgdt2[celli] > 0.0 && alpha2[celli] < 1.0)
                {
                    Sp[celli] -= dgdt2[celli]*(1.0 - alpha2[celli]);
                    Su[celli] += dgdt2[celli]*alpha[celli];
                }
                else if (dgdt2[celli] < 0.0 && alpha2[celli] > 0.0)
                {
                    Sp[celli] += dgdt2[celli]*alpha2[celli];
                }
            }
        }

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            alphaPhi,
            Sp,
            Su
        );

        rhoPhi_ += fvc::interpolate(alpha.thermo().rho())*alphaPhi;

        Info<< alpha.name() << " volume fraction, min, max = "
            << alpha.weightedAverage(mesh_.V()).value()
            << ' ' << min(alpha).value()
            << ' ' << max(alpha).value()
            << endl;

        sumAlpha += alpha;
    }

    Info<< "Phase-sum volume fraction, min, max = "
        << sumAlpha.weightedAverage(mesh_.V()).value()
        << ' ' << min(sumAlpha).value()
        << ' ' << max(sumAlpha).value()
        << endl;

    calcAlphas();
}


// ************************************************************************* //
