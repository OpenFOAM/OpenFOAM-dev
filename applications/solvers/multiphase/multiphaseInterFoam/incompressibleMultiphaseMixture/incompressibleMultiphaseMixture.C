/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "incompressibleMultiphaseMixture.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "correctContactAngle.H"
#include "Time.H"
#include "subCycle.H"
#include "MULES.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcFlux.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::incompressibleMultiphaseMixture::calcAlphas()
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

Foam::incompressibleMultiphaseMixture::incompressibleMultiphaseMixture
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

    phases_(lookup("phases"), incompressiblePhase::iNew(U.mesh())),

    mesh_(U.mesh()),
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

    nu_
    (
        IOobject
        (
            "nu",
            mesh_.time().name(),
            mesh_
        ),
        mu()/rho()
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
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::incompressibleMultiphaseMixture::rho() const
{
    tmp<volScalarField> trho
    (
        phases_[0]*phases_[0].rho()
    );
    volScalarField& rho = trho.ref();

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        rho += phases_[phasei]*phases_[phasei].rho();
    }

    return trho;
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleMultiphaseMixture::rho(const label patchi) const
{
    tmp<scalarField> trho
    (
        phases_[0].boundaryField()[patchi]*phases_[0].rho().value()
    );
    scalarField& rho = trho.ref();

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        rho +=
            phases_[phasei].boundaryField()[patchi]
           *phases_[phasei].rho().value();
    }

    return trho;
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleMultiphaseMixture::mu() const
{
    tmp<volScalarField> tmu
    (
        phases_[0]*phases_[0].rho()*phases_[0].nu()
    );
    volScalarField& mu = tmu.ref();

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        mu += phases_[phasei]*phases_[phasei].rho()*phases_[phasei].nu();
    }

    return tmu;
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleMultiphaseMixture::mu(const label patchi) const
{
    tmp<scalarField> tmu
    (
        phases_[0].boundaryField()[patchi]
       *phases_[0].rho().value()
       *phases_[0].nu(patchi)
    );
    scalarField& mu = tmu.ref();

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        mu +=
            phases_[phasei].boundaryField()[patchi]
           *phases_[phasei].rho().value()
           *phases_[phasei].nu(patchi);
    }

    return tmu;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleMultiphaseMixture::muf() const
{
    tmp<surfaceScalarField> tmuf
    (
        fvc::interpolate(phases_[0])
       *phases_[0].rho()*fvc::interpolate(phases_[0].nu())
    );
    surfaceScalarField& muf = tmuf.ref();

    for (label phasei=1; phasei<phases_.size(); phasei++)
    {
        muf +=
            fvc::interpolate(phases_[phasei])
           *phases_[phasei].rho()*fvc::interpolate(phases_[phasei].nu());
    }

    return tmuf;
}


Foam::tmp<Foam::volScalarField>
Foam::incompressibleMultiphaseMixture::nu() const
{
    return nu_;
}


Foam::tmp<Foam::scalarField>
Foam::incompressibleMultiphaseMixture::nu(const label patchi) const
{
    return nu_.boundaryField()[patchi];
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleMultiphaseMixture::nuf() const
{
    return muf()/fvc::interpolate(rho());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleMultiphaseMixture::surfaceTensionForce() const
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
        const incompressiblePhase& alpha1 = phases_[phasei];

        for (label phasej = phasei+1; phasej<phases_.size(); phasej++)
        {
            const incompressiblePhase& alpha2 = phases_[phasej];

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


void Foam::incompressibleMultiphaseMixture::solve()
{
    correct();

    const Time& runTime = mesh_.time();

    const dictionary& alphaControls = mesh_.solution().solverDict("alpha");
    label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));
    scalar cAlpha(alphaControls.lookup<scalar>("cAlpha"));

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.name(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(rhoPhi_.dimensions(), 0)
        );

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

    // Update the mixture kinematic viscosity
    nu_ = mu()/rho();
}


void Foam::incompressibleMultiphaseMixture::correct()
{
    forAll(phases_, phasei)
    {
        phases_[phasei].correct();
    }
}


Foam::tmp<Foam::surfaceVectorField>
Foam::incompressibleMultiphaseMixture::nHatfv
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


Foam::tmp<Foam::surfaceScalarField> Foam::incompressibleMultiphaseMixture::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


Foam::tmp<Foam::volScalarField> Foam::incompressibleMultiphaseMixture::K
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
Foam::incompressibleMultiphaseMixture::nearInterface() const
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


void Foam::incompressibleMultiphaseMixture::solveAlphas
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
        const incompressiblePhase& alpha = phases_[phasei];

        alphas.set(phasei, &alpha);

        alphaPhis.set
        (
            phasei,
            new surfaceScalarField
            (
                "phi" + alpha.name() + "Corr",
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
            incompressiblePhase& alpha2 = phases_[phasej];

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

    forAll(phases_, phasei)
    {
        incompressiblePhase& alpha = phases_[phasei];
        surfaceScalarField& alphaPhi = alphaPhis[phasei];

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            alphaPhi
        );

        rhoPhi_ += alphaPhi*alpha.rho();

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

    // Correct the sum of the phase-fractions to avoid 'drift'
    volScalarField sumCorr(1.0 - sumAlpha);
    forAll(phases_, phasei)
    {
        incompressiblePhase& alpha = phases_[phasei];
        alpha += alpha*sumCorr;
    }

    calcAlphas();
}


bool Foam::incompressibleMultiphaseMixture::read()
{
    if (regIOobject::read())
    {
        lookup("sigmas") >> sigmas_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
