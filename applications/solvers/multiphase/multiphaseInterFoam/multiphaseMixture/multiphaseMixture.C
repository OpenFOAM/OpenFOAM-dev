/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "multiphaseMixture.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "Time.H"
#include "subCycle.H"
#include "MULES.H"
#include "surfaceInterpolate.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "fvcFlux.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::multiphaseMixture::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiphaseMixture::calcAlphas()
{
    scalar level = 0.0;
    alphas_ == 0.0;

    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        alphas_ += level*iter();
        level += 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseMixture::multiphaseMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    phases_(lookup("phases"), phase::iNew(U, phi)),

    mesh_(U.mesh()),
    U_(U),
    phi_(phi),

    rhoPhi_
    (
        IOobject
        (
            "rhoPhi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("rhoPhi", dimMass/dimTime, 0.0)
    ),

    alphas_
    (
        IOobject
        (
            "alphas",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("alphas", dimless, 0.0)
    ),

    nu_
    (
        IOobject
        (
            "nu",
            mesh_.time().timeName(),
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
Foam::multiphaseMixture::rho() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> trho = iter()*iter().rho();
    volScalarField& rho = trho.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        rho += iter()*iter().rho();
    }

    return trho;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::rho(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> trho = iter().boundaryField()[patchi]*iter().rho().value();
    scalarField& rho = trho.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        rho += iter().boundaryField()[patchi]*iter().rho().value();
    }

    return trho;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixture::mu() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<volScalarField> tmu = iter()*iter().rho()*iter().nu();
    volScalarField& mu = tmu.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        mu += iter()*iter().rho()*iter().nu();
    }

    return tmu;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::mu(const label patchi) const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<scalarField> tmu =
        iter().boundaryField()[patchi]
       *iter().rho().value()
       *iter().nu(patchi);
    scalarField& mu = tmu.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        mu +=
            iter().boundaryField()[patchi]
           *iter().rho().value()
           *iter().nu(patchi);
    }

    return tmu;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixture::muf() const
{
    PtrDictionary<phase>::const_iterator iter = phases_.begin();

    tmp<surfaceScalarField> tmuf =
        fvc::interpolate(iter())*iter().rho()*fvc::interpolate(iter().nu());
    surfaceScalarField& muf = tmuf.ref();

    for (++iter; iter != phases_.end(); ++iter)
    {
        muf +=
            fvc::interpolate(iter())*iter().rho()*fvc::interpolate(iter().nu());
    }

    return tmuf;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixture::nu() const
{
    return nu_;
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseMixture::nu(const label patchi) const
{
    return nu_.boundaryField()[patchi];
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixture::nuf() const
{
    return muf()/fvc::interpolate(rho());
}


Foam::tmp<Foam::surfaceScalarField>
Foam::multiphaseMixture::surfaceTensionForce() const
{
    tmp<surfaceScalarField> tstf
    (
        new surfaceScalarField
        (
            IOobject
            (
                "surfaceTensionForce",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "surfaceTensionForce",
                dimensionSet(1, -2, -2, 0, 0),
                0.0
            )
        )
    );

    surfaceScalarField& stf = tstf.ref();

    forAllConstIter(PtrDictionary<phase>, phases_, iter1)
    {
        const phase& alpha1 = iter1();

        PtrDictionary<phase>::const_iterator iter2 = iter1;
        ++iter2;

        for (; iter2 != phases_.end(); ++iter2)
        {
            const phase& alpha2 = iter2();

            sigmaTable::const_iterator sigma =
                sigmas_.find(interfacePair(alpha1, alpha2));

            if (sigma == sigmas_.end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << " in list of sigma values"
                    << exit(FatalError);
            }

            stf += dimensionedScalar("sigma", dimSigma_, sigma())
               *fvc::interpolate(K(alpha1, alpha2))*
                (
                    fvc::interpolate(alpha2)*fvc::snGrad(alpha1)
                  - fvc::interpolate(alpha1)*fvc::snGrad(alpha2)
                );
        }
    }

    return tstf;
}


void Foam::multiphaseMixture::solve()
{
    correct();

    const Time& runTime = mesh_.time();

    volScalarField& alpha = phases_.first();

    const dictionary& alphaControls = mesh_.solverDict("alpha");
    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    scalar cAlpha(readScalar(alphaControls.lookup("cAlpha")));

    if (nAlphaSubCycles > 1)
    {
        surfaceScalarField rhoPhiSum
        (
            IOobject
            (
                "rhoPhiSum",
                runTime.timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", rhoPhi_.dimensions(), 0)
        );

        dimensionedScalar totalDeltaT = runTime.deltaT();

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha, nAlphaSubCycles);
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


void Foam::multiphaseMixture::correct()
{
    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        iter().correct();
    }
}


Foam::tmp<Foam::surfaceVectorField> Foam::multiphaseMixture::nHatfv
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


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseMixture::nHatf
(
    const volScalarField& alpha1,
    const volScalarField& alpha2
) const
{
    // Face unit interface normal flux
    return nHatfv(alpha1, alpha2) & mesh_.Sf();
}


// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::multiphaseMixture::correctContactAngle
(
    const phase& alpha1,
    const phase& alpha2,
    surfaceVectorField::Boundary& nHatb
) const
{
    const volScalarField::Boundary& gbf
        = alpha1.boundaryField();

    const fvBoundaryMesh& boundary = mesh_.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(gbf[patchi]))
        {
            const alphaContactAngleFvPatchScalarField& acap =
                refCast<const alphaContactAngleFvPatchScalarField>(gbf[patchi]);

            vectorField& nHatPatch = nHatb[patchi];

            vectorField AfHatPatch
            (
                mesh_.Sf().boundaryField()[patchi]
               /mesh_.magSf().boundaryField()[patchi]
            );

            alphaContactAngleFvPatchScalarField::thetaPropsTable::
                const_iterator tp =
                acap.thetaProps().find(interfacePair(alpha1, alpha2));

            if (tp == acap.thetaProps().end())
            {
                FatalErrorInFunction
                    << "Cannot find interface " << interfacePair(alpha1, alpha2)
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            bool matched = (tp.key().first() == alpha1.name());

            scalar theta0 = convertToRad*tp().theta0(matched);
            scalarField theta(boundary[patchi].size(), theta0);

            scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > small)
            {
                scalar thetaA = convertToRad*tp().thetaA(matched);
                scalar thetaR = convertToRad*tp().thetaR(matched);

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    U_.boundaryField()[patchi].patchInternalField()
                  - U_.boundaryField()[patchi]
                );
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall
                (
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch
                );

                // Normalise nWall
                nWall /= (mag(nWall) + small);

                // Calculate Uwall resolved normal to the interface parallel to
                // the interface
                scalarField uwall(nWall & Uwall);

                theta += (thetaA - thetaR)*tanh(uwall/uTheta);
            }


            // Reset nHatPatch to correspond to the contact angle

            scalarField a12(nHatPatch & AfHatPatch);

            scalarField b1(cos(theta));

            scalarField b2(nHatPatch.size());

            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatPatch = a*AfHatPatch + b*nHatPatch;

            nHatPatch /= (mag(nHatPatch) + deltaN_.value());
        }
    }
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseMixture::K
(
    const phase& alpha1,
    const phase& alpha2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(alpha1, alpha2);

    correctContactAngle(alpha1, alpha2, tnHatfv.ref().boundaryFieldRef());

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseMixture::nearInterface() const
{
    tmp<volScalarField> tnearInt
    (
        new volScalarField
        (
            IOobject
            (
                "nearInterface",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("nearInterface", dimless, 0.0)
        )
    );

    forAllConstIter(PtrDictionary<phase>, phases_, iter)
    {
        tnearInt.ref() =
            max(tnearInt(), pos0(iter() - 0.01)*pos0(0.99 - iter()));
    }

    return tnearInt;
}


void Foam::multiphaseMixture::solveAlphas
(
    const scalar cAlpha
)
{
    static label nSolves=-1;
    nSolves++;

    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic(mag(phi_/mesh_.magSf()));
    phic = min(cAlpha*phic, max(phic));

    PtrList<surfaceScalarField> alphaPhiCorrs(phases_.size());
    int phasei = 0;

    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        phase& alpha = iter();

        alphaPhiCorrs.set
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

        surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];

        forAllIter(PtrDictionary<phase>, phases_, iter2)
        {
            phase& alpha2 = iter2();

            if (&alpha2 == &alpha) continue;

            surfaceScalarField phir(phic*nHatf(alpha, alpha2));

            alphaPhiCorr += fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha,
                alpharScheme
            );
        }

        MULES::limit
        (
            1.0/mesh_.time().deltaT().value(),
            geometricOneField(),
            alpha,
            phi_,
            alphaPhiCorr,
            zeroField(),
            zeroField(),
            oneField(),
            zeroField(),
            true
        );

        phasei++;
    }

    MULES::limitSum(alphaPhiCorrs);

    rhoPhi_ = dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0);

    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("sumAlpha", dimless, 0)
    );

    phasei = 0;

    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        phase& alpha = iter();

        surfaceScalarField& alphaPhi = alphaPhiCorrs[phasei];
        alphaPhi += upwind<scalar>(mesh_, phi_).flux(alpha);

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

        phasei++;
    }

    Info<< "Phase-sum volume fraction, min, max = "
        << sumAlpha.weightedAverage(mesh_.V()).value()
        << ' ' << min(sumAlpha).value()
        << ' ' << max(sumAlpha).value()
        << endl;

    // Correct the sum of the phase-fractions to avoid 'drift'
    volScalarField sumCorr(1.0 - sumAlpha);
    forAllIter(PtrDictionary<phase>, phases_, iter)
    {
        phase& alpha = iter();
        alpha += alpha*sumCorr;
    }

    calcAlphas();
}


bool Foam::multiphaseMixture::read()
{
    if (transportModel::read())
    {
        bool readOK = true;

        PtrList<entry> phaseData(lookup("phases"));
        label phasei = 0;

        forAllIter(PtrDictionary<phase>, phases_, iter)
        {
            readOK &= iter().read(phaseData[phasei++].dict());
        }

        lookup("sigmas") >> sigmas_;

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
