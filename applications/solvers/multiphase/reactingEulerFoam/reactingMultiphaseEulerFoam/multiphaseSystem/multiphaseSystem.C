/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

#include "multiphaseSystem.H"
#include "alphaContactAngleFvPatchScalarField.H"

#include "MULES.H"
#include "subCycle.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcMeshPhi.H"
#include "fvcSup.H"

#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseSystem, 0);
    defineRunTimeSelectionTable(multiphaseSystem, dictionary);
}

const Foam::scalar Foam::multiphaseSystem::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::multiphaseSystem::calcAlphas()
{
    scalar level = 0.0;
    alphas_ == 0.0;

    forAll(phases(), i)
    {
        alphas_ += level*phases()[i];
        level += 1.0;
    }

    alphas_.correctBoundaryConditions();
}


void Foam::multiphaseSystem::solveAlphas()
{
    PtrList<surfaceScalarField> alphaPhiCorrs(phases().size());
    forAll(phases(), phasei)
    {
        phaseModel& phase = phases()[phasei];
        volScalarField& alpha1 = phase;

        phase.alphaPhi() =
            dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0);

        alphaPhiCorrs.set
        (
            phasei,
            new surfaceScalarField
            (
                "phi" + alpha1.name() + "Corr",
                fvc::flux
                (
                    phi_,
                    phase,
                    "div(phi," + alpha1.name() + ')'
                )
            )
        );

        surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];

        forAll(phases(), phasej)
        {
            phaseModel& phase2 = phases()[phasej];
            volScalarField& alpha2 = phase2;

            if (&phase2 == &phase) continue;

            surfaceScalarField phir(phase.phi() - phase2.phi());

            cAlphaTable::const_iterator cAlpha
            (
                cAlphas_.find(phasePairKey(phase.name(), phase2.name()))
            );

            if (cAlpha != cAlphas_.end())
            {
                surfaceScalarField phic
                (
                    (mag(phi_) + mag(phir))/mesh_.magSf()
                );

                phir += min(cAlpha()*phic, max(phic))*nHatf(phase, phase2);
            }

            word phirScheme
            (
                "div(phir," + alpha2.name() + ',' + alpha1.name() + ')'
            );

            alphaPhiCorr += fvc::flux
            (
                -fvc::flux(-phir, phase2, phirScheme),
                phase,
                phirScheme
            );
        }

        // Ensure that the flux at inflow BCs is preserved
        forAll(alphaPhiCorr.boundaryField(), patchi)
        {
            fvsPatchScalarField& alphaPhiCorrp =
                alphaPhiCorr.boundaryField()[patchi];

            if (!alphaPhiCorrp.coupled())
            {
                const scalarField& phi1p = phase.phi().boundaryField()[patchi];
                const scalarField& alpha1p = alpha1.boundaryField()[patchi];

                forAll(alphaPhiCorrp, facei)
                {
                    if (phi1p[facei] < 0)
                    {
                        alphaPhiCorrp[facei] = alpha1p[facei]*phi1p[facei];
                    }
                }
            }
        }

        MULES::limit
        (
            1.0/mesh_.time().deltaT().value(),
            geometricOneField(),
            phase,
            phi_,
            alphaPhiCorr,
            zeroField(),
            zeroField(),
            phase.alphaMax(),
            0,
            true
        );
    }

    MULES::limitSum(alphaPhiCorrs);

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


    volScalarField divU(fvc::div(fvc::absolute(phi_, phases().first().U())));

    forAll(phases(), phasei)
    {
        phaseModel& phase = phases()[phasei];
        volScalarField& alpha = phase;

        surfaceScalarField& alphaPhic = alphaPhiCorrs[phasei];
        alphaPhic += upwind<scalar>(mesh_, phi_).flux(phase);

        volScalarField::DimensionedInternalField Sp
        (
            IOobject
            (
                "Sp",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("Sp", divU.dimensions(), 0.0)
        );

        volScalarField::DimensionedInternalField Su
        (
            IOobject
            (
                "Su",
                mesh_.time().timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            divU*min(alpha, scalar(1))
        );

        if (phase.divU().valid())
        {
            const scalarField& dgdt = phase.divU()();

            forAll(dgdt, celli)
            {
                if (dgdt[celli] > 0.0)
                {
                    Sp[celli] -= dgdt[celli];
                    Su[celli] += dgdt[celli];
                }
                else if (dgdt[celli] < 0.0)
                {
                    Sp[celli] +=
                        dgdt[celli]
                       *(1.0 - alpha[celli])/max(alpha[celli], 1e-4);
                }
            }
        }

        forAll(phases(), phasej)
        {
            const phaseModel& phase2 = phases()[phasej];
            const volScalarField& alpha2 = phase2;

            if (&phase2 == &phase) continue;

            if (phase2.divU().valid())
            {
                const scalarField& dgdt2 = phase2.divU()();

                forAll(dgdt2, celli)
                {
                    if (dgdt2[celli] < 0.0)
                    {
                        Sp[celli] +=
                            dgdt2[celli]
                           *(1.0 - alpha2[celli])/max(alpha2[celli], 1e-4);

                        Su[celli] -=
                            dgdt2[celli]
                           *alpha[celli]/max(alpha2[celli], 1e-4);
                    }
                    else if (dgdt2[celli] > 0.0)
                    {
                        Sp[celli] -= dgdt2[celli];
                    }
                }
            }
        }

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            alphaPhic,
            Sp,
            Su
        );

        phase.alphaPhi() += alphaPhic;

        Info<< phase.name() << " volume fraction, min, max = "
            << phase.weightedAverage(mesh_.V()).value()
            << ' ' << min(phase).value()
            << ' ' << max(phase).value()
            << endl;

        sumAlpha += phase;
    }

    Info<< "Phase-sum volume fraction, min, max = "
        << sumAlpha.weightedAverage(mesh_.V()).value()
        << ' ' << min(sumAlpha).value()
        << ' ' << max(sumAlpha).value()
        << endl;
}


Foam::tmp<Foam::surfaceVectorField> Foam::multiphaseSystem::nHatfv
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


Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseSystem::nHatf
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

void Foam::multiphaseSystem::correctContactAngle
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    surfaceVectorField::GeometricBoundaryField& nHatb
) const
{
    const volScalarField::GeometricBoundaryField& gbf
        = phase1.boundaryField();

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
                    acap.thetaProps()
                   .find(phasePairKey(phase1.name(), phase2.name()));

            if (tp == acap.thetaProps().end())
            {
                FatalErrorIn
                (
                    "multiphaseSystem::correctContactAngle"
                    "(const phaseModel& phase1, const phaseModel& phase2, "
                    "fvPatchVectorFieldField& nHatb) const"
                )   << "Cannot find interface "
                    << phasePairKey(phase1.name(), phase2.name())
                    << "\n    in table of theta properties for patch "
                    << acap.patch().name()
                    << exit(FatalError);
            }

            bool matched = (tp.key().first() == phase1.name());

            scalar theta0 = convertToRad*tp().theta0(matched);
            scalarField theta(boundary[patchi].size(), theta0);

            scalar uTheta = tp().uTheta();

            // Calculate the dynamic contact angle if required
            if (uTheta > SMALL)
            {
                scalar thetaA = convertToRad*tp().thetaA(matched);
                scalar thetaR = convertToRad*tp().thetaR(matched);

                // Calculated the component of the velocity parallel to the wall
                vectorField Uwall
                (
                    phase1.U()().boundaryField()[patchi].patchInternalField()
                  - phase1.U()().boundaryField()[patchi]
                );
                Uwall -= (AfHatPatch & Uwall)*AfHatPatch;

                // Find the direction of the interface parallel to the wall
                vectorField nWall
                (
                    nHatPatch - (AfHatPatch & nHatPatch)*AfHatPatch
                );

                // Normalise nWall
                nWall /= (mag(nWall) + SMALL);

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


Foam::tmp<Foam::volScalarField> Foam::multiphaseSystem::K
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    tmp<surfaceVectorField> tnHatfv = nHatfv(phase1, phase2);

    correctContactAngle(phase1, phase2, tnHatfv().boundaryField());

    // Simple expression for curvature
    return -fvc::div(tnHatfv & mesh_.Sf());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseSystem::multiphaseSystem
(
    const fvMesh& mesh
)
:
    phaseSystem(mesh),

    alphas_
    (
        IOobject
        (
            "alphas",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alphas", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    cAlphas_(lookup("interfaceCompression")),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(mesh.V()), 1.0/3.0)
    )
{
    forAll(phases(), phasei)
    {
        volScalarField& alphai = phases()[phasei];
        mesh.setFluxRequired(alphai.name());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseSystem::~multiphaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::multiphaseSystem::surfaceTension
(
    const phaseModel& phase1
) const
{
    tmp<surfaceScalarField> tSurfaceTension
    (
        new surfaceScalarField
        (
            IOobject
            (
                "surfaceTension",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar
            (
                "surfaceTension",
                dimensionSet(1, -2, -2, 0, 0),
                0
            )
        )
    );

    forAll(phases(), phasej)
    {
        const phaseModel& phase2 = phases()[phasej];

        if (&phase2 != &phase1)
        {
            phasePairKey key12(phase1.name(), phase2.name());

            cAlphaTable::const_iterator cAlpha(cAlphas_.find(key12));

            if (cAlpha != cAlphas_.end())
            {
                tSurfaceTension() +=
                    fvc::interpolate(sigma(key12)*K(phase1, phase2))
                   *(
                        fvc::interpolate(phase2)*fvc::snGrad(phase1)
                      - fvc::interpolate(phase1)*fvc::snGrad(phase2)
                    );
            }
        }
    }

    return tSurfaceTension;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseSystem::nearInterface() const
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

    forAll(phases(), phasei)
    {
        tnearInt() = max
        (
            tnearInt(),
            pos(phases()[phasei] - 0.01)*pos(0.99 - phases()[phasei])
        );
    }

    return tnearInt;
}


void Foam::multiphaseSystem::solve()
{
    const fvMesh& mesh = this->mesh();
    const Time& runTime = mesh.time();

    const dictionary& alphaControls = mesh_.solverDict("alpha");
    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));

    bool LTS = fv::localEulerDdt::enabled(mesh);

    if (nAlphaSubCycles > 1)
    {
        tmp<volScalarField> trSubDeltaT;

        if (LTS)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
        }

        dimensionedScalar totalDeltaT = runTime.deltaT();

        PtrList<volScalarField> alpha0s(phases().size());
        PtrList<surfaceScalarField> phiSums(phases().size());

        forAll(phases(), phasei)
        {
            phaseModel& phase = phases()[phasei];
            volScalarField& alpha = phase;

            alpha0s.set
            (
                phasei,
                new volScalarField(alpha.oldTime())
            );

            phiSums.set
            (
                phasei,
                new surfaceScalarField
                (
                    IOobject
                    (
                        "phiSum" + alpha.name(),
                        runTime.timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("0", dimensionSet(0, 3, -1, 0, 0), 0)
                )
            );
        }

        for
        (
            subCycleTime alphaSubCycle
            (
                const_cast<Time&>(runTime),
                nAlphaSubCycles
            );
            !(++alphaSubCycle).end();
        )
        {
            solveAlphas();

            forAll(phases(), phasei)
            {
                phiSums[phasei] += phases()[phasei].phi();
            }
        }

        forAll(phases(), phasei)
        {
            phaseModel& phase = phases()[phasei];
            volScalarField& alpha = phase;

            phase.phi() = phiSums[phasei]/nAlphaSubCycles;

            // Correct the time index of the field
            // to correspond to the global time
            alpha.timeIndex() = runTime.timeIndex();

            // Reset the old-time field value
            alpha.oldTime() = alpha0s[phasei];
            alpha.oldTime().timeIndex() = runTime.timeIndex();
        }
    }
    else
    {
        solveAlphas();
    }

    forAll(phases(), phasei)
    {
        phaseModel& phase = phases()[phasei];
        phase.alphaRhoPhi() = fvc::interpolate(phase.rho())*phase.alphaPhi();
    }

    calcAlphas();
}


// ************************************************************************* //
