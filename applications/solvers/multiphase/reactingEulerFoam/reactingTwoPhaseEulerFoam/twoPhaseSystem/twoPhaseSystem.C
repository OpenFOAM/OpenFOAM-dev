/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

#include "twoPhaseSystem.H"
#include "dragModel.H"
#include "virtualMassModel.H"

#include "MULES.H"
#include "subCycle.H"
#include "UniformField.H"

#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "fvcSup.H"

#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseSystem, 0);
    defineRunTimeSelectionTable(twoPhaseSystem, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::twoPhaseSystem
(
    const fvMesh& mesh
)
:
    phaseSystem(mesh),
    phase1_(phaseModels_[0]),
    phase2_(phaseModels_[1])
{
    phase2_.volScalarField::operator=(scalar(1) - phase1_);

    volScalarField& alpha1 = phase1_;
    mesh.setFluxRequired(alpha1.name());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::twoPhaseSystem::~twoPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::sigma() const
{
    return sigma
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::Kd() const
{
    return Kd
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Kdf() const
{
    return Kdf
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::Vm() const
{
    return Vm
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


void Foam::twoPhaseSystem::solve
(
    const PtrList<volScalarField>& rAUs,
    const PtrList<surfaceScalarField>& rAUfs
)
{
    volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;

    const dictionary& alphaControls = mesh_.solverDict(alpha1.name());

    label nAlphaSubCycles(alphaControls.lookup<label>("nAlphaSubCycles"));
    label nAlphaCorr(alphaControls.lookup<label>("nAlphaCorr"));

    bool LTS = fv::localEulerDdt::enabled(mesh_);

    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    const surfaceScalarField& phi1 = phase1_.phi();
    const surfaceScalarField& phi2 = phase2_.phi();

    // Construct the dilatation rate source term
    tmp<volScalarField::Internal> tdgdt;

    if (phase1_.divU().valid() && phase2_.divU().valid())
    {
        tdgdt =
        (
            alpha2()
           *phase1_.divU()()()
          - alpha1()
           *phase2_.divU()()()
        );
    }
    else if (phase1_.divU().valid())
    {
        tdgdt =
        (
            alpha2()
           *phase1_.divU()()()
        );
    }
    else if (phase2_.divU().valid())
    {
        tdgdt =
        (
          - alpha1()
           *phase2_.divU()()()
        );
    }

    alpha1.correctBoundaryConditions();

    surfaceScalarField phir("phir", phi1 - phi2);

    tmp<surfaceScalarField> alphaPhiDbyA0;
    if (implicitPhasePressure() && (rAUs.size() || rAUfs.size()))
    {
        alphaPhiDbyA0 =
            this->DByAfs(rAUs, rAUfs)[phase1_.index()]
           *fvc::snGrad(alpha1, "bounded")*mesh_.magSf();
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimless/dimTime, 0)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                mesh_.time().timeName(),
                mesh_
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi_)*min(alpha1, scalar(1))
        );

        if (tdgdt.valid())
        {
            scalarField& dgdt = tdgdt.ref();

            forAll(dgdt, celli)
            {
                if (dgdt[celli] > 0.0)
                {
                    Sp[celli] -= dgdt[celli]/max(1 - alpha1[celli], 1e-4);
                    Su[celli] += dgdt[celli]/max(1 - alpha1[celli], 1e-4);
                }
                else if (dgdt[celli] < 0.0)
                {
                    Sp[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
                }
            }
        }

        tmp<volScalarField> trSubDeltaT;

        if (LTS && nAlphaSubCycles > 1)
        {
            trSubDeltaT =
                fv::localEulerDdt::localRSubDeltaT(mesh_, nAlphaSubCycles);
        }

        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
            surfaceScalarField alphaPhi1
            (
                fvc::flux
                (
                    phi_,
                    alpha1,
                    alphaScheme
                )
              + fvc::flux
                (
                    -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                    alpha1,
                    alpharScheme
                )
            );

            phase1_.correctInflowOutflow(alphaPhi1);

            if (alphaPhiDbyA0.valid())
            {
                alphaPhi1 +=
                    fvc::interpolate(max(alpha1, scalar(0)))
                   *fvc::interpolate(max(scalar(1) - alpha1, scalar(0)))
                   *alphaPhiDbyA0();
            }

            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi_,
                alphaPhi1,
                Sp,
                Su,
                UniformField<scalar>(phase1_.alphaMax()),
                zeroField()
            );

            if (alphaSubCycle.index() == 1)
            {
                phase1_.alphaPhiRef() = alphaPhi1;
            }
            else
            {
                phase1_.alphaPhiRef() += alphaPhi1;
            }

            if (alphaPhiDbyA0.valid())
            {
                const surfaceScalarField alphaDbyA
                (
                    fvc::interpolate(max(alpha1, scalar(0)))
                   *fvc::interpolate(max(scalar(1) - alpha1, scalar(0)))
                   *this->DByAfs(rAUs, rAUfs)[phase1_.index()]
                );

                fvScalarMatrix alpha1Eqn
                (
                    fvm::ddt(alpha1) - fvc::ddt(alpha1)
                  - fvm::laplacian(alphaDbyA, alpha1, "bounded")
                );

                alpha1Eqn.solve();

                phase1_.alphaPhiRef() += alpha1Eqn.flux();
            }
        }

        if (nAlphaSubCycles > 1)
        {
            phase1_.alphaPhiRef() /= nAlphaSubCycles;
        }

        phase1_.alphaRhoPhiRef() =
            fvc::interpolate(phase1_.rho())*phase1_.alphaPhi();

        phase2_.alphaPhiRef() = phi_ - phase1_.alphaPhi();
        phase2_.correctInflowOutflow(phase2_.alphaPhiRef());
        phase2_.alphaRhoPhiRef() =
            fvc::interpolate(phase2_.rho())*phase2_.alphaPhi();

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh_.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;

        // Ensure the phase-fractions are bounded
        alpha1.maxMin(0, 1);

        // Update the phase-fraction of the other phase
        alpha2 = scalar(1) - alpha1;
    }
}


// ************************************************************************* //
