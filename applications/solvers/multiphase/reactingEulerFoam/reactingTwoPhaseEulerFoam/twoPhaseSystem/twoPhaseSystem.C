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

#include "twoPhaseSystem.H"
#include "dragModel.H"
#include "virtualMassModel.H"

#include "MULES.H"
#include "subCycle.H"

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


const Foam::dragModel& Foam::twoPhaseSystem::drag(const phaseModel& phase) const
{
    return lookupSubModel<dragModel>(phase, otherPhase(phase));
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


const Foam::virtualMassModel&
Foam::twoPhaseSystem::virtualMass(const phaseModel& phase) const
{
    return lookupSubModel<virtualMassModel>(phase, otherPhase(phase));
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::Vm() const
{
    return Vm
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Vmf() const
{
    return Vmf
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


Foam::tmp<Foam::volVectorField>
Foam::twoPhaseSystem::F() const
{
    return F
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::twoPhaseSystem::Ff() const
{
    return Ff
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::D() const
{
    return D
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


bool Foam::twoPhaseSystem::transfersMass() const
{
    return transfersMass(phase1());
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseSystem::dmdt() const
{
    return dmdt
    (
        phasePairKey(phase1().name(), phase2().name())
    );
}


void Foam::twoPhaseSystem::solve()
{
    const fvMesh& mesh = this->mesh();
    const Time& runTime = mesh.time();

    volScalarField& alpha1 = phase1_;
    volScalarField& alpha2 = phase2_;

    const dictionary& alphaControls = mesh.solverDict(alpha1.name());

    label nAlphaSubCycles(readLabel(alphaControls.lookup("nAlphaSubCycles")));
    label nAlphaCorr(readLabel(alphaControls.lookup("nAlphaCorr")));

    bool LTS = fv::localEulerDdt::enabled(mesh);

    word alphaScheme("div(phi," + alpha1.name() + ')');
    word alpharScheme("div(phir," + alpha1.name() + ')');

    const surfaceScalarField& phi = this->phi();
    const surfaceScalarField& phi1 = phase1_.phi();
    const surfaceScalarField& phi2 = phase2_.phi();

    // Construct the dilatation rate source term
    tmp<volScalarField::DimensionedInternalField> tdgdt;

    if (phase1_.divU().valid() && phase2_.divU().valid())
    {
        tdgdt =
        (
            alpha2.dimensionedInternalField()
           *phase1_.divU()().dimensionedInternalField()
          - alpha1.dimensionedInternalField()
           *phase2_.divU()().dimensionedInternalField()
        );
    }
    else if (phase1_.divU().valid())
    {
        tdgdt =
        (
            alpha2.dimensionedInternalField()
           *phase1_.divU()().dimensionedInternalField()
        );
    }
    else if (phase2_.divU().valid())
    {
        tdgdt =
        (
          - alpha1.dimensionedInternalField()
           *phase2_.divU()().dimensionedInternalField()
        );
    }

    alpha1.correctBoundaryConditions();
    surfaceScalarField alpha1f(fvc::interpolate(max(alpha1, scalar(0))));

    surfaceScalarField phic("phic", phi);
    surfaceScalarField phir("phir", phi1 - phi2);

    tmp<surfaceScalarField> alphaDbyA;

    if (notNull(phase1_.DbyA()) && notNull(phase2_.DbyA()))
    {
        surfaceScalarField DbyA(phase1_.DbyA() + phase2_.DbyA());

        alphaDbyA =
            fvc::interpolate(max(alpha1, scalar(0)))
           *fvc::interpolate(max(alpha2, scalar(0)))
           *DbyA;

        phir += DbyA*fvc::snGrad(alpha1, "bounded")*mesh_.magSf();
    }

    for (int acorr=0; acorr<nAlphaCorr; acorr++)
    {
        volScalarField::DimensionedInternalField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar("Sp", dimless/dimTime, 0.0)
        );

        volScalarField::DimensionedInternalField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            fvc::div(phi)*min(alpha1, scalar(1))
        );

        if (tdgdt.valid())
        {
            scalarField& dgdt = tdgdt();

            forAll(dgdt, celli)
            {
                if (dgdt[celli] > 0.0)
                {
                    Sp[celli] -= dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
                    Su[celli] += dgdt[celli]/max(1.0 - alpha1[celli], 1e-4);
                }
                else if (dgdt[celli] < 0.0)
                {
                    Sp[celli] += dgdt[celli]/max(alpha1[celli], 1e-4);
                }
            }
        }

        surfaceScalarField alphaPhic1
        (
            fvc::flux
            (
                phic,
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

        // Ensure that the flux at inflow BCs is preserved
        forAll(alphaPhic1.boundaryField(), patchi)
        {
            fvsPatchScalarField& alphaPhic1p =
                alphaPhic1.boundaryField()[patchi];

            if (!alphaPhic1p.coupled())
            {
                const scalarField& phi1p = phi1.boundaryField()[patchi];
                const scalarField& alpha1p = alpha1.boundaryField()[patchi];

                forAll(alphaPhic1p, facei)
                {
                    if (phi1p[facei] < 0)
                    {
                        alphaPhic1p[facei] = alpha1p[facei]*phi1p[facei];
                    }
                }
            }
        }

        if (nAlphaSubCycles > 1)
        {
            tmp<volScalarField> trSubDeltaT;

            if (LTS)
            {
                trSubDeltaT =
                    fv::localEulerDdt::localRSubDeltaT(mesh, nAlphaSubCycles);
            }

            for
            (
                subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
                !(++alphaSubCycle).end();
            )
            {
                surfaceScalarField alphaPhic10(alphaPhic1);

                MULES::explicitSolve
                (
                    geometricOneField(),
                    alpha1,
                    phi,
                    alphaPhic10,
                    (alphaSubCycle.index()*Sp)(),
                    (Su - (alphaSubCycle.index() - 1)*Sp*alpha1)(),
                    phase1_.alphaMax(),
                    0
                );

                if (alphaSubCycle.index() == 1)
                {
                    phase1_.alphaPhi() = alphaPhic10;
                }
                else
                {
                    phase1_.alphaPhi() += alphaPhic10;
                }
            }

            phase1_.alphaPhi() /= nAlphaSubCycles;
        }
        else
        {
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                alphaPhic1,
                Sp,
                Su,
                phase1_.alphaMax(),
                0
            );

            phase1_.alphaPhi() = alphaPhic1;
        }

        if (alphaDbyA.valid())
        {
            fvScalarMatrix alpha1Eqn
            (
                fvm::ddt(alpha1) - fvc::ddt(alpha1)
              - fvm::laplacian(alphaDbyA, alpha1, "bounded")
            );

            alpha1Eqn.relax();
            alpha1Eqn.solve();

            phase1_.alphaPhi() += alpha1Eqn.flux();
        }

        phase1_.alphaRhoPhi() =
            fvc::interpolate(phase1_.rho())*phase1_.alphaPhi();

        phase2_.alphaPhi() = phi - phase1_.alphaPhi();
        alpha2 = scalar(1) - alpha1;
        phase2_.alphaRhoPhi() =
            fvc::interpolate(phase2_.rho())*phase2_.alphaPhi();

        Info<< alpha1.name() << " volume fraction = "
            << alpha1.weightedAverage(mesh.V()).value()
            << "  Min(alpha1) = " << min(alpha1).value()
            << "  Max(alpha1) = " << max(alpha1).value()
            << endl;
    }
}


// ************************************************************************* //
