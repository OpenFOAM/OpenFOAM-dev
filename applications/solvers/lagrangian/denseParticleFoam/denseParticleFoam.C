/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2021 OpenFOAM Foundation
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

Application
    denseParticleFoam

Description
    Transient solver for the coupled transport of particle clouds including the
    effect of the volume fraction of particles on the continuous phase, with
    optional mesh motion and mesh topology changes.

\*---------------------------------------------------------------------------*/

#include "NamedEnum.H"

namespace Foam
{
    enum class cloudForceSplit
    {
        faceExplicitCellImplicit, // Implicit part of the cloud force added to
                                  // the cell momentum equation. Explicit part
                                  // to the face momentum equation. This is the
                                  // least likely to create staggering patterns
                                  // in the velocity field, but it can create
                                  // unphysical perturbations in cell
                                  // velocities even when particles and flow
                                  // have the similar velocities.

        faceExplicitCellLagged,   // Entire cloud force evaluated explicitly
                                  // and added to the face momentum equation.
                                  // Lagged correction (i.e.,
                                  // fvm::Sp(cloudSU.diag(), Uc) -
                                  // cloudSU.diag()*Uc) added to the cell
                                  // momentum equation. This creates physical
                                  // cell velocities when particles and flow
                                  // have the same velocity, but can also
                                  // result in staggering patterns in packed
                                  // beds. Unsuitable for MPPIC.

        faceImplicit              // Implicit and explicit parts of the force
                                  // both added to the face momentum equation.
                                  // Behaves somewhere between the other two.
    };

    template<>
    const char* NamedEnum<cloudForceSplit, 3>::names[] =
    {
        "faceExplicitCellImplicit",
        "faceExplicitCellLagged",
        "faceImplicit"
    };

    const NamedEnum<cloudForceSplit, 3> cloudForceSplitNames;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "phaseKinematicMomentumTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "parcelCloudList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUcfIfPresent.H"
    #include "initContinuityErrs.H"

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Store the particle positions
        clouds.storeGlobalPositions();

        mesh.update();

        if (mesh.changing())
        {
            if (correctPhi)
            {
                #include "correctPhic.H"
            }

            if (checkMeshCourantNo)
            {
                #include "meshCourantNo.H"
            }
        }

        continuousPhaseTransport.correct();
        muc = rhoc*continuousPhaseTransport.nu();

        clouds.evolve();

        // Update continuous phase volume fraction field
        alphac = max(1.0 - clouds.theta(), alphacMin);
        alphac.correctBoundaryConditions();
        alphacf = fvc::interpolate(alphac);
        alphaPhic = alphacf*phic;

        // Cloud forces
        fvVectorMatrix cloudSU(clouds.SU(Uc));
        volVectorField cloudSUu
        (
            IOobject
            (
                "cloudSUu",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedVector(dimForce/dimVolume, Zero),
            zeroGradientFvPatchVectorField::typeName
        );
        volScalarField cloudSUp
        (
            IOobject
            (
                "cloudSUp",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar(dimForce/dimVelocity/dimVolume, Zero),
            zeroGradientFvPatchVectorField::typeName
        );

        const cloudForceSplit cloudSUSplit =
            pimple.dict().found("cloudForceSplit")
          ? cloudForceSplitNames.read(pimple.dict().lookup("cloudForceSplit"))
          : cloudForceSplit::faceExplicitCellImplicit;

        switch (cloudSUSplit)
        {
            case cloudForceSplit::faceExplicitCellImplicit:
                cloudSUu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
                cloudSUu.correctBoundaryConditions();
                cloudSUp.primitiveFieldRef() = Zero;
                cloudSUp.correctBoundaryConditions();
                //cloudSU.diag() = cloudSU.diag();
                cloudSU.source() = Zero;
                break;

            case cloudForceSplit::faceExplicitCellLagged:
                cloudSUu.primitiveFieldRef() =
                    (cloudSU.diag()*Uc() - cloudSU.source())/mesh.V();
                cloudSUu.correctBoundaryConditions();
                cloudSUp.primitiveFieldRef() = Zero;
                cloudSUp.correctBoundaryConditions();
                //cloudSU.diag() = cloudSU.diag();
                cloudSU.source() = cloudSU.diag()*Uc();
                break;

            case cloudForceSplit::faceImplicit:
                cloudSUu.primitiveFieldRef() = -cloudSU.source()/mesh.V();
                cloudSUu.correctBoundaryConditions();
                cloudSUp.primitiveFieldRef() = cloudSU.diag()/mesh.V();
                cloudSUp.correctBoundaryConditions();
                cloudSU.diag() = Zero;
                cloudSU.source() = Zero;
                break;
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fvModels.correct();

            #include "UcEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                continuousPhaseTurbulence->correct();
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
